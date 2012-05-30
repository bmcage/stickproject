#
# Copyright (C) 2010  B. Malengier
# Copyright (C) 2010  P.Li
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

""" 
Module holding a generic diffusion model for a yarn. 
"""
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
from __future__ import division
import os.path
import sys
import const
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

HAVE_ODES = False
try:
    from scikits.odes import ode as sc_ode
    HAVE_ODES = True
except:
    print 'Could not load scikits.odes, odes solver not available'
    sys.exit(0)

from scikits.odes.sundials.common_defs import CV_RootFunction

CV_SUCCESS = 0
CV_TSTOP_RETURN = 1
CV_ROOT_RETURN = 2
#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils

#-------------------------------------------------------------------------
#
# PCMState class 
#
#-------------------------------------------------------------------------
class PCMState(object):
    """
    A PCM can be in different states: solid, liquid, in-between. Depending
    on the state characteristics change, and the solution method to model 
    the PCM is different
    """
    #different states we consider, first state at r=0, then at r=L
    INVALID = -2
    UNSET = -1
    SOLID = 0
    LIQUID = 1
    SOLID_LIQUID = 2
    LIQUID_SOLID = 3
    
    def __init__(self, config):
        self.state = PCMState.UNSET
        self.meltpoint = config.get("pcm.melting_point")
        self.L = config.get("pcm.radius")  # in mm !
        self.epsilon = 1e-4
        self.latent_heat_fusion = 1000 * config.get("pcm.latent_heat_fusion")
        #units in mm instead of m
        self.solid = {
            'state' : PCMState.SOLID,
            'rho' : config.get("pcm.density") * 10**(-9),
            'C'   : 1000 * config.get("pcm.specific_heat_solid"),
            'K'   : config.get("pcm.thermal_cond_solid") * 10**(-3)
            }
        self.liquid = {
            'state' : PCMState.LIQUID,
            'rho' : config.get("pcm.density") * 10**(-9),
            'C'   : 1000 * config.get("pcm.specific_heat_liquid"),
            'K'   : config.get("pcm.thermal_cond_liquid") * 10**(-3)
            }

        #we store computational grid in the state, as we need to swap grids
        #when state changes
        self.n_edge = config.get('discretization.n_edge') #discretize the PCM radius
        self.grid = None
        self.grid_edge = None
        self.outer_data = None
        self.inner_data = None
    
    def set_state(self, init_cond):
        """
        init_cond should be a function over (0, L) giving initial temperature
        """
        grid = np.linspace(0., self.L, 100)
        state = PCMState.INVALID
        if init_cond(grid[0]) < self.meltpoint:
            state = PCMState.SOLID
        elif init_cond(grid[0]) > self.meltpoint:
            state = PCMState.LIQUID
        else:
            self.state = state
            return

        self.R = 0.
        for rpos in grid[1:]:
            if state == PCMState.SOLID and init_cond(rpos) <= self.meltpoint:
                continue
            elif state == PCMState.LIQUID and init_cond(rpos) >= self.meltpoint:
                continue
            elif state == PCMState.LIQUID and init_cond(rpos) < self.meltpoint:
                self.R = rpos
                state = PCMState.LIQUID_SOLID
            elif state == PCMState.SOLID and init_cond(rpos) > self.meltpoint:
                self.R = rpos
                state = PCMState.SOLID_LIQUID
            elif state == PCMState.SOLID_LIQUID and init_cond(rpos) > self.meltpoint:
                continue
            elif state == PCMState.LIQUID_SOLID and init_cond(rpos) < self.meltpoint:
                continue
            elif state == PCMState.LIQUID_SOLID and init_cond(rpos) > self.meltpoint:
                # this is not supported in the model, two interfaces
                self.state = PCMState.INVALID
                break
            elif state == PCMState.SOLID_LIQUID and init_cond(rpos) < self.meltpoint:
                # this is not supported in the model, two interfaces
                self.state = PCMState.INVALID
                break
            else:
                print 'unexpected data, cannot determine state PCM, stopping'
                sys.exit(0)
        self.state = state
        #we set inner data, and outer data
        if self.state == PCMState.SOLID:
            self.outer_data = self.solid
            self.inner_data = self.solid
        elif self.state == PCMState.LIQUID:
            self.outer_data = self.liquid
            self.inner_data = self.liquid
        elif self.state == PCMState.LIQUID_SOLID:
            self.outer_data = self.solid
            self.inner_data = self.liquid
        elif self.state == PCMState.SOLID_LIQUID:
            self.outer_data = self.liquid
            self.inner_data = self.solid
        else:
            self.outer_data = None
            self.inner_data = None

    def single_phase(self):
        if self.state in [PCMState.SOLID, PCMState.LIQUID]:
            return True
        return False

    def calc_grid(self):
        """
        Determines the fixed grid on the unit inteval for the current state
        """
        self.inner_gridx_edge = sp.linspace(0., 1., self.n_edge)
        self.outer_gridx_edge = sp.linspace(0., 1., self.n_edge)
        #construct cell centers from this
        self.inner_gridx = (self.inner_gridx_edge[:-1] + self.inner_gridx_edge[1:])/2.
        #obtain cell sizes
        self.inner_delta_x = self.inner_gridx_edge[1:] - self.inner_gridx_edge[:-1]
        self.inner_delta_x_avg = (self.inner_delta_x[:-1] + self.inner_delta_x[1:])/2
        #construct cell centers from this
        self.outer_gridx = (self.outer_gridx_edge[:-1] + self.outer_gridx_edge[1:])/2.
        #obtain cell sizes
        self.outer_delta_x = self.outer_gridx_edge[1:] - self.outer_gridx_edge[:-1]
        self.outer_delta_x_avg = (self.outer_delta_x[:-1] + self.outer_delta_x[1:])/2

    def outer_x_to_r(self, x, R=None):
        #conversion function
        if R == None:
            R = self.R
        return self.L - (self.L - R) * x

    def inner_x_to_r(self, x, R=None):
        #conversion function
        if R == None:
            R = self.R
        return R * x 

    def reset_state(self, temp_outer, Rintf):
        """
        Reset the state based on new outer temperature and new interface
        position Rintf. New interface position will be stored in self.R (which
        will be Rintf or 0)!
        There are the following possibilities
        1. single phase, en temp_outer is at edges at melting point 
            ==> becomes two phases
        2. two phases, but R becomes close to edge 
            ==> becomes one phase
            
        return value: tuple (changed, switch)
            changed = True if state changed
            switch = True if inner and outer must be switched
        """
        self.R = Rintf
        changed = False
        switch = False
        if self.state == PCMState.SOLID:
            if temp_outer[0] >= self.meltpoint:
                self.state = PCMState.SOLID_LIQUID
                self.R = (1. - self.epsilon) * self.L
                changed = True
                switch = True
        elif self.state == PCMState.LIQUID:
            if temp_outer[0] <= self.meltpoint:
                self.state = PCMState.LIQUID_SOLID
                self.R = (1. - self.epsilon) * self.L
                changed = True
                switch = True
        elif self.state == PCMState.SOLID_LIQUID:
            if abs(self.R-self.L) < self.epsilon *self.L:
                self.state = PCMState.SOLID
                self.R = 0.
                changed = True
                switch = True
            elif self.R < self.epsilon * self.L:
                self.R = 0.
                self.state = PCMState.LIQUID
                changed = True
                switch = False
        elif self.state == PCMState.LIQUID_SOLID:
            if abs(self.R-self.L) < self.epsilon *self.L:
                self.R = 0.
                self.state = PCMState.LIQUID
                changed = True
                switch = True
            elif self.R < self.epsilon * self.L:
                self.R = 0.
                self.state = PCMState.SOLID
                changed = True
                switch = False

        #we set inner data, and outer data
        if self.state == PCMState.SOLID:
            self.outer_data = self.solid
            self.inner_data = self.solid
        elif self.state == PCMState.LIQUID:
            self.outer_data = self.liquid
            self.inner_data = self.liquid
        elif self.state == PCMState.LIQUID_SOLID:
            self.outer_data = self.solid
            self.inner_data = self.liquid
        elif self.state == PCMState.SOLID_LIQUID:
            self.outer_data = self.liquid
            self.inner_data = self.solid
        else:
            self.outer_data = None
            self.inner_data = None

        return (changed, switch)

class RootFnsp(CV_RootFunction):
    '''single phase rootfunction to identify when it becomes double phase'''
    def set_data(self, xend, meltpoint):
        self.L = xend
        self.meltpoint = meltpoint

    def evaluate(self, t, u, out, userdata):
        out[0] = u[0]/self.L - self.meltpoint
        #print 'test root', t, u[0]/self.L , self.meltpoint, out[0]
        return 0

class RootFndp(CV_RootFunction):
    '''double phase rootfunction to identify when it becomes single phase'''
    def set_data(self, L, pos_s, eps=1e-4):
        self.L = L
        self.pos_s = pos_s
        self.eps = eps
        self.minval = eps/2.
        self.maxval = L-eps/2.

    def evaluate(self, t, u, out, userdata):
        out[0] = u[self.pos_s] - self.minval
        out[1] = self.maxval - u[self.pos_s]
        print 'test root', u[self.pos_s], out[0], out[1]
        return 0

#-------------------------------------------------------------------------
#
# PCMModel class 
#
#-------------------------------------------------------------------------
class PCMModel(object):
    """
    pcm.PCMModel is a special diffusion model for a single spherical
    PCM which is composed of a specific material that can undergo melting.
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        
        """
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')

        self.solver = None
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        #set correct delta_t
        self.delta_t = self.times[1]-self.times[0]
        if self.verbose:
            print "Timestep used in pcm model:", self.delta_t

        self.state = PCMState(self.cfg)
        #unit in mm instead of m
        self.hT = self.cfg.get('boundary.heat_transfer_coeff') * 10**(-6)
        self.lambda_m = self.cfg.get('pcm.latent_heat_fusion') * 10**3
        self.Tout = eval(self.cfg.get('boundary.T_out'))

        self.initialized = False

        self.plotevery = self.cfg.get("plot.plotevery")

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh (in mm) on which we project results
        grid: the space position of each central point for every cell element;
        """
        self.init_temp = eval(self.cfg.get('init.init_temp'))
        self.state.set_state(self.init_temp)
        self.state.calc_grid()

    def initial_PCM(self):
        """
        Initialize PCM data
        """
        self.initial_T_in = sp.empty(self.state.n_edge-1, float)
        self.initial_T_out = sp.empty(self.state.n_edge-1, float)
        # we have a state, now we construct an initial condition over the grid
        for i, pos in enumerate(self.state.outer_gridx):
            self.initial_T_out[i] = self.init_temp(self.state.outer_x_to_r(pos))
        for i, pos in enumerate(self.state.inner_gridx):
            self.initial_T_in[i] = self.init_temp(self.state.inner_x_to_r(pos))

        self.volume = self.calc_volume()
        if self.verbose:
            temps = [0., 10, 20, 25, 30, 35, 40]
            eno = self.calc_energy(self.initial_T_in, 
                            self.initial_T_out, temps)
            print 'initial energy = ', eno[0], 'J from base 0 degree Celcius'
            olden = 0.;
            for tt, en in zip(temps, eno[1]):
                print '  for temp ', tt, 'C, energy =', en , 'diff en =',en-olden
                olden = en

        self.unknowns_single = sp.empty(self.state.n_edge-1, float)
        self.unknowns_double = sp.empty(2*(self.state.n_edge-1)+1, float)
        if self.state.single_phase():
            #only outer
            self.unknowns_single[:] = self.initial_T_out[:]*self.state.outer_x_to_r(self.state.outer_gridx[:])
            self.unknowns = self.unknowns_single
        else:
            self.unknowns_double[0:self.state.n_edge-1] = \
                        self.initial_T_out[:]*self.state.outer_x_to_r(self.state.outer_gridx[:])
            self.unknowns_double[self.state.n_edge-1] = self.state.R
            self.unknowns_double[self.state.n_edge:] = \
                        self.initial_T_in[::-1]*self.state.inner_x_to_r(self.state.inner_gridx[::-1])
            self.unknowns = self.unknowns_double

    def calc_volume(self):
        """ volume in mm^3 """
        return 4/3 * np.pi * self.state.L**3 

    def calc_energy(self, inner_T, outer_T, base=None):
        """ calculate the energy in the PCM based on temperature over the 
            inner x grid, and outer x grid 
            In J, with 0 degree Celcius equal to 0 J! """
        Cl = self.state.liquid['C'] # inner specific heat J/kg K
        rhol = self.state.liquid['rho']
        Cs = self.state.solid['C'] # inner specific heat J/kg K
        rhos = self.state.solid['rho']
        meltpoint = self.state.meltpoint
        innerE = 0.
        prevedgex = 0.
        for xpos, temp in zip(self.state.inner_gridx_edge[1:], self.initial_T_in):
            volshell = 4/3*np.pi * \
                (self.state.inner_x_to_r(xpos)**3 
                 - self.state.inner_x_to_r(prevedgex)**3)
            prevedgex = xpos
            if temp > meltpoint:
                innerE += Cl * (temp-meltpoint) * rhol * volshell \
                            + Cs * meltpoint * rhos * volshell\
                            + self.lambda_m * rhos * volshell
            else:
                innerE += Cs * temp * rhos * volshell
        outerE = 0.
        prevedgex = 0.
        for xpos, temp in zip(self.state.outer_gridx_edge[1:], self.initial_T_out):
            volshell = 4/3*np.pi * \
                (self.state.outer_x_to_r(prevedgex)**3
                 - self.state.outer_x_to_r(xpos)**3 )
            prevedgex = xpos
            if temp > meltpoint:
                outerE += Cl * (temp-meltpoint) * rhol * volshell \
                            + Cs * meltpoint * rhos * volshell \
                            + self.lambda_m * rhos * volshell
            else:
                outerE += Cs * temp * rhos * volshell
        #print 'test in - out energy', innerE, outerE
        out = []
        for temp in base or []:
            vol = self.calc_volume()
            if temp > meltpoint:
                heatmelt = self.lambda_m * rhos * vol
                out += [Cl * (temp - meltpoint) * rhol * vol 
                        + Cs * meltpoint* rhos * vol + heatmelt]
            else:
                out += [Cs * temp * rhos * vol]
        if base:
            return innerE + outerE, out
        else:
            return innerE + outerE

    def f_odes_sph(self, t, u_rep, diff_u_t):
        """ RHS function for the cvode solver for single phase energy diffusion
        The grid is in x from 0 to 1, with only outer data, so 0 is at L and
        1 is at center of the sample
        """
        grid = self.state.outer_gridx
        grid_edge = self.state.outer_gridx_edge
        Dx = self.state.outer_delta_x
        Dxavg = self.state.outer_delta_x_avg
        n_cell = len(grid)
        flux_edge = self.__tmp_flux_edge[:n_cell+1]
        L = self.state.L
        Cv = self.state.outer_data['C']
        K = self.state.outer_data['K']
        rho = self.state.outer_data['rho']
        dxdr =  -1./L

        #print difference in temperature! 
        #temp = u_rep[:]/self.state.outer_x_to_r(self.state.outer_gridx)
        #print 'temp', temp
        Rlast = self.state.outer_x_to_r(self.state.outer_gridx[0])

        flux_edge[-1] = -K/rho/Cv * dxdr * u_rep[-1] \
                        /self.state.outer_x_to_r(self.state.outer_gridx[-1])
        flux_edge[0] = (self.hT / rho / Cv * dxdr*L*
                         (u_rep[0]/Rlast - self.Tout(t))
                        -K/rho/Cv * dxdr*u_rep[0] /Rlast)
        flux_edge[1:-1] = -K/rho/Cv * (u_rep[1:]-u_rep[:-1])/ Dxavg[:]*dxdr**2
        diff_u_t[:] = (flux_edge[:-1]-flux_edge[1:])/Dx[:]
        #print 'flux', t, flux_edge, diff_u_t[0:2]
        #raw_input()

    def f_odes_dph(self, t, u_rep, diff_u_t):
        """ RHS function for the cvode solver for double phase energy diffusion
        The grid is in x with first x from 0 to 1 with outer data, so 0 is 
        at L and 1 is at interface, then the interface position, then x from 1
        to 0 with innter data, with 1 the interface, and 0 the center
        """
        pass

    def solve_odes_reinit(self):
        """
        Reinitialize the cvode solver to start again
        """
        self.initial_t = self.times[0]

        if self.state.single_phase():
            self.solverunknowns = self.pos_s
            rootfn=RootFnsp()
            rootfn.set_data(self.state.outer_x_to_r(self.state.outer_gridx[0]), 
                            self.state.meltpoint)
            self.solver = sc_ode('cvode', self.f_odes_sph,
                                 max_steps=50000, lband=1, uband=1,
                                 nr_rootfns=1, rootfn=rootfn)
        else:
            self.solverunknowns = self.nrunknowns
            rootfn=RootFndp()
            rootfn.set_data(self.state.L, self.pos_s, self.state.epsilon)
            self.solver = sc_ode('cvode', self.f_odes_dph,
                                 max_steps=50000, lband=1, uband=1,
                                 nr_rootfns=2, rootfn=rootfn)
        self.solver.init_step(self.step_old_time, self.step_old_sol)

    def solve_init(self):
        """
        Initialize the solver so they can be solved stepwize
        """
        self.nrunknowns = 2*(self.state.n_edge-1)+1
        if not HAVE_ODES:
            raise Exception, 'Not possible to solve with given method, scikits.odes not available'
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        self.step_old_sol = self.unknowns
        
        self.pos_s = self.state.n_edge-1
        #data storage
        self.all_sol_u = np.empty((len(self.times), self.nrunknowns), float)
        self.ret_sol = sp.empty(self.nrunknowns, float)

        if self.state.single_phase():
            self.all_sol_u[0][self.pos_s] = 0.
            self.all_sol_u[0][:self.pos_s] = self.unknowns[:]
        else:
            self.all_sol_u[0][:] = self.unknowns[:]
        self.__tmp_diff_sol_t = sp.empty(self.nrunknowns, float)
        self.__tmp_flux_edge = sp.empty(self.nrunknowns+1, float)

        self.tstep = 0
        self.solve_odes_reinit()
        self.initialized = True

    def solve(self):
        """
        Solve the PCM model
        """
        if not self.initialized:
            print 'ERROR, solver not initialized'
            sys.exit()

        tstep = 0
        single_phase = self.state.single_phase()
        changed = False
        cont = True
        while cont:
            if single_phase:
                time = self.times[tstep+1]
                if changed:
                    prev_time = self.step_old_time
                    #retrieve stored start point
                    initval = self.step_old_sol
                    changed = False
                else:
                    prev_time = self.times[tstep]
                    initval = self.all_sol_u[tstep][:self.solverunknowns]
                self.solver.set_tcrit(time)
                flag, t_retn, u_retn, t_out, u_last = self.solver.solve(
                                [prev_time,time], 
                                initval,
                                )
                #print flag, t_retn, time, t_out
                if flag == CV_ROOT_RETURN:
                    print 'At time', t_out, 'no longer single phase'
                    #store the last solution
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_last
                    self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_out
                    single_phase = False
                    changed = True
                elif flag == CV_TSTOP_RETURN:
                    #here crit time is just the last good time we want
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_last
                    self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_out
                elif flag == CV_SUCCESS:
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_retn[-1]
                    self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_retn[-1]
                else:
                    print 'ERROR: unable to compute solution, flag', flag
                    break
                    
                if self.verbose:
                    print 'INFO: pcmmodel at t = ', tret
                if single_phase:
                    tstep += 1
            else:
                #double phase model to run
                time = self.times[tstep+1]
                if changed:
                    prev_time = self.step_old_time
                    #retrieve stored start point
                    initval = self.step_old_sol
                    changed = False
                else:
                    prev_time = self.times[tstep]
                    initval = self.all_sol_u[tstep][:self.solverunknowns]
                self.solver.set_tcrit(time)
                flag, t_retn, u_retn, t_out, u_last = self.solver.solve(
                                [prev_time,time], 
                                initval,
                                )
                
                if flag == CV_ROOT_RETURN:
                    print 'At time', t_out, 'no longer double phase'
                    #store the last solution
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_last
                    self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_out
                    single_phase = True
                    changed = True
                elif flag == CV_TSTOP_RETURN:
                    #here crit time is just the last good time we want
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_last
                    self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_out
                elif flag == CV_SUCCESS:
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_retn[-1]
                    self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_retn[-1]
                else:
                    print 'ERROR: unable to compute solution, flag', flag
                    break
                    
                if self.verbose:
                    print 'INFO: pcmmodel at t = ', tret
                if not single_phase:
                    tstep += 1

            if tstep == len(self.times)-1:
                break
            if changed:
                # we need to swap the solver we use!
                self.step_old_time = tret
                if single_phase:
                    #we went from double phase to single phase
                    #1. determine what phase
                    #2. change solver
                    #3. update data
                    if abs(self.all_sol_u[tstep+1][self.pos_s]) < self.state.epsilon:
                        # now only outer state
                        self.state.state = self.state.outer_data['state']
                        self.state.inner_data = self.state.outer_data
                        self.step_old_sol = self.project_outdp_sp(self.all_sol_u[tstep+1])
                    elif abs(self.state.L - 
                             self.all_sol_u[tstep+1][self.pos_s]) \
                            < self.state.epsilon:
                        # now only inner state
                        self.state.state = self.state.inner_data['state']
                        self.state.outer_data = self.state.inner_data
                        self.step_old_sol = self.project_indp_sp(self.all_sol_u[tstep+1])
                    else:
                        raise NotImplementedError, 'We should not reach this'
                    # no more interface:
                    self.state.R = 0.
                else:
                    #we went from single phase to double phase, interface
                    #arises at R=L-eps/2
                    #1. determine what phase
                    #2. change solver
                    #3. update data
                    self.state.R = self.state.L - self.state.epsilon
                    if self.state.outer_data['state'] == PCMState.LIQUID:
                        self.state.state = PCMState.LIQUID_SOLID
                        self.state.inner_data = self.state.outer_data
                        self.state.outer_data = self.state.solid
                    elif self.state.outer_data['state'] == PCMState.SOLID:
                        self.state.state = PCMState.SOLID_LIQUID
                        self.state.inner_data = self.state.outer_data
                        self.state.outer_data = self.state.liquid
                    else:
                        raise NotImplementedError, 'We should not reach this'
                    self.step_old_sol = self.project_outsp_dp(self.all_sol_u[tstep+1],
                                            self.state.R)
                    
                self.solve_odes_reinit()
                    
        self.last_sol_tstep = tstep
        
        self.view_sol(self.times, self.all_sol_u)

    def project_outsp_dp(self, sol, R):
        """ project a single phase solution (so only outer) to a double phase
            solution, with interface at R """
        print 'introducing interface'
        #The outer solution is set at meltingtemp
        newdataC = np.empty(len(sol), float)
        newrout = self.state.outer_x_to_r(self.state.outer_gridx, R=R)
        for i in xrange(self.pos_s):
            newdataC[:self.pos_s] = self.state.meltpoint * newrout
        #new interface is set
        newdataC[self.pos_s] = R
        #and we project the old outer solution to inner solution, which we 
        # store inverted
        # SIMPLE: assume projection can be neglected
        newdataC[self.pos_s+1:] = sol[:self.pos_s]
        return newdataC

    def project_outdp_sp(self, sol):
        """ project a double phase solution sol, to a single phase solution
            over the outer grid, assuming interface goes to 0 """
        print 'interface at', sol[self.pos_s], ' projecting'
        origr = self.state.outer_x_to_r(self.state.outer_gridx, R=sol[self.pos_s])
        origredge = self.state.outer_x_to_r(self.state.outer_gridx_edge, R=sol[self.pos_s])
        # conserved is volumeshell T = volumeshell U/r
        origdataC = sol[:self.pos_s] / origr * \
                        4/3*np.pi * (origredge[:-1]**3 - origredge[1:]**3)
        newr = self.state.outer_x_to_r(self.state.outer_gridx, R=0.)
        newredge = self.state.outer_x_to_r(self.state.outer_gridx_edge, R=0.)
        newdataC = np.empty(self.pos_s, float)
        pos = 0
        posorig=0
        #remaining piece goes to meltpoint temperature
        newC = self.state.meltpoint * 4/3*np.pi * sol[self.pos_s]**3
        for pos in np.arange(self.pos_s):
            while origredge[posorig] > newredge[pos+1]:
                if origredge[posorig+1] <= newredge[pos+1]:
                    nextC = -sol[posorig]/origr[posorig] * \
                        4/3*np.pi * (origredge[posorig+1]**3 - newredge[pos+1]**3)
                    newC += -sol[posorig]/origr[posorig] * \
                        4/3*np.pi * (newredge[pos+1]**3 - origredge[posorig]**3)
                else:
                    nextC = 0.
                    newC += -sol[posorig]/origr[posorig] * \
                        4/3*np.pi * (origredge[posorig+1]**3 - origredge[posorig]**3)
                posorig += 1
            newdataC[pos] = newC
            newC = nextC
            nextC = 0.
        #now we derive the new u = r T_avg
        newdataC = newdataC * newr \
                    / (4/3*np.pi * (newredge[:-1]**3 - newredge[1:]**3))
        return newdataC


    def project_indp_sp(self, sol):
        """ project a double phase solution sol, to a single phase solution
            over the outer grid, assuming interface goes to L, so inner part
            of double phase needs to be considered
        """
        print 'interface at', sol[self.pos_s], ' projecting'
        origr = self.state.inner_x_to_r(self.state.inner_gridx, R=sol[self.pos_s])
        origredge = self.state.inner_x_to_r(self.state.inner_gridx_edge, R=sol[self.pos_s])
        # conserved is volumeshell T = volumeshell U/r
        origdataC = sol[self.pos_s+1:][::-1] / origr * \
                        4/3*np.pi * (origredge[:-1]**3 - origredge[1:]**3)
        newr = self.state.outer_x_to_r(self.state.outer_gridx, R=0.)
        newredge = self.state.outer_x_to_r(self.state.outer_gridx_edge, R=0.)
        newdataC = np.empty(self.pos_s, float)
        pos = 0
        posorig=self.pos_s-1
        #remaining piece goes to meltpoint temperature
        newC = self.state.meltpoint * 4/3*np.pi * (self.state.L**3-sol[self.pos_s]**3)
        
        for pos in np.arange(self.pos_s):
            while origredge[posorig] > newredge[pos+1]:
                if origredge[posorig-1] <= newredge[pos+1]:
                    nextC = -sol[posorig]/origr[posorig] * \
                        4/3*np.pi * (origredge[posorig-1]**3 - newredge[pos+1]**3)
                    newC += -sol[posorig]/origr[posorig] * \
                        4/3*np.pi * (newredge[pos+1]**3 - origredge[posorig]**3)
                else:
                    nextC = 0.
                    newC += -sol[posorig]/origr[posorig] * \
                        4/3*np.pi * (origredge[posorig-1]**3 - origredge[posorig]**3)
                posorig -= 1
            newdataC[pos] = newC
            newC = nextC
            nextC = 0.
        #now we derive the new u = r T_avg
        newdataC = newdataC * newr \
                    / (4/3*np.pi * (newredge[:-1]**3 - newredge[1:]**3))
        return newdataC


    def view_sol(self, times, rTemp):
        """
        Show the solution in Temp with times.
        rTemp[i][:] contains u = r*T at time times[i]
        """
        from fipy import CellVariable, Matplotlib1DViewer, Grid1D
        if rTemp[0][self.pos_s] == 0.:
            meshr = self.state.outer_x_to_r(self.state.outer_gridx, 
                                                        rTemp[0][self.pos_s])
            dr = -(self.state.outer_x_to_r(
                self.state.outer_delta_x, rTemp[0][self.pos_s])-self.state.L)
            mesh_PCM = Grid1D(dx=tuple(dr[::-1]))
            value = rTemp[0][:self.pos_s] / meshr[:]
            solution_view = CellVariable(name = "PCM temperature", 
                                        mesh = mesh_PCM, value = value[::-1])
        else:
            raise NotImplementedError
        if self.plotevery:
            self.viewer =  Matplotlib1DViewer(vars = solution_view, datamin=0., 
                                datamax=45.)
            self.viewer.plot()
        self.viewerplotcount = 0
        for time, rT in zip(times[1:self.last_sol_tstep], rTemp[1:self.last_sol_tstep][:]):
            if rT[self.pos_s] == 0.:
                #meshr = self.state.outer_x_to_r(self.state.outer_gridx, rT[self.pos_s])
                #print' sol', rT[:self.pos_s]/meshr
                value = rT[:self.pos_s]/meshr
                if self.plotevery and self.viewerplotcount == 0:
                    solution_view.setValue(value[::-1])
                    self.viewer.axes.set_title('PCM Temp vs radius at time %s' %str(time))
                    self.viewer.axes.set_xlabel('Radius')
                    self.viewer.axes.set_ylabel('Temp')
                    self.viewer.plot()
                    #self.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'PCM%sconc%08.4f.png' % (name,time))
                #else:
                #    self.viewer.plot()
                if self.plotevery:
                    self.viewerplotcount += 1
                    self.viewerplotcount = self.viewerplotcount % self.plotevery
            else:
                raise NotImplementedError  

    def dump_solution(self):
        pass

    def run_init(self):
        self.create_mesh()
        self.initial_PCM()

    def run(self, wait=False, output=False):
        self.run_init()
        if not self.initialized:
            self.solve_init()
        self.solve()

        if output:
            self.dump_solution()
        if wait:
            raw_input("Finished PCM run")

    def __del__(self):
        #remove the memory
        if self.solver:
            print 'del self.solver'
            del self.solver
        self.solver = None
