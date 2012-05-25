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
        self.latent_heat_fusion = config.get("pcm.latent_heat_fusion")
        self.solid = {
            'rho' : config.get("pcm.density"),
            'C'   : 1000 * config.get("pcm.specific_heat_solid"),
            'K'   : config.get("pcm.thermal_cond_solid")
            }
        self.liquid = {
            'rho' : config.get("pcm.density"),
            'C'   : 1000 * config.get("pcm.specific_heat_liquid"),
            'K'   : config.get("pcm.thermal_cond_liquid")
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
        print 'test root', t, u[0]/self.L , self.meltpoint, out[0]
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
        self.hT = self.cfg.get('boundary.heat_transfer_coeff')
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
            print 'initial energy = ', self.calc_energy(self.initial_T_in, 
                            self.initial_T_out), 'J from base 0 degree Celcius'

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
        """ volume in m^3 """
        return 4/3 * np.pi * self.state.L**3 * 10**(-9)

    def calc_energy(self, inner_T, outer_T):
        """ calculate the energy in the PCM based on temperature over the 
            inner x grid, and outer x grid 
            In J, with 0 degree Celcius equal to 0 J! """
        Ci = self.state.inner_data['C'] * 1000 # inner specific heat J/kg K
        rhoi = self.state.inner_data['rho']
        Co = self.state.outer_data['C'] * 1000 # inner specific heat J/kg K
        rhoo = self.state.outer_data['rho']
        innerE = 0.
        prevedgex = 0.
        for xpos, temp in zip(self.state.inner_gridx_edge[1:], self.initial_T_in):
            volshell = 4/3*np.pi *10**(-9) * \
                (self.state.inner_x_to_r(xpos)**3 
                 - self.state.inner_x_to_r(prevedgex)**3)
            prevedgex = xpos
            innerE += Ci * temp * rhoi * volshell
        outerE = 0.
        prevedgex = 0.
        for xpos, temp in zip(self.state.outer_gridx_edge[1:], self.initial_T_out):
            volshell = 4/3*np.pi *10**(-9) * \
                (self.state.outer_x_to_r(prevedgex)**3
                 - self.state.outer_x_to_r(xpos)**3 )
            prevedgex = xpos
            outerE += Ci * temp * rhoi * volshell
        #print 'test in - out energy', innerE, outerE
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
        temp = u_rep[:]/self.state.outer_x_to_r(self.state.outer_gridx)
        print 'temp', temp

        # TODO TODO why L*L ??
        flux_edge[-1] = -L*L*K/rho/Cv * dxdr * u_rep[-1] \
                        /self.state.outer_x_to_r(self.state.outer_gridx[-1])**2 #self.state.outer_x_to_r(self.state.outer_gridx[-1])
        flux_edge[0] = (self.hT / rho / Cv * dxdr*L*(u_rep[0]/L - self.Tout(t))
                        -K/rho/Cv * dxdr*u_rep[0] /L)
        flux_edge[1:-1] = -K/rho/Cv * (u_rep[1:]-u_rep[:-1])/ Dxavg[:]*dxdr*dxdr
        print 'flux', flux_edge
        #print 'test flux', u_rep[-2], u_rep[-1], flux_edge[-2], flux_edge[-1]
        diff_u_t[:] = -(flux_edge[:-1]-flux_edge[1:])/Dx[:]
        raw_input()

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
            rootfn=RootFnsp()
            rootfn.set_data(self.state.outer_x_to_r(self.state.outer_gridx[0]), 
                            self.state.meltpoint)
            self.solver = sc_ode('cvode', self.f_odes_sph,
                                 max_steps=50000, lband=1, uband=1,
                                 nr_rootfns=1, rootfn=rootfn)
        else:
            rootfn=RootFndp()
            rootfn.set_data(self.state.L, self.pos_s, self.epsilon)
            self.solver = sc_ode('cvode', self.f_odes_dph,
                                 max_steps=50000, lband=1, uband=1,
                                 nr_rootfns=2, rootfn=rootfn)
        self.solver.init_step(self.step_old_time, self.step_old_sol)

    def solve_init(self):
        """
        Initialize the solver so they can be solved stepwize
        """
        nrunknowns = 2*(self.state.n_edge-1)+1
        if not HAVE_ODES:
            raise Exception, 'Not possible to solve with given method, scikits.odes not available'
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        self.step_old_sol = self.unknowns
        
        self.pos_s = self.state.n_edge-1
        #data storage
        self.all_sol_u = np.empty((len(self.times), nrunknowns), float)
        self.ret_sol = sp.empty(nrunknowns, float)

        if self.state.single_phase():
            self.solverunknowns = self.pos_s
            self.all_sol_u[0][self.pos_s] = 0.
            self.all_sol_u[0][:self.pos_s] = self.unknowns[:]
        else:
            self.solverunknowns = nrunknowns
            self.all_sol_u[0][:] = self.unknowns[:]
        self.__tmp_diff_sol_t = sp.empty(nrunknowns, float)
        self.__tmp_flux_edge = sp.empty(nrunknowns+1, float)

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
        for time in self.times[1:]:
            if single_phase:
                self.solver.set_tcrit(time)
                flag, t_retn, u_retn, t_out, u_last = self.solver.solve(
                                    [self.times[tstep],time], 
                                    self.all_sol_u[tstep][:self.solverunknowns],
                                    )
                print flag, t_retn, time
                if flag == CV_ROOT_RETURN:
                    print 'At time', t_out, 'no longer single phase'
                    single_phase = False
                    break
                elif flag < 0:
                    print 'ERROR: unable to compute solution, flag', flag
                    break
                else:
                    self.all_sol_u[tstep+1][:self.solverunknowns] = \
                                                            u_retn[-1]
                    self.all_sol_u[tstep+1][self.pos_s] = 0.
                    
                if self.verbose:
                    print 'INFO: pcmmodel at t = ', t_retn[-1]
                tstep += 1
            else:
                raise NotImplementedError
        self.last_sol_tstep = tstep
        
        #self.view_sol(self.times, self.all_sol_u)

    def view_sol(self, times, rTemp):
        """
        Show the solution in Temp with times.
        rTemp[i][:] contains u = r*T at time times[i]
        """
        from fipy import CellVariable, Matplotlib1DViewer, Grid1D
        if rTemp[0][self.pos_s] == 0.:
            meshr = self.state.outer_x_to_r(self.state.outer_gridx, rTemp[0][self.pos_s])
            
            mesh_PCM = Grid1D(dx=tuple(self.state.outer_x_to_r(
                            self.state.outer_delta_x, rTemp[0][self.pos_s])))
            value = rTemp[0][:self.pos_s] / meshr[:]
            solution_view = CellVariable(name = "PCM temperature", 
                mesh = mesh_PCM, value = value[::-1])
            #print 'init', rTemp[0][:self.pos_s] / meshr[:]
        else:
            raise NotImplementedError
            
        self.viewer =  Matplotlib1DViewer(vars = solution_view, datamin=0.)#, datamax=conc.max()+0.20*conc.max())
        self.viewer.plot()
        self.plotevery = 1
        self.viewerplotcount = 0
        for time, rT in zip(times[1:self.last_sol_tstep], rTemp[1:self.last_sol_tstep][:]):
            if rT[self.pos_s] == 0.:
                #meshr = self.state.outer_x_to_r(self.state.outer_gridx, rT[self.pos_s])
                #print' sol', rT[:self.pos_s]/meshr
                value = rT[:self.pos_s]/meshr
                solution_view.setValue(value[::-1])
                self.viewer.axes.set_title('PCM Temp vs radius at time %s' %str(time))
                self.viewer.axes.set_xlabel('Radius')
                self.viewer.axes.set_ylabel('Temp')
                if self.viewerplotcount == 0:
                    self.viewer.plot()
                    #self.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'PCM%sconc%08.4f.png' % (name,time))
                #else:
                #    self.viewer.plot()
                    
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
