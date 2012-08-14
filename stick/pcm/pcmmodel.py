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
from copy import copy

HAVE_ODES = False
try:
    from scikits.odes import ode as sc_ode
    HAVE_ODES = True
except:
    print 'Could not load scikits.odes, odes solver not available'
    sys.exit(0)

from scikits.odes.sundials.cvode import CV_RootFunction

CV_SUCCESS = 0
CV_TSTOP_RETURN = 1
CV_ROOT_RETURN = 2
#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
from lib.utils.utilsbm import deriv133, inter3, inter2

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
                state = PCMState.INVALID
                break
            elif state == PCMState.SOLID_LIQUID and init_cond(rpos) < self.meltpoint:
                # this is not supported in the model, two interfaces
                state = PCMState.INVALID
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
        Determines the fixed grid on the unit interval for the current state
        """
        #nonuniform partition
        UNIFORM = 1 
        REFINTF = 2
        REFBOTH = 3
        METHOD = REFBOTH
        if METHOD == UNIFORM:
            #uniform x grid
            self.inner_gridx_edge = sp.linspace(0., 1., self.n_edge)
            self.outer_gridx_edge = sp.linspace(0., 1., self.n_edge)
        elif METHOD == REFINTF:
            #refined grid at interface
            self.inner_gridx_edge = np.empty(self.n_edge, float)
            self.outer_gridx_edge = np.empty(self.n_edge, float)
            dh = np.empty(self.n_edge-1, float)
            N = self.n_edge - 1
            for j in range(N):
                dh[j] = 2/(N+1) * (N-j)/N
                self.inner_gridx_edge[j+1] = self.inner_gridx_edge[j] + dh[j]
                self.outer_gridx_edge[j+1] = self.outer_gridx_edge[j] + dh[j]
        elif METHOD == REFBOTH:
            #Edge and interface refined
            self.inner_gridx_edge = np.empty(self.n_edge, float)
            self.outer_gridx_edge = np.empty(self.n_edge, float)
            if self.n_edge // 2 == 0:
                print 'ERROR: n_edge should be odd !'
                sys.exit()
            dh = np.empty(self.n_edge-1, float)
            N = self.n_edge - 1
            Nh = N//2
            for j in range(Nh):
                dh[Nh-1-j] = (2/(Nh+1) * (Nh-j)/Nh)/2.
                dh[Nh+j] = dh[Nh-1-j]
            self.inner_gridx_edge[0] = 0.
            self.outer_gridx_edge[0] = 0.
            for j in range(N):
                self.inner_gridx_edge[j+1] = self.inner_gridx_edge[j] + dh[j]
                self.outer_gridx_edge[j+1] = self.outer_gridx_edge[j] + dh[j]
        else:
            print 'ERROR: Unknown GRID type'
            sys.exit()
        #make sure 1 is 1!
        self.inner_gridx_edge[-1] = 1.
        self.outer_gridx_edge[-1] = 1.
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
        self.minval = L*eps/2.  # used also elsewhere in the code !!
        self.maxval = L*(1-eps/2.)

    def evaluate(self, t, u, out, userdata):
        out[0] = u[self.pos_s] - self.minval
        out[1] = self.maxval - u[self.pos_s]
        #print 'test root double ph', u[self.pos_s]/self.L, out[0], out[1]
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

        self.__cont = False
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
        self.lambda_m = self.state.latent_heat_fusion
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

        Rlast = self.state.outer_x_to_r(self.state.outer_gridx[0])
        #print 't', t, u_rep, u_rep/self.state.outer_x_to_r(self.state.outer_gridx)

        #center point of PCM, is stored at last point
        flux_edge[-1] = -K/rho/Cv * dxdr * u_rep[-1] \
                        /self.state.outer_x_to_r(self.state.outer_gridx[-1])
        #edge of PCM, transfer cond, stored at first point
        #print u_rep[0]/Rlast - self.Tout(t), u_rep[0]
        flux_edge[0] = (self.hT / rho / Cv * dxdr*L*
                         (u_rep[0]/Rlast - self.Tout(t))
                        -K/rho/Cv * dxdr*u_rep[0] /Rlast)
        flux_edge[1:-1] = -K/rho/Cv * (u_rep[1:]-u_rep[:-1])/ Dxavg[:]*dxdr**2
        diff_u_t[:] = (flux_edge[:-1]-flux_edge[1:])/Dx[:]

    def f_odes_dph_FD(self, t, u_rep, diff_u_t):
        """ RHS function for the cvode solver for double phase energy diffusion
        The grid is in x with first x from 0 to 1 with outer data, so 0 is 
        at L and 1 is at interface, then the interface position, then x from 1
        to 0 with innter data, with 1 the interface, and 0 the center
        """
        gridi = self.state.inner_gridx
        gridi_edge = self.state.inner_gridx_edge
        ui = u_rep[:self.pos_s]
        s = u_rep[self.pos_s]
        #we store outer from highest x to lowest x!
        uo =  u_rep[self.pos_s+1:][::-1]
        Dxi = self.state.inner_delta_x
        Dxiavg = self.state.inner_delta_x_avg
        grido = self.state.outer_gridx
        grido_edge = self.state.outer_gridx_edge
        Dxo = self.state.outer_delta_x
        Dxoavg = self.state.outer_delta_x_avg

        n_celli = len(gridi) 
        n_cello = len(grido) 
        flux_edgei = self.__tmp_flux_edge[:n_celli+1]
        flux_edgeo = self.__tmp_flux_edge[n_celli+1:n_celli+n_cello+2]
        L = self.state.L
        Cvo = self.state.outer_data['C']
        Ko = self.state.outer_data['K']
        rhoo = self.state.outer_data['rho']
        Cvi = self.state.inner_data['C']
        Ki = self.state.inner_data['K']
        rhoi = self.state.inner_data['rho']
        
        Tmelt = self.state.meltpoint
        if uo[0]/self.state.outer_x_to_r(self.state.outer_gridx[0],s) > self.Tout(t):
            print 'WRONG', t, uo[0]/self.state.outer_x_to_r(self.state.outer_gridx[0],s), '>',  self.Tout(t)
        ao = Ko/rhoo/Cvo
        ai = Ki/rhoi/Cvi

        dxdri =  1./s
        dxdro = -1./(L-s)
        
        #at interface
        dudxi_ats = deriv133(gridi[-2],ui[-2],gridi[-1],ui[-1],1.,s*Tmelt)
        dudxo_ats = deriv133(grido[-2],uo[-2],grido[-1],uo[-1],1.,s*Tmelt)

        #equation sdot
        ##sdot = 1./rhoo/self.lambda_m *(
        ##        Ki*(dudxi_ats*dxdri/s - Tmelt/s)
        ##      - Ko*(dudxo_ats*dxdro/s - Tmelt/s))
        ##sdot2 = 1./rhoo/self.lambda_m*\
        ##         (Ki*dxdri*(Tmelt-ui[-1]/s/gridi[-1])/(1-gridi[-1]) 
        ##        - Ko*dxdro*(Tmelt-uo[-1]/(L-(L-s)*grido[-1]))/(1-grido[-1]))
        sdot = 1/self.lambda_m /s*(ai*Cvi*dxdri*dudxi_ats-ao*Cvo*dxdro*dudxo_ats) \
                -Tmelt/self.lambda_m/ s**3 *(ai*Cvi - ao*Cvo)
        ##sdot = sdot3

        #average value of dxdt at the gridpoints of centers
        dxdtic = - gridi / s * sdot
        dxdtoc = grido / (L-s) * sdot
        dxdti_ats = -sdot / s
        dxdto_ats = sdot / (L-s)

        # at center point of PCM, dTdr=0, so x=0 for inner
        r_center = gridi[0]
        Tcenter =  inter2(gridi[0], ui[0] /self.state.inner_x_to_r(r_center, s),
                            gridi[1], ui[1] /self.state.inner_x_to_r(gridi[1], s),
                            0.)
        flux_edgei[0] = -ai * dxdri * Tcenter
        ## alternative form:
        ##flux_edgei[0] = -Ki/rhoi/Cvi * dxdri * ui[0] \
        ##                /self.state.inner_x_to_r(self.state.inner_gridx[0],s)

        #inner part
        flux_edgei[1:self.pos_s] = -ai * (ui[1:]-ui[:-1])/ Dxiavg[:]*dxdri**2
        #at interface
        flux_edgei[self.pos_s] = -ai * dudxi_ats *dxdri**2
        diff_u_t[:self.pos_s] = (flux_edgei[:-1]-flux_edgei[1:])/Dxi[:]
        #FD part
        # 1. use flux_edgei to store u_edgei
        u_edgei = flux_edgei
        u_edgei[0] = 0. # at r=0, u =r*T is 0! ui[0]
        u_edgei[-1] = s * Tmelt
        u_edgei[1:-1] = inter2(gridi[:-1], ui[:-1], gridi[1:], ui[1:], gridi_edge[1:-1])#(ui[:-1] + ui[:-1])/2. ##TODO should we not average T ??
        # 2. add the FD part of dxdt term
        diff_u_t[:self.pos_s] += sdot * dxdri * \
                    ((gridi_edge[1:] * u_edgei[1:] - gridi_edge[:-1] * u_edgei[:-1])/Dxi[:]
                    - (u_edgei[:-1]/2 + u_edgei[1:]/2))
        ## alternative 1
        ##diff_u_t[:self.pos_s] += dxdtic[:]/Dxi[:] * (u_edgei[:-1] - u_edgei[1:])
        ## alternative 2
        ##diff_u_t[:self.pos_s] += sdot * dxdri * \
        ##            ((gridi_edge[1:] * u_edgei[1:] - gridi_edge[:-1] * u_edgei[:-1])/Dxi[:]
        ##            - ui[:])

        #interface equation
        diff_u_t[self.pos_s] = sdot
        
        # at edge of PCM, transfer cond, so x=0 for outer
        r_edge = self.state.outer_x_to_r(self.state.outer_gridx[0], s)
        flux_edgeo[0] = (self.hT / rhoo / Cvo * dxdro*L*
                         (uo[0]/r_edge - self.Tout(t))
                        -ao * dxdro * uo[0] /r_edge)
    
        if flux_edgeo[0] < 0.:
            print  'flux_edge < 0,', flux_edgeo[0], uo[0]/r_edge, self.Tout(t)
        #outer part
        flux_edgeo[1:self.pos_s] = -ao * (uo[1:]-uo[:-1])/ Dxoavg[:]*dxdro**2
        flux_edgeo[self.pos_s] = -ao * dudxo_ats *dxdro**2
        tmp = (flux_edgeo[:-1]-flux_edgeo[1:])/Dxo[:]
        #FD part
        # 1. use flux_edgeo to store u_edgeo
        u_edgeo = flux_edgeo
        u_edgeo[0] = uo[0] # at r=L, u =L*T 
        u_edgeo[-1] = s * Tmelt
        u_edgeo[1:-1] = (uo[:-1] + uo[:-1])/2. ##TODO should we not average T ??
        # 2. add the FD part of dxdt term
        tmp[:self.pos_s] += dxdtoc[:]/Dxo[:] * (u_edgeo[:-1] - u_edgeo[1:])
        #we store it inverse, like the physical PCM in order !
        diff_u_t[self.pos_s+1:] = tmp[::-1]

    def f_odes_dph_FV(self, t, u_rep, diff_u_t):
        """ RHS function for the cvode solver for double phase energy diffusion
        The grid is in x with first x from 0 to 1 with outer data, so 0 is 
        at L and 1 is at interface, then the interface position, then x from 1
        to 0 with innter data, with 1 the interface, and 0 the center
        """
        print "ERROR , Finite Volume way not working"
        sys.exit(1)

        gridi = self.state.inner_gridx
        gridi_edge = self.state.inner_gridx_edge
        ui = u_rep[:self.pos_s]
        s = u_rep[self.pos_s]
        #we store outer from highest x to lowest x!
        uo =  u_rep[self.pos_s+1:][::-1]
        Dxi = self.state.inner_delta_x
        Dxiavg = self.state.inner_delta_x_avg
        grido = self.state.outer_gridx
        grido_edge = self.state.outer_gridx_edge
        Dxo = self.state.outer_delta_x
        Dxoavg = self.state.outer_delta_x_avg

        n_celli = len(gridi) 
        n_cello = len(grido) 
        flux_edgei = self.__tmp_flux_edge[:n_celli+1]
        flux_edgeo = self.__tmp_flux_edge[n_celli+1:n_celli+n_cello+2]
        L = self.state.L
        Cvo = self.state.outer_data['C']
        Ko = self.state.outer_data['K']
        rhoo = self.state.outer_data['rho']
        Cvi = self.state.inner_data['C']
        Ki = self.state.inner_data['K']
        rhoi = self.state.inner_data['rho']
        
        Tmelt = self.state.meltpoint
        if uo[0]/self.state.outer_x_to_r(self.state.outer_gridx[0],s) > self.Tout(t):
            print 'WRONG', t, uo[0]/self.state.outer_x_to_r(self.state.outer_gridx[0],s), '>',  self.Tout(t)
        ao = Ko/rhoo/Cvo
        ai = Ki/rhoi/Cvi

        dxdri =  1./s
        dxdro = -1./(L-s)
        
        #at interface
        dudxi_ats = deriv133(gridi[-2],ui[-2],gridi[-1],ui[-1],1.,s*Tmelt)
        dudxo_ats = deriv133(grido[-2],uo[-2],grido[-1],uo[-1],1.,s*Tmelt)

        #equation sdot
        sdot = 1./rhoo/self.lambda_m *(
                Ki*(dudxi_ats*dxdri/s - Tmelt/s)
              - Ko*(dudxo_ats*dxdro/s - Tmelt/s))
        sdot2 = 1./rhoo/self.lambda_m*\
                 (Ki*dxdri*(Tmelt-ui[-1]/s/gridi[-1])/(1-gridi[-1]) 
                - Ko*dxdro*(Tmelt-uo[-1]/(L-(L-s)*grido[-1]))/(1-grido[-1]))
        sdot3 = 1/self.lambda_m /s*(ai*Cvi*dxdri*dudxi_ats-ao*Cvo*dxdro*dudxo_ats) \
                -Tmelt/self.lambda_m/ s**3 *(ai*Cvi - ao*Cvo)
        sdot = sdot3

        #average value of dxdt at the edge, take value at edge point
        dxdti = - gridi_edge[1:-1] / s * sdot
        dxdto = grido_edge[1:-1] / (L-s) * sdot
        dxdti_ats = -sdot / s
        dxdto_ats = sdot / (L-s)

        # at center point of PCM, dTdr=0, so x=0 for inner
        r_center = gridi[0]
        Tcenter =  (ui[0] /self.state.inner_x_to_r(r_center,s))
        flux_edgei[0] = -ai * dxdri * Tcenter
        #inner part
        flux_edgei[1:self.pos_s] = -ai * (ui[1:]-ui[:-1])/ Dxiavg[:]*dxdri**2\
                    + dxdti * (ui[:-1]+ui[1:])/2
        #at interface
        flux_edgei[self.pos_s] = -ai * dudxi_ats *dxdri**2\
                    + dxdti_ats * s * Tmelt
        diff_u_t[:self.pos_s] = (flux_edgei[:-1]-flux_edgei[1:])/Dxi[:]
        
        #interface equation
        diff_u_t[self.pos_s] = sdot
        
        # at edge of PCM, transfer cond, so x=0 for outer
        r_edge = self.state.outer_x_to_r(self.state.outer_gridx[0], s)
        flux_edgeo[0] = (self.hT / rhoo / Cvo * dxdro*L*
                         (uo[0]/r_edge - self.Tout(t))
                        -ao * dxdro * uo[0] /r_edge)
    
        if flux_edgeo[0] < 0.:
            print  'flux_edge < 0,', flux_edgeo[0], uo[0]/r_edge, self.Tout(t)
        #outer part
        flux_edgeo[1:self.pos_s] = -ao * (uo[1:]-uo[:-1])/ Dxoavg[:]*dxdro**2\
                    + dxdto * (uo[:-1]+uo[1:])/2
        flux_edgeo[self.pos_s] = -ao * dudxo_ats *dxdro**2\
                    + dxdto_ats * s * Tmelt
        tmp = (flux_edgeo[:-1]-flux_edgeo[1:])/Dxo[:]
        #we store it inverse, like the physical PCM in order !
        diff_u_t[self.pos_s+1:] = tmp[::-1]

    def solve_odes_reinit(self):
        """
        Reinitialize the cvode solver to start again
        """
        self.initial_t = self.times[0]
        print 'REINIT SOLVER'
        del self.solver
        if self.state.single_phase():
            print ' ... SINGLE'
            self.solverunknowns = self.pos_s
            rootfn=RootFnsp()
            rootfn.set_data(self.state.outer_x_to_r(self.state.outer_gridx[0]), 
                            self.state.meltpoint)
            self.solver = sc_ode('cvode', self.f_odes_sph,
                                 max_steps=50000, lband=1, uband=1,
                                 nr_rootfns=1, rootfn=rootfn,
                                 first_step_size=1e-12)
        else:
            print ' ... DOUBLE'
            self.solverunknowns = self.nrunknowns
            rootfn=RootFndp()
            rootfn.set_data(self.state.L, self.pos_s, self.state.epsilon)
            self.solver = sc_ode('cvode', self.f_odes_dph_FD,
                                 max_steps=50000, 
                                 #lband=1, uband=1,
                                 nr_rootfns=2, rootfn=rootfn,
                                 atol=1e-12, rtol=1e-12,
                                 first_step_size=1e-12)

##        self.solver.init_step(self.step_old_time, self.step_old_sol)

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
        self.__tmp_flux_edge = sp.empty(self.nrunknowns+2, float)

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
                self.solver.set_options(tcrit=time)
                flag, t_retn, u_retn, t_out, u_last = self.solver.solve(
                                [prev_time,time], 
                                initval
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
                self.solver.set_options(tcrit=time)
                flag, t_retn, u_retn, t_out, u_last = self.solver.solve(
                                [prev_time,time], 
                                initval,
                                )
                
                if flag == CV_ROOT_RETURN:
                    print 'At time', t_out, 'no longer double phase, s at', u_last[self.pos_s]/self.state.L
                    #store the last solution
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_last
                    #self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_out
                    single_phase = True
                    changed = True
                elif flag == CV_TSTOP_RETURN:
                    #here crit time is just the last good time we want
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_last
                    #self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_out
                elif flag == CV_SUCCESS:
                    self.all_sol_u[tstep+1][:self.solverunknowns] = u_retn[-1]
                    #self.all_sol_u[tstep+1][self.pos_s] = 0.
                    tret = t_retn[-1]
                else:
                    print 'ERROR: unable to compute solution, flag', flag
                    break
                    
                if self.verbose:
                    print 'INFO: pcmmodel at t = ', tret
                if not single_phase:
                    tstep += 1

            if tstep == len(self.times)-1:
                print 'End reached'
                break
            if changed:
                # we need to swap the solver we use!
                self.step_old_time = tret
                if single_phase:
                    #we went from double phase to single phase
                    #1. determine what phase
                    #2. change solver
                    #3. update data
                    if abs(self.all_sol_u[tstep+1][self.pos_s]) \
                        < self.state.L*self.state.epsilon/2. * 1.001:
                        # now only outer state
                        self.state.state = self.state.outer_data['state']
                        self.state.inner_data = self.state.outer_data
                        ## obscure error with passing a np.array after double phase
                        ## doing a copy fixes it
                        self.step_old_sol = copy(self.project_outdp_sp(self.all_sol_u[tstep+1]))
                    elif abs(self.state.L - 
                             self.all_sol_u[tstep+1][self.pos_s]) \
                            < self.state.epsilon:
                        # now only inner state
                        self.state.state = self.state.inner_data['state']
                        self.state.outer_data = self.state.inner_data
                        self.step_old_sol = copy(self.project_indp_sp(self.all_sol_u[tstep+1]))
                    else:
                        raise NotImplementedError, 'We should not reach this'
                    # no more interface:
                    self.state.R = 0.
                else:
                    #we went from single phase to double phase, interface
                    #arises at R=L*(1-eps)
                    #1. determine what phase
                    #2. change solver
                    #3. update data
                    self.state.R = self.state.L*(1 - self.state.epsilon)
                    if self.state.outer_data['state'] == PCMState.LIQUID:
                        self.state.state = PCMState.LIQUID_SOLID
                        self.state.inner_data = self.state.outer_data
                        self.state.outer_data = self.state.solid
                    elif self.state.outer_data['state'] == PCMState.SOLID:
                        self.state.state = PCMState.SOLID_LIQUID
                        self.state.inner_data = self.state.outer_data
                        self.state.outer_data = self.state.liquid
                        print 'to SOLID_LIQUID', self.state.R
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
        print 'introducing interface, project out sp to dp'
        #The outer solution is set at meltingtemp
        newdataC = np.zeros(len(self.state.outer_gridx)+1+len(self.state.inner_gridx), float)
        newrout = self.state.outer_x_to_r(self.state.outer_gridx, R=R)
        for i in xrange(self.pos_s):
            newdataC[self.pos_s+1:] = self.state.meltpoint \
                                    * (1+self.state.epsilon) * newrout[::-1]
        #new interface is set
        newdataC[self.pos_s] = R
        #and we project the old outer solution to inner solution, which we 
        # store inverted
##        # SIMPLE: assume projection can be neglected
##        newdataC[:self.pos_s] = sol[:self.pos_s][::-1]/self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=0.) \
##                                * self.state.inner_x_to_r(self.state.inner_gridx, R=R)
        
        origr = self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=0.)
        origredge = self.state.outer_x_to_r(self.state.outer_gridx_edge[::-1], R=0.)
        #in one phase mode, sol is in inverse r order
        uo = sol[:self.pos_s][::-1]
        # conserved is volumeshell T = volumeshell U/r
        origdataC = uo / origr * \
                        4/3*np.pi * (origredge[:-1]**3 - origredge[1:]**3)
        newr = self.state.inner_x_to_r(self.state.inner_gridx[:], R=R)
        newredge = self.state.inner_x_to_r(self.state.inner_gridx_edge[:], R=R)
        pos = 0
        posorig=0
        for pos in np.arange(self.pos_s):
            while origredge[posorig] < newredge[pos+1]:
                if origredge[posorig+1] >= newredge[pos+1]:
                    if origredge[posorig] >= newredge[pos]:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (newredge[pos+1]**3 - origredge[posorig]**3)
                    else:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (newredge[pos+1]**3 - newredge[pos]**3**3)
                    ind = 0
                    while (pos+1+ind < self.pos_s) and \
                            (origredge[posorig+1] > newredge[pos+1+ind]):
                        if origredge[posorig+1] > newredge[pos+1+ind+1]:
                            newdataC[pos+1+ind] = uo[posorig]/origr[posorig] * \
                                4/3*np.pi * (newredge[pos+1+ind+1]**3 - newredge[pos+1+ind]**3)
                        else:
                            newdataC[pos+1+ind] = uo[posorig]/origr[posorig] * \
                                4/3*np.pi * (origredge[posorig+1]**3 - newredge[pos+1+ind]**3)
                        ind += 1
                else:
                    nextC = 0.
                    if origredge[posorig] >= newredge[pos]:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (origredge[posorig+1]**3 - origredge[posorig]**3)
                    else:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (newredge[pos+1]**3 - newredge[pos]**3**3)
                posorig += 1
        #now we derive the new u = r T_avg
        newdataC[:self.pos_s] = newdataC[:self.pos_s] * newr \
                    / (4/3*np.pi * (newredge[1:]**3 - newredge[:-1]**3))
##        print 'before T', sol[:self.pos_s][::-1]/self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=0.)
##        print '  over', self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=0.) / self.state.L
##        print 'after T', R, newdataC[:self.pos_s] / self.state.inner_x_to_r(self.state.inner_gridx, R=R)
##        print '  over', self.state.inner_x_to_r(self.state.inner_gridx, R=R) / self.state.L
##        print '      T', newdataC[self.pos_s+1:]/self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=R)
##        print '  over', self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=R) / self.state.L
        return newdataC

    def project_outdp_sp(self, sol):
        """ project a double phase solution sol, to a single phase solution
            over the outer grid, assuming interface goes to 0 """
        # x order is inverse of the r order
        origr = self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=sol[self.pos_s])
        origredge = self.state.outer_x_to_r(self.state.outer_gridx_edge[::-1], R=sol[self.pos_s])
        #in two phase mode, sol is in correct r order
        uo = sol[self.pos_s+1:]
        # conserved is volumeshell T = volumeshell U/r
        origdataC = uo / origr * \
                        4/3*np.pi * (origredge[:-1]**3 - origredge[1:]**3)
        newr = self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=0.)
        newredge = self.state.outer_x_to_r(self.state.outer_gridx_edge[::-1], R=0.)
        newdataC = np.empty(self.pos_s, float)
        pos = 0
        posorig=0
        #remaining piece goes to meltpoint temperature
        newdataC[pos] = self.state.meltpoint * 4/3*np.pi * sol[self.pos_s]**3
        for pos in np.arange(self.pos_s):
            while origredge[posorig] < newredge[pos+1]:
                if origredge[posorig+1] >= newredge[pos+1]:
                    if origredge[posorig] >= newredge[pos]:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (newredge[pos+1]**3 - origredge[posorig]**3)
                    else:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (newredge[pos+1]**3 - newredge[pos]**3**3)
                    ind = 0
                    while (pos+1+ind < self.pos_s) and \
                            (origredge[posorig+1] > newredge[pos+1+ind]):
                        if origredge[posorig+1] > newredge[pos+1+ind+1]:
                            newdataC[pos+1+ind] = uo[posorig]/origr[posorig] * \
                                4/3*np.pi * (newredge[pos+1+ind+1]**3 - newredge[pos+1+ind]**3)
                        else:
                            newdataC[pos+1+ind] = uo[posorig]/origr[posorig] * \
                                4/3*np.pi * (origredge[posorig+1]**3 - newredge[pos+1+ind]**3)
                        ind += 1
                else:
                    nextC = 0.
                    if origredge[posorig] >= newredge[pos]:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (origredge[posorig+1]**3 - origredge[posorig]**3)
                    else:
                        newdataC[pos] += uo[posorig]/origr[posorig] * \
                            4/3*np.pi * (newredge[pos+1]**3 - newredge[pos]**3**3)
                posorig += 1
        #now we derive the new u = r T_avg
        newdataC = newdataC * newr \
                    / (4/3*np.pi * (newredge[1:]**3 - newredge[:-1]**3))
        #we store outer data in one phase in the x order, so inverse from r
        newdataC = newdataC[::-1]
##        print 'before T', sol[:self.pos_s]/self.state.inner_x_to_r(self.state.inner_gridx[:], R=sol[self.pos_s])
##        print '  over', self.state.inner_x_to_r(self.state.inner_gridx[:], R=sol[self.pos_s]) / self.state.L
##        print ' and   T', sol[self.pos_s+1:]/self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=sol[self.pos_s])
##        print '  over', self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=sol[self.pos_s]) / self.state.L
##        print 'after T', newdataC[::-1]/self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=0.)
##        print '  over', self.state.outer_x_to_r(self.state.outer_gridx[::-1], R=0.) / self.state.L
        return newdataC

    def project_indp_sp(self, sol):
        """ project a double phase solution sol, to a single phase solution
            over the outer grid, assuming interface goes to L, so inner part
            of double phase needs to be considered
        """
        raise NotImplementedError, 'project double to single in'
        print 'interface at', sol[self.pos_s], ' projecting'
        origr = self.state.inner_x_to_r(self.state.inner_gridx, R=sol[self.pos_s])
        origredge = self.state.inner_x_to_r(self.state.inner_gridx_edge, R=sol[self.pos_s])
        # conserved is volumeshell T = volumeshell U/r
        ui = sol[:self.pos_s] 
        origdataC = ui / origr * \
                        4/3*np.pi * (origredge[:-1]**3 - origredge[1:]**3)
        newr = self.state.outer_x_to_r(self.state.outer_gridx, R=0.)
        newredge = self.state.outer_x_to_r(self.state.outer_gridx_edge, R=0.)
        newdataC = np.empty(self.pos_s, float)
        pos = 0
        posorig=self.pos_s-1
        #remaining piece goes to meltpoint temperature
        newC = self.state.meltpoint * 4/3*np.pi * (self.state.L**3-sol[self.pos_s]**3)
        nextC = 0.
        
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
        import pylab
        
        pylab.ion()
        figaspect=1.0
        w, h = pylab.figaspect(figaspect)
        self.__fig = pylab.figure(figsize=(w, h))
        self.__figaxes = pylab.gca()
        self.__figid = 1
        self.__figaxes.set_title('Temp over PCM radius')
        
        if rTemp[0][self.pos_s] == 0.:
            meshr = self.state.outer_x_to_r(self.state.outer_gridx[::-1],
                                                        rTemp[0][self.pos_s])
            solr = rTemp[0][:self.pos_s][::-1] / meshr[:]
            
        else:
            raise NotImplementedError
        if self.plotevery:
            self.__figaxes.set_ylim(ymin=0., ymax=45.)
            self.__figaxes.axes.set_xlabel('Radius')
            self.__figaxes.axes.set_ylabel('Temp')
            self.__figaxes.plot(meshr, solr)
        self.viewerplotcount = 0
        for time, rT in zip(times[1:self.last_sol_tstep], rTemp[1:self.last_sol_tstep][:]):
            if rT[self.pos_s] == 0.:
                meshr = self.state.outer_x_to_r(self.state.outer_gridx[::-1], 
                                                        rT[self.pos_s])
                if self.plotevery and self.viewerplotcount == 0:
                    solr = rT[:self.pos_s][::-1]/meshr
                    self.__figaxes.clear()
                    self.__figaxes.set_ylim(ymin=0., ymax=45.)
                    self.__figaxes.axes.set_xlabel('Radius')
                    self.__figaxes.axes.set_ylabel('Temp')
                    self.__figaxes.plot(meshr, solr, 'r')
                    self.__figaxes.set_title('PCM Temp vs radius at time %s' %str(time))
                    pylab.draw()
                if self.plotevery:
                    self.viewerplotcount += 1
                    self.viewerplotcount = self.viewerplotcount % self.plotevery
            else:
                if self.plotevery and self.viewerplotcount == 0:
                    meshr1 = self.state.inner_x_to_r(self.state.inner_gridx, 
                                                        rT[self.pos_s])
                    meshr2 = self.state.outer_x_to_r(self.state.outer_gridx[::-1], 
                                                        rT[self.pos_s])
                    solr1 = rT[:self.pos_s]/meshr1
                    solr2 = rT[self.pos_s+1:]/meshr2
                    self.__figaxes.clear()
                    self.__figaxes.set_ylim(ymin=0., ymax=45.)
                    self.__figaxes.axes.set_xlabel('Radius')
                    self.__figaxes.axes.set_ylabel('Temp')
                    self.__figaxes.plot(meshr1, solr1, 'r', meshr2, solr2, 'r')
                    self.__figaxes.set_title('PCM Temp vs radius at time %s' %str(time))
                    pylab.draw()
                if self.plotevery:
                    self.viewerplotcount += 1
                    self.viewerplotcount = self.viewerplotcount % self.plotevery

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
