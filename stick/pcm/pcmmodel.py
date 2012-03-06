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
        self.L = config.get("pcm.radius")
        self.epsilon = 1e-4
        self.latent_heat_fusion = config.get("pcm.latent_heat_fusion")
        self.solid = {
            'rho' : config.get("pcm.density"),
            'C'   : config.get("pcm.specific_heat_solid"),
            'K'   : config.get("pcm.thermal_cond_solid")
            }
        self.liquid = {
            'rho' : config.get("pcm.density"),
            'C'   : config.get("pcm.specific_heat_liquid"),
            'K'   : config.get("pcm.thermal_cond_liquid")
            }

        #we store computational grid in the state, as we need to swap grids
        #when state changes
        self.n_edge = self.cfg.get('discretization.n_edge') #discretize the PCM radius
        self.grid = None
        self.grid_edge = None
    
    def set_state(self, init_cond):
        """
        init_cond should be a function over (0, R) giving initial temperature
        """
        grid = np.linspace(0., self.R, 100)
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
            self.inner_data = None
        elif self.state == PCMState.LIQUID:
            self.outer_data = self.liquid
            self.inner_data = None
        elif self.state == PCMState.LIQUID_SOLID:
            self.outer_data = self.solid
            self.inner_data = self.liquid
        elif self.state == PCMState.SOLID_LIQUID:
            self.outer_data = self.liquid
            self.inner_data = self.solid

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
        #construct cell centers from this
        self.outer_gridx = (self.outer_gridx_edge[:-1] + self.outer_gridx_edge[1:])/2.
        #obtain cell sizes
        self.outer_delta_x = self.outer_gridx_edge[1:] - self.outer_gridx_edge[:-1]
        
        #conversion functions
        self.outer_x_to_r = lambda x: self.L - (self.L-self.R) * x 
        self.inner_x_to_r = lambda x: self.R * x 

    def reset_state(self, temp_outer):
        """
        reset state will update the state. There are the following possibilities
        1. single phase, en temp_outer is at edges at melting point 
            ==> becomes two phases
        2. two phases, but R becomes close to edge 
            ==> becomes one phase
            
        return value: tuple (changed, switch)
            changed = True if state changed
            switch = True if inner and outer must be switched
        """
        changed = False
        switch = False
        if self.state == PCMState.SOLID:
            if temp_outer[0] >= self.meltpoint:
                self.state = PCMState.SOLID_LIQUID
                changed = True
                switch = True
        elif 

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
        pass

    def solve_init(self):
        """
        Initialize the solver so they can be solved stepwize
        """
        pass

    def solve(self):
        """
        Initialize the solver so they can be solved stepwize
        """
        pass

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
