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

        self.n_edge = self.cfg.get('discretization.n_edge') #discretize the PCM radius

        self.initialized = False

        self.plotevery = self.cfg.get("plot.plotevery")

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh (in mm) on which we project results
        grid: the space position of each central point for every cell element;
        """
        self.grid_edge = sp.linspace(0.,
            self.cfg.get('pcm.radius'), self.cfg.get('discretization.n_edge'))
        self.init_temp = eval(self.cfg.get('init.init_temp'))
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]

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
