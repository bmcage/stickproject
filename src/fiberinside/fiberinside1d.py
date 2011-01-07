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
Module holding a generic diffusion model in a single fiber. 
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
from scipy.integrate import odeint, ode, trapz
import matplotlib.pyplot as plt
import sets
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
import yarn2d.config as conf
from mycorrection import MyDiffusionTermNoCorrection
from yarn2dfiber import Yarn2dFiber
from fiberfipy.config import FiberfipyConfigManager
from fiberfipy.fibermodel import FiberModel

#-------------------------------------------------------------------------
#
#Fipy Imports
#
#-------------------------------------------------------------------------
from fipy import *

class FiberInside(object):
    """
    fiberinside1d.FiberInside is a special diffusion model for single radial fiber
    which is composed of a specific material (cotton, polyester,.etc). In this model
    the diffusion in the fiber is considered
    """
    def __init__(self,config):
        self.datatime = []
        self.cfg = config
        self.method = self.cfg.get('general.method')
        if not (self.method in ['FVM']):
            print 'unknown solution method'
            sys.exit(0)
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = self.time_period / self.delta_t
        self.n_edge_in = self.cfg.get('fiber.n_edge_in') #discretize the fiber radius
        self.beginning_surf = self.cfg.ge('fiber.radius_pure_fiber')
        self.diff_coeff_in = self.cfg.get('fiberin.diffusion_inside')
    
    def create_mesh():
        n_edge_in = [0]
        self.total_edges_in = 0
        first = True
        n_edge_in += [self.n_edge]
        for nr in n_edge_in:
            if nr == 0 and self.total_edges_in == 0 and not first:
                print "Error, no discretization is given"
                sys.exit(0)
            if not (nr == 0):
                first = False 
                if self.total_edges_in:
                    self.total_edges_in += nr -1
                else:
                    self.total_edges_in = nr
        self.grid_edge_in = sp.empty(self.total_edges_in, float)
        left_in = 0.
        totnr_in = 0
        first = True
        for nr, right_in in zip(n_edge_in, self.beginning_surf):
            if nr:
                if first:
                    self.grid_edge_in[totnr_in:totnr_in + nr] = sp.linspace(left_in, right_in, nr)
                    totnr_in += nr
        self.grid_in = (self.grid_edge_in[:-1] + self.grid_edge_in[1:])/2.
        self.delta_r_in = self.grid_in[:-1]-self.grid_in[1:]
        self.mesh_fiber_in = CylindericalGrid1D(dr=tuple(self.delta_r_in))
        self.mesh_fiber_in.periodicBC = False
        
    def initial_fiber_in(self):
        """Initial concentration over the domain"""
        self.initial_c_in = sp.empty(self.total_edges_in - 1, float)
        self.initial_c_in[:] = self.cfg.get('fiberin.initial_in')
        
    def calc_mass(self,conc_r_in):
        """calculate mass of component present"""
    
    def solve_in_fipy(self):
        #using fipy to solve the 1D problem in fiber (DEET diffusion in fiber)
        discretization_t = self.steps + 1
        self.times = sp.linspace(0, self.time_period, discretization_t)
        self.solution_in_fiber = CellVariable(name = "Conc. Active component in fiber",
                                    mesh = self.mesh_fiber, 
                                    value = self.initial_c_in, hasOld = 1)
        self.viewer = Viewer(vars = self.solution_in_fiber, datamin = 0, datamax = 1.0)
        self.conc_c_in = sp.empty(len(self.times), len(self.initial_c_in), float)
        print 'the length of solution in fiber', len(self.solution_in_fiber)
        print 'the length of grid', len(self.gird)
        
        #set up the boundary condition
        if self.boundary_in_right:
            self.BCs_fiber_in = (FixedValue(faces = self.mesh_fiber_in.getFacesRight(),
                                    value = self.boundary_in_right), 
                                FixedValue(faces = self.mesh_fiber_in.getFacesLeft(),
                                    value = self.boundary_in_left))
        self.eqX_fiber_in = TransientTerm() == DiffusionTerm(coeff = self.diff_coeff_in)
        tsetp = 0
        self.conc_c_in[tstep][:] = self.initial_c_in
        self.solve_fipy_in_step()
        if self.viewer is not None:
            self.viewer.plot()
            #raw_input("please<return>.....")
        tstep += 1
        self.conc_c_in[tstep][:] = self.solution_in_fiber.getValue()
        

    
    def run(self):
        self.creat_mesh()
        self.initial_c_in()
        self.solution_in_fiber()