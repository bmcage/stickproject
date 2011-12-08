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
    upscaling to fabric1d domain
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
from scipy.integrate import ode
from scipy.integrate import vode
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
import bednet.config as conf
from mycorrection import MyDiffusionTermNoCorrection
from yarn2dgrid import Yarn2dGrid
from Grid2D import Grid2D
from yarn2dfiber import Yarn2dFiber
from yarn2d.config import *
from fiberfipy.config import FiberfipyConfigManager
from fiberfipy.fibermodel import FiberModel

#-------------------------------------------------------------------------
#import fipy
#-------------------------------------------------------------------------
from fipy import *

#-------------------------------------------------------------------------
#
# DiffusionModel-fabric1d class 
#
#-------------------------------------------------------------------------
class Bednet(object):
    """
    Fabirc2dUpscaling is a special diffusion model to calculate the concentration of 
    DEET outside the fabric. This upscaling method is simple with the equation:
    C_outside = C_yarn * R_yarn/d_distance
    """
    def __init__(self,cfg):
        self.datatime = []
        self.cfg = config
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = self.time_period / self.delta_t
        self.domain_size = self.cfg.get('domain.domain_size')
        self.dx = self.cfg.get('domain.dx')
        self.dy  = self.cfg.get('domain.dy')
        self.size_sample = self.cfg.get('sample.size_sample')
        #self.thickness_sample = self.cfg.get('sample.thickness_sample')
        self.boundary_up = self.cfg.get('boundary.boundary_up')
        self.boundary_bottom = self.cfg.get('boundary.boundary_bottom')
        self.boundary_left = self.cfg.get('boundary.boundary_left')
        self.boundary_right = self.cfg.get('boundary.boundary_right')
        self.diffusion_coef_DEET = self.cfg.get('diffusion_DEET.diffusion_coef_DEET')
        self.saturation_conc = self.cfg.get('saturation.saturation_conc')

        """
        self.distance_yarn = self.cfg.get('domain.distance_yarn')
        self.grid_space_vertical = self.cfg.get('domain.grid_space_vertical')
        self.number_nodes = (self.domain_size[0] / self.distance_yarn + 1) * \
                            (self.domain_size[1] / self.grid_space_vertical + 1)
        """
        self.cfg_fiber = []        
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                        os.path.dirname(self.cfg.filename), filename))
            self.cfg_fiber.append(Fiber1dConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.time_period)
            self.cfg_fiber[-1].set("time.dt", self.delta_t)  
            
        #create yarn models
        self.yarn_models = []
        for cfg in self.cfg_yarn:
            self.yarn_model.append(Yarn2DModel(cfg))
        self.nr_models = len(self.yarn_models)
        self.verbose = self.cfg.get('general.verbose')
                
        #some memory
        self.cache_index_t_bednet = 0
        self.cache_index_t_yarn = [0] * self.nr_models
        
        #plot the result every few seconds so outcome becomes visible during calculations
        self.plotevery = self.cfg.get("plot.plotevery")
    
    def create_mesh_2d(self):
        self.nx_domain = int(self.domain_size[0] / self.dx)
        self.ny_domain = int(self.domain_size[1] / self.dy)
        self.mesh2d = Grid2D(dx = self.dx, dy = self.dy, nx = self.nx_domain, 
                        ny = self.ny_domain)
    
    def initial_boundary(self):
        self.initial_conc_DEET = self.cfg.get('initial.initial_conc_DEET')
        self.solution_DEET_void = CellVariable(name = "Active Component Conc in Void Space",
                                    mesh = self.mesh2d, value = self.initial_conc_DEET)
        self.viewer = None
        self.viewer = Viewer(vars = self.solution_DEET_void, datamin = 0., datamax = 0.0005)
        self.viewer.plot()
        self.viewerplotcount = 1
        
    def solve_yarn(self):
        for model in self.yarn_models:
            model.run()
        
    def solution(self):
        self.solution_DEET_void = self.solution_DEET_void
          
    def conc_out(self):
        self.eq = TransientTerm() == DiffusionTerm(coeff = D)
        x_center, y_center = self.mesh2d.getFaceCenters()
        boundary_zero_1 = (self.domain_size[0] - self.size_sample[0]) / 2.
        boundary_zero_2 = (self.domain_size[0] + self.size_sample[1]) / 2. 
        facesTop = ((self.mesh2d.getFacesLeft()) | (self.mesh2d.getFacesTop()) |
                        (self.mesh2d.getFacesRight()) | self.mesh2d.getFacesBottom()) 
        faceBottom = self.mesh2d.getFacesBottom() & (x_center >= boundary_zero_1 &
                        x_center <= boundary_zero_2)
        BCs = (FixedValue(faces = facesTop, value = self.boundary_up)) 
                #FixedValue(faces = faceBottom, value = self.release_source))
        self.steps = self.time_period / self.delta_t
        self.initial_t = 0
        for i in sp.arange(self.steps):
            self.initial_t += self.delta_t
            print 'the time step is', self.intial_t
            self.eq.solve(var = self.solution_DEET_void, boundaryConditions = BCs,
                        dt = self.delta_t)
            if self.viewer is not None:
                self.viewer.plot()
        raw_input("Finished <return>.....")
            
            
                
    def run(self):
        self.create_mesh_2d()
        self.solve_yarn()
        self.initial_boundary()
        self.conc_out()
    
        
        
        