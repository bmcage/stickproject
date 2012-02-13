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
from scipy import integrate
#from scipy.integrate import ode
#from scipy.integrate import vode
import matplotlib.pyplot as plt
import sets
import time
import mpmath

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
import bednet.config as conf
#from mycorrection import MyDiffusionTermNoCorrection
#from yarn2dgrid import Yarn2dGrid
#from numMesh import Grid2D
#from yarn2dfiber import Yarn2dFiber
from yarn1d.config import *
from yarn1d.yarn1dmodel import *
#from fiberfipy.config import FiberfipyConfigManager
#from fiberfipy.fibermodel import FiberModel

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
    a volatile outside the fabric. This upscaling method is simple with the equation:
    C_outside = C_yarn * R_yarn/d_distance
    """
    def __init__(self, config):
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.timesteps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        #set correct delta_t
        self.delta_t = self.times[1]-self.times[0]
        if self.verbose:
            print "Timestep used in bednet model:", self.delta_t

        #self.n = self.cfg.get('domain.nr_vert_yarns')
        #self.m = self.cfg.get('domain.nr_hor_yarns')   
        self.domain_size = self.cfg.get('domain.domain_size')
        self.dx = self.cfg.get('domain.dx')
        self.dy  = self.cfg.get('domain.dy')
        self.size_sample = self.cfg.get('sample.size_sample')
        #self.thickness_sample = self.cfg.get('sample.thickness_sample')
        self.boundary_up = self.cfg.get('boundary.boundary_up')
        self.boundary_bottom = self.cfg.get('boundary.boundary_bottom')
        self.boundary_left = self.cfg.get('boundary.boundary_left')
        self.boundary_right = self.cfg.get('boundary.boundary_right')
        self.diff_coef = self.cfg.get('diffusion.diff_coef')
        self.saturation_conc = self.cfg.get('saturation.saturation_conc')
        self.x0 = self.cfg.get('observer.x0')
        """
        self.distance_yarn = self.cfg.get('domain.distance_yarn')
        self.grid_space_vertical = self.cfg.get('domain.grid_space_vertical')
        self.number_nodes = (self.domain_size[0] / self.distance_yarn + 1) * \
                            (self.domain_size[1] / self.grid_space_vertical + 1)
        """
        self.cfg_yarn = []        
        for filename in self.cfg.get('sample.yarn_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                        os.path.dirname(self.cfg.filename), filename))
            self.cfg_yarn.append(YarnConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            print 'time', self.time_period
            self.cfg_yarn[-1].set("time.time_period", self.time_period)
            #self.cfg_yarn[-1].set("time.dt", self.delta_t)  
            
        #create yarn models
        self.yarn_models = []
        for cfg in self.cfg_yarn:
            self.yarn_models.append(Yarn1DModel(cfg))
        self.nr_models = len(self.yarn_models)

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
        ##xfc, yfc = self.mesh2d.getFaceCenters()
        ##xcc, ycc = self.mesh2d.getCellCenters()
        ##self.cell_volume = self.mesh2d.getCellVolumes()
    
    def initial_void_conc(self):
        """ initial concentration over the domain"""
        self.init_void = self.cfg.get('initial.init_void')
        if self.plotevery:
            self.solution_void = CellVariable(name = "Active Component Conc in Void Space",
                                        mesh = self.mesh2d, value = self.init_void)
            self.viewer = None
            self.viewer = Viewer(vars = self.solution_void, datamin = 0., datamax = 0.0005)
            self.viewer.plot()
            self.viewerplotcount = 1
        
    def init_yarn(self):
        self.yarn_surface = [0] * self.nr_models
        self.yarn_surface_time = [0] * self.nr_models
        
        for ind, model in enumerate(self.yarn_models):
            self.yarn_surface[ind] = []
            self.yarn_surface_time[ind] = []
            model.create_mesh()
            model.initial_yarn1d()
            model.solve_fiber_init()
            model.solve_ode_init()
    
    def initial_boundary_conc(self):
        self.initconc = self.cfg.get('initial.init_conc')
        self.initvoidconc = self.cfg.get('initial.init_void')

    def solve_timestep(self, t):
        
        for ind, model in enumerate(self.yarn_models):
            rt, rety = model.do_yarn_step(t)
            self.yarn_surface[ind].append(rety[-1])
            self.yarn_surface_time[ind].append(t)
        
        # we need to determine value conc volatile at x0 where observer is
        # then, for next timestep, we need to set correct boundary condition
        # on the yarn level
        x0 = self.x0
        Heaviside = lambda x: (sp.sign(x)+1)/2
        print "TODO, CHECK FOLLOWING, WHAT is delta_t??"
        # concentration is a consequence of all previous releases, so sum 
        # over all times, and compute contribution of that moment.
        C = self.diff_coef*(t-self.index_t_bednet*delta_t)
        AH = math.pow(dx,2)/4*C
        AW = math.pow(dy,2)/4*C
        conc = yarn_surface[:][self.index_t_yarn] * (Heaviside(t-self.index_t_bednet*delta_t)/(8 * math.pi * C)) *exp(-math.pow(x0,2)/4*C)* \
                    (mpmath.elliptic.jtheta(3,0,exp(-AH))+mpmath.elliptic.jtheta(3,0,exp(-AW))+2)
                    
        self.solution_void = conc

    def run(self):
        self.create_mesh_2d()
        self.init_yarn()
        self.initial_boundary()
        for t in self.times[1:]:
            self.solve_timestep(t)
