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
import sets
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
import test.config as conf
from mycorrection import MyDiffusionTermNoCorrection
#from yarn2dgrid import Yarn2dGrid
#from yarn2dgridnew import Yarn2dNewGrid
#from yarn2d_overlap import Yarn2DOverlapping
from fiber1d.config import Fiber1dConfigManager
from fiber1d.fibermodel import FiberModel

#-------------------------------------------------------------------------
#
#Fipy Imports
#-------------------------------------------------------------------------
from fipy import *

#-------------------------------------------------------------------------
#
# DiffusionModel class 
#
#-------------------------------------------------------------------------

class TestFiber(object):
    """
    To test the solver of DEET diffusion in fiber, especailly the method which
    uses the yarn-fiber time period. 
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        """
        self.datatime = []
        self.cfg = config
        self.time_period = self.cfg.get('time.time_period')
        self.step = self.cfg.get('time.step')
        self.steps = (self.time_period*(1.+self.step*1e-6)) // self.step #self.delta_t
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.delta_t = self.cfg.get('time.dt')
##        self.steps = (self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t
##        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.delta_t = self.times[1] - self.times[0]
        self.Ry = self.cfg.get('domain.yarnradius')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        self.Rf = self.cfg.get('fiber.radius_fiber')
        self.radius_fiber =  [self.scaleL * rad for rad in self.Rf]
        self.eps_value = self.cfg.get('fiber.eps_value')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        self.nrtypefiber = self.cfg.get('fiber.number_type')
        assert self.nrtypefiber == len(self.blend) == len(self.Rf)
        #computational radius
        self.cfg_fiber = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            self.cfg_fiber.append(Fiber1dConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.cfg.get("time.time_period"))
            if self.cfg_fiber[-1].get("time.dt") > self.cfg.get("time.time_period"):
                self.cfg_fiber[-1].set("time.dt", self.cfg.get("time.time_period"))
        
        #create fiber models
        self.fiber_models = []
        for cfg in self.cfg_fiber:
            self.fiber_models.append(FiberModel(cfg))

        self.verbose = self.cfg.get('general.verbose')
        
    def solve_fiber_step(self):
        """
        Solve the diffusion process on the fiber for a single step, starting
        from where we where last. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        step = 0.1
        for model in self.fiber_models:
            model.run_init()
            print 'self.times', self.times
            raw_input("continue <Enter>")
            for i in self.times:
            #for model in self.fiber_models:
                model.run_step(step)
                print 'times', len(self.times)
            #times of the models must coincide at the moment as later on we do
            #conc_on_fib[nyfib] = self.fiber_models[nyfib].fiber_surface[i+1]
            #we should interpolate to avoid that
##            assert len(self.times) == len(model.times), 'fiber and yarn model need same time steps'
##            test = (sp.array(self.times) == sp.array(model.times))
##            assert test.all(), 'fiber and yarn model need same time steps'
            
    def solve_fiber(self):
        """
        Solve the diffusion process on the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        for model in self.fiber_models:
            model.run()
            #times of the models must coincide at the moment as later on we do
            #conc_on_fib[nyfib] = self.fiber_models[nyfib].fiber_surface[i+1]
            #we should interpolate to avoid that
            assert len(self.times) == len(model.times), 'fiber and yarn model need same time steps'
            test = (sp.array(self.times) == sp.array(model.times))
            assert test.all(), 'fiber and yarn model need same time steps'

    def solve_single_component(self):
        """
        The DEET diffusion process is divided into two parts:
        (1) DEET diffuses through the layer containing permithrine on the fiber 
        and reaches the surface;
        (2) DEET begins to diffuse in the void space of  yarn
        So it means that the boundary condition of fiber has two steps:
        (1) When the DEET does not reach the surface of fiber, the inner and out
        boundaries are no flux boundary condition;
        (2) When the DEET reaches surface, the evaporation happens. So the boundaries 
        of fiber and yarn are changed to constant flux (Neumann boundary condition)
        """
        self.solve_fiber_step()

    def run(self):
        self.solve_single_component()