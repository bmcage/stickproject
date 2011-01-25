#
# Copyright (C) 2010  B. Malengier
# Copyright (C) 2010  P.Li
# Copyright (C) 2010 T. Goessens
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
Module holding a 1D cylindrical diffusion model for a yarn. 
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
import yarn2d.config as conf
from mycorrection import MyDiffusionTermNoCorrection
from yarn2dgrid import Yarn2dGrid
from fiberfipy.config import FiberfipyConfigManager
from fiberfipy.fibermodel import FiberModel

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
class Yarn1DModel(object):
    """
    Yarn1DModel is a special diffusion model for a single yarn which is composed 
    by a certain amount of fibers. A cross-section of a fiber is 
    generated. The domain is a line from the center of the yarn to the surface. 
    On this line there are some fibers distributed as in the module yarn1Dgrid.
    Only diffusion processes in a single fiber and yarn are considered. 
    ODE of scipy solve the diffusion process in the layers of DEET and permithrine
    which are on the fiber
    Fipy solve the transient diffusion problem in the whole domain
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        """
        self.datatime = []
        self.cfg = config
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = self.time_period / self.delta_t
        self.cfg_fiber = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            self.cfg_fiber.append(FiberfipyConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.cfg.get("time.time_period"))
        
        #create fiber models
        self.fiber_models = []
        for cfg in self.cfg_fiber:
            self.fiber_models.append(FiberModel(cfg))
        self.verbose = self.cfg.get('general.verbose')
        self.diffusion_coeff = sp.empty(self.tot_edges-1, float)

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh!
        
        grid: the space position of each central point for every cell element (r-coordinate);
        """
        self.beginning_point = 0 #center of the yarn, r=0, with r the distance from center yarn.
        self.end_point = self.cfg.get('domain.yarnradius')
        self.nr_edge = 10
        self.diff_coef = [0.]
        self.init_conc = lambda x: 0.0
        #we now construct the full edge grid
        self.grid_edge = sp.linspace(self.beginning_point, self.end_point, self.nr_edge)
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]
        #create cylindrical 1D grid over domain.
        self.mesh_yarn = CylindricalGrid1D(dr=tuple(self.delta_r))
        self.mesh_yarn.periodicBC = False
        self.mesh_yarn = self.mesh_yarn + (self.beginning_point,)
                     
    def initial_yarn1d(self):
        """ initial concentration over the domain"""
        self.init_conc = self.cfg.get('initial.init_conc')
        self.conc = CellVariable(name = "", 
                    mesh = self.mesh_yarn, value = self.init_conc)
        #self.viewer = None
        self.viewer = Viewer(vars = self.conc, datamin = 0., datamax = none)
#Viewer(vars = self.conc, datamin = 0., datamax =0.005)
        self.initial_c1 = sp.empty(self.tot_edges-1, float)
        st = 0
        surf = self.surf[st]
        for i, pos in enumerate(self.grid):
                while pos > surf:
                    st += 1
                    surf = self.surf[st]
                self.initial_c1[i] = self.init_conc[st](pos)
                self.diffusion_coeff[i] = self.diff_coef[st]

    def solve_fiber(self):
        """
        Solve the diffusion process for a repellent on the fiber at radial position r in the yarn.
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform.
        
        """
        for model in self.fiber_models:
            model.run()
            
   
    def _set_bound_flux(self, flux_edge, conc_r):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, conc_r contains solution in the cell centers
        """
        self.boundary_transf_right = self.cfg.get('boundary.transfer_conc1')
        flux_edge[-1] = -self.boundary_transf_right * conc_r[-1]
            
            
    def f_conc1(self, conc_r):
        grid = self.grid
        n_cellcenters = len(grid)
        #Initialize the left side of ODE equations
        diff_u_t = sp.empty(n_cellcenters, float)
        #solve a fiber model on radial position r in the yarn
## TODO: update initial concentration for the fiber model on position r in yarn as the result of the last yarnsolution.
## for now: 1 fibermodel is solved. The same fibersurface result is used in every cell.        
        self.solve_fiber()
        self.endconc=sp.empty((self.steps,n_cellcenters),float)
        for i in grid:
                self.endconc[i][:]=self.fiber_surface[:]
                
        #Initialize the flux rate on the edges
        flux_edge = sp.empty(n_cellcenters+1, float)
        source=sp.empty(self.steps,float)
        self._set_bound_flux(flux_edge, conc_r)
        #Initialize the tortuosity
        self.tortuosity= self.cfg.get('tortuosity')
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = (self.diffusion_coeff[:-1] * sp.exp(-self.diff_exp_fact * conc_r[:-1]) \
                         + self.diffusion_coeff[1:] * sp.exp(-self.diff_exp_fact * conc_r[1:]))\
                    * (self.grid_edge[1:-1]/self.tortuosity) \
                    * (conc_r[1:] - conc_r[:-1])\
                    / ((self.delta_r[:-1] + self.delta_r[1:]))
        for i in grid:
            source[:]=self.endconc[i][:-1]+self.endconc[i][1:]/(4*math.pi)
        diff_u_t[:]=(flux_edge[1:]-flux_edge[:-1])/(2*self.grid_edge[:]*self.delta_r[:]+self.delta_r[:]**2)+source[:]
        return diff_u_t
    
    def f_conc1_ode(self, conc_r):
        return self.f_conc1(conc_r)
    
    def solve_odeint(self):
        self.conc1=odeint(self.f_conc1, initial_c1, self.times)
        self.view_sol(self.times, self.conc1)
    
    
    def solve_ode(self):
        self.delta_t = self.times[1]-self.times[0]
        self.initial_t = self.times[0]
        endT = self.times[-1]
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        r = ode(self.f_conc1_ode).set_integrator('vode', method = 'bdf')
        r.set_initial_value(initial_c1, self.initial_t)#.set_f_params(2.0)
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        while r.successful() and r.t < endT - self.delta_t /10.:
            r.integrate(r.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] = r.y 
            self.view_sol(self.times, self.conc1)
  
    def run(self):        
        self.create_mesh()
        self.initial_yarn2d()
        self.solve_fiber()
        self.solve_ode()
        
