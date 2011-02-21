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
from scipy.integrate import ode
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
        self.times=empty(self.steps,float)
        i=0
        while i<=self.steps:
            self.times[i]+=self.delta_t
            
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
        #print 'delta_r[:-1]', self.delta_r[:-1].shape
        #print 'delta_r[1:]', self.delta_r[1:].shape
        #print 'som van delta_r', self.delta_r[1:]+self.delta_r[:-1]
                     
    def initial_yarn1d(self):
        """ initial concentration over the domain"""
        init_conc = self.cfg.get('initial.init_conc')
        self.init_conc = sp.ones(self.nr_edge-1, float)
        self.init_conc *= init_conc
        self.conc = CellVariable(name = "", 
                    mesh = self.mesh_yarn, value = self.init_conc)
        #self.viewer = None
        self.viewer = Viewer(vars = self.conc, datamin = 0., datamax = None)
#Viewer(vars = self.conc, datamin = 0., datamax =0.005)
        #self.initial_c1 = sp.empty(self.tot_edges-1, float)
        #st = 0
        #surf = self.surf[st]
        #for i, pos in enumerate(self.grid):
         #       while pos > surf:
          #          st += 1
           #         surf = self.surf[st]
            #    self.initial_c1[i] = self.init_conc[st](pos)
             #   self.diffusion_coeff[i] = self.diff_coef[st]

    def solve_fiber(self):
        """
        Solve the diffusion process for a repellent on the fiber at radial position r in the yarn.
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform.
        
        """
        for ind,model in enumerate(self.fiber_models):
            model.run()
            self.nr_models = len(self.fiber_models)
            self.nr_timesteps = len(model.times)
            self.timesteps= sp.empty((self.nr_models,self.nr_timesteps),float)
            self.timesteps[ind][:]=model.times            
            self.fiber_surface=sp.empty((self.nr_models,self.nr_timesteps),float)
            self.fiber_surface[ind][:] = model.fiber_surface          

    def _set_bound_flux(self, flux_edge, conc_r):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, conc_r contains solution in the cell centers
        """
        self.boundary_transf_right = self.cfg.get('boundary.transfer_conc1')
        flux_edge[-1] = -self.boundary_transf_right * conc_r[-1]
    
    def get_source(self,t):
        #find right index of t in model.times and in self.times
        i=0
        while i <= self.steps:
            if self.times[i] <= t and t< self.times[i+1] :
               self.index_t_yarn=i
            print i
            i+=1
            
        nr=0           
        j=0  
        while nr <= self.nr_models:       
            while j <= len(self.timesteps[nr][:]):
                if self.timesteps[nr][j] <= t and t < self.timesteps[nr][j+1]:
                    self.index_t_fiber=j
                print j
                j+=1
            print nr
            nr+=1
            
        source=sp.empty(self.steps,float)        
        #source term is n*Cf(R,r_i+,t)/2pi=(m*delta(r**2)_i/Ry)*Cf(R,r_i+,t)/2pi with n the number of fibers in a shell,
        #m the number of fibers per yarn.
        self.nr_fibers = self.cfg.get('fiber.number_fiber')
        grid_square=self.grid_edge**2
        self.delta_rsquare=grid_square[1:]-grid_square[:-1]
        n = self.nr_fibers*self.delta_rsquare/(self.end_point**2)
        self.blend=self.cfg.get('fiber.blend')
        self.interpolated_fibersurf=sp.empty(self.steps,float)
        for i in enumerate(self.blend):
            self.interpolated_fibersurf[self.index_t_yarn]+=self.fiber_surface[i][self.index_t_fiber]*self.blend[i]/100
        source[self.index_t_yarn]=n*self.interpolated_fibersurf[self.index_t_yarn]/(2*math.pi)
        
        
    def f_conc1(self, conc_r, t):
        print 'concr', t, conc_r
        grid = self.grid
        n_cellcenters = len(grid)
        #get the sourceterm on time t
        source=self.get_source(t);
        #solve a fiber model on radial position r in the yarn
## TODO: update initial concentration for the fiber model on position r in yarn as the result of the last yarnsolution.
## for now: 1 fibermodel is solved. The same fibersurface result is used in every cell.        
        self.solve_fiber()
        #self.endconc=sp.empty((self.steps,n_cellcenters),float)
        #for i in grid:
                #self.endconc[i][:]=self.fiber_surface[:]            
        #Initialize the flux rate on the edges
        flux_edge = sp.empty(n_cellcenters+1,float)
        self._set_bound_flux(flux_edge, conc_r)
        #Initialize the tortuosity
        self.tortuosity= self.cfg.get('yarn.tortuosity')
        #constant diffusion coefficient
        self.diffusioncoeff = self.cfg.get('diffusion.diffusion_conc')
        
        #calculate flux rate in each edge of the domain
        flux_edge[0]=0
        concdiff=conc_r[1:]-conc_r[:-1]
        deel1=self.grid_edge[1:-1]*concdiff
        flux_edge[1:-1] = (2*self.diffusioncoeff*deel1)\
                          /((self.delta_r[:-1]+self.delta_r[1:])*self.tortuosity)
        diff_u_t=(flux_edge[1:]-flux_edge[:-1])/(2*self.grid_edge[:]*self.delta_r[:]+self.delta_r[:]**2)+source
        return diff_u_t
    
    def f_conc1_ode(self, t,conc_r):
        return self.f_conc1(conc_r,t)
    
    #def solve_odeint(self):
        #self.conc1=odeint(self.f_conc1, initial_c1, self.times)
        #self.view_sol(self.times, self.conc1)
    
    
    def solve_ode(self):
        self.initial_t = 0
        endT = self.time_period
        self.conc1 = np.empty((self.steps, self.nr_edge-1), float)
        r = ode(self.f_conc1_ode).set_integrator('vode', method = 'bdf')
        r.set_initial_value(self.init_conc, self.initial_t)#.set_f_params(2.0)
        tstep = 0
        self.conc1[tstep][:] = self.init_conc[:]
        while r.successful() and r.t < endT - self.delta_t /10.:
            r.integrate(r.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] = r.y 
            self.view_sol(self.times, self.conc1)
  
    def run(self):        
        self.create_mesh()
        self.initial_yarn1d()
        self.solve_fiber()
        self.solve_ode()
        
