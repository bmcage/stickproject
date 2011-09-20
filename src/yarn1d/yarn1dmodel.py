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
from scipy.integrate import vode
import matplotlib.pyplot as plt
import sets
import time
from copy import copy

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
import yarn2d.config as conf
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
        self.steps = int(self.time_period / self.delta_t)
        self.times = empty(self.steps+1, float)
        i=1
        self.times[0] = 0.
        while i <= self.steps:
            self.times[i] = self.times[i-1] + self.delta_t
            i += 1
            
        self.cfg_fiber = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            self.cfg_fiber.append(Fiber1dConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.time_period)
            #if self.cfg_fiber[-1].get("time.dt") > self.cfg.get("time.time_period"):
                #self.cfg_fiber[-1].set("time.dt", self.cfg.get("time.time_period"))
            self.cfg_fiber[-1].set("time.dt", self.delta_t)  
            
                       
        #create fiber models
        self.fiber_models = []
        for cfg in self.cfg_fiber:
            self.fiber_models.append(FiberModel(cfg))
        self.nr_models = len(self.fiber_models)
        self.verbose = self.cfg.get('general.verbose')
                
        #some memory
        self.cache_index_t_yarn = 0
        self.cache_index_t_fiber = [0] * self.nr_models
        
        #Initialize the tortuosity
        self.tortuosity= self.cfg.get('yarn.tortuosity')
        #constant diffusion coefficient
        self.diffusioncoeff = self.cfg.get('diffusion.diffusion_conc')
        self.boundary_transf_right = self.cfg.get('boundary.transfer_conc1')
        self.evap_equilibrium = self.cfg.get('boundary.evap_equilibrium')
        self.nr_fibers = self.cfg.get('fiber.number_fiber')
        self.blend=self.cfg.get('fiber.blend')
        

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
        self.init_conc = sp.ones(self.nr_edge-1, float)
        init_conc = self.cfg.get('initial.init_conc')
        self.init_conc *= init_conc
        self.conc = CellVariable(name = "Conc. Active Component", 
                   mesh = self.mesh_yarn, value = self.init_conc)
        self.viewer = None
        #self.viewer = Viewer(vars = self.conc, datamin = 0., datamax = None)

    def solve_fiber(self):
        """
        Solve the diffusion process for a repellent on the fiber at radial position r in the yarn.
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform.
        """        
        
        self.nr_timesteps = np.empty((self.nr_models),int)
        self.timesteps = [0]*self.nr_models
        self.fiber_surface = [0] * self.nr_models
        
        for ind, model in enumerate(self.fiber_models):
            print 'solving fibermodel', ind
            model.run()
            self.nr_timesteps[ind] = len(model.times)
            self.timesteps[ind] = copy(model.times)            
            self.fiber_surface[ind] = copy(model.fiber_surface)
            
            
    def _set_bound_flux(self, flux_edge, conc_r):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, conc_r contains solution in the cell centers
        """
        flux_edge[0]=0.
        #flux_edge[-1] = -self.boundary_transf_right * (self.evap_equilibrium-conc_r[-1])
        flux_edge[-1] = self.boundary_transf_right * conc_r[-1]
   
    def get_source(self, t):
        #find right index of interval for t in model.times and in self.times
        if self.times[self.cache_index_t_yarn] <= t and \
                t< self.times[self.cache_index_t_yarn+1] :
            self.index_t_yarn = self.cache_index_t_yarn
        else:
            self.index_t_yarn = None
            i = max([self.cache_index_t_yarn - 1,0])
            while i < self.steps-1:
                if self.times[i] <= t and t< self.times[i+1] :
                    self.index_t_yarn = i
                    break
                i += 1
            if self.index_t_yarn is None:
                #backward in time, so reducing timestep it seems
                i = self.cache_index_t_yarn-1
                while i > 0:
                    if self.times[i] <= t and t< self.times[i+1] :
                        self.index_t_yarn = i
                        break
                    i -= 1
            if self.index_t_yarn is None:
                #no interval found
                if t > self.times[-1]:
                    self.index_t_yarn = self.steps-1
                    print 'time over endtime', t, '>', self.times[-1], ", set index t to max", self.index_t_yarn
                else:
                    self.index_t_yarn = self.steps-1
                    print "endtime,", t, self.times, ", set index t to max", self.index_t_yarn
                    
                    #raise exception, 'something wrong'
        self.cache_index_t_yarn = self.index_t_yarn
        
        #the same now for the time of the fiber models
        nr = 0
        self.index_t_fiber = [0]*self.nr_models
        
        while nr < self.nr_models:
            if self.timesteps[nr][self.cache_index_t_fiber[nr]] <= t and \
                t< self.timesteps[nr][self.cache_index_t_fiber[nr]+1] :
                self.index_t_fiber[nr] = self.cache_index_t_fiber[nr]
                #print "interval found in loop 1"
            else:
                #self.index_t_fiber[nr] = None
                i = max([self.cache_index_t_fiber[nr] - 1,0])
                while i < self.nr_timesteps[nr] - 1:
                        if self.timesteps[nr][i] <= t and t< self.timesteps[nr][i+1] :
                            self.index_t_fiber[nr] = i
                            break
                        #print 'i', i, "index_t", self.index_t_fiber, "interval found in loop 2"
                        i += 1                                    
                if self.index_t_fiber[nr] is None:
                    #backward in time, so reducing timestep it seems
                    i = self.cache_index_t_fiber[nr]-1
                    while i > 0:
                        if self.timesteps[nr][i] <= t and t< self.timesteps[nr][i+1] :
                            self.index_t_fiber[nr] = i
                            break
                        #print "interval found in loop 3"
                        i -= 1                    
                if self.index_t_fiber[nr] is None:
                    #no interval found
                    if t > self.timesteps[nr][-1]:
                        print 'ERROR: time over endtime', t, '>', self.timesteps[nr][-1]
                        self.index_t_fiber[nr] = self.nr_timesteps[nr] - 1
                        break
                    else:
                        print nr, t, self.timesteps
                        raise Exception, 'something wrong'
            self.cache_index_t_fiber[nr] = self.index_t_fiber[nr]
            nr+=1
            
        #source term is n*Cf(R,r_i+,t)/2pi=(m*delta(r**2)_i/Ry**2)*Cf(R,r_i+,t)/2pi with n the number of fibers in a shell,
        #m the number of fibers per yarn.
        grid_square = np.power(self.grid_edge, 2)
        self.delta_rsquare = grid_square[1:] - grid_square[:-1]
        n = self.nr_fibers*self.delta_rsquare/(self.end_point**2)

        fibersurf = 0.
        for ind, blend in enumerate(self.blend):
            #print "size fiber_surf", size(self.fiber_surface), "size index_t", size(self.index_t_fiber), "ind", ind, "t", self.index_t_fiber[ind]
            fiber_surf_t = self.fiber_surface[ind][self.index_t_fiber[ind]]  +\
                    (self.fiber_surface[ind][self.index_t_fiber[ind]+1] 
                      - self.fiber_surface[ind][self.index_t_fiber[ind]])\
                    /(self.timesteps[ind][self.index_t_fiber[ind]+1]
                       -self.timesteps[ind][self.index_t_fiber[ind]])\
                    *(t-self.timesteps[ind][self.index_t_fiber[ind]])
            fibersurf = fibersurf + fiber_surf_t * blend/100
        self.source[self.index_t_yarn,:]=n*fibersurf/(2*np.pi)
        
    def f_conc1_ode(self, t,conc_r):
        return self.f_conc1(conc_r,t)
        
    def f_conc1(self, conc_r, t):
        print 't, concr', t, conc_r
        grid = self.grid
        n_cellcenters = len(grid)
        #get the sourceterm on time t
        self.get_source(t)
        #solve a fiber model on radial position r in the yarn
## TODO: update initial concentration for the fiber model on position r in yarn as the result of the last yarnsolution.
## for now: 1 fibermodel is solved. The same fibersurface result is used in every cell.        
        #self.solve_fiber()
        #self.endconc=sp.empty((self.steps,n_cellcenters),float)
        #for i in grid:
                #self.endconc[i][:]=self.fiber_surface[:]            
        #Initialize the flux rate on the edges
        flux_edge = sp.empty(n_cellcenters+1,float)
        self._set_bound_flux(flux_edge, conc_r)

        #calculate flux rate in each edge of the domain
        diff_u_t = sp.empty(n_cellcenters, float)
        flux_edge[1:-1] = 2*(self.diffusioncoeff/self.tortuosity)*(self.grid_edge[1:-1]*(conc_r[1:]-conc_r[:-1]))\
                          /(self.delta_r[:-1]+self.delta_r[1:])
        diff_u_t=((flux_edge[1:]-flux_edge[:-1])+ self.source[self.index_t_yarn,:])/ (self.grid_edge[:-1]*self.delta_r + np.power(self.delta_r,2)/2) \
                    
        print 'diff', diff_u_t
        #raw_input("cont")
        return diff_u_t    
    #def solve_odeint(self):
        #self.conc1=odeint(self.f_conc1, initial_c1, self.times)
        #self.view_sol(self.times, self.conc1)
        
    #def jac(self):
        #n=self.nr_edge-2
        #m=2*(self.nr_edge-2)+1
        #jaco = zeros((n,n))
        #i=1
        #j=0
        #while i <= m-1 and j<=n:
           # jaco[j,i-1]= self.grid_edge[(i+1)/2]/(self.delta_r[(i+1)/2]+self.delta_r[(i-1)/2])
           # jaco[j,i+1]= - self.grid_edge[(i+3)/2]/(self.delta_r[(i+3)/2]+self.delta_r[(i+1)/2]) - self.grid_edge[(i+1)/2]/(self.delta_r[(i+1)/2]+self.delta_r[(i-1)/2])
            #jaco[j,i+3]= self.grid_edge[(i+3)/2]/(self.delta_r[(i+3)/2]+self.delta_r[(i+1)/2])
            #i+=1
            #j+=1
            #print i,j
        #jaco = jaco*2*(self.diffusioncoeff/self.tortuosity)    
        #return jaco 
       
    def solve_ode(self):
        self.initial_t = 0.
        endT = self.time_period
        self.conc1 = np.empty((len(self.times), len(self.init_conc)), float)
        self.source = sp.empty((self.steps, self.nr_edge-1), float)
        #r = ode(self.f_conc1_ode,self.jac).set_integrator('vode', method = 'bdf', with_jacobian= True)
        r = ode(self.f_conc1_ode).set_integrator('vode', method = 'bdf', lband=1, uband=1)
        r.set_initial_value(self.init_conc, self.initial_t)
        tstep = 0
        self.conc1[tstep][:] = self.init_conc[:]
        while r.successful() and r.t < endT - self.delta_t /10.:
            r.integrate(r.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] =  r.y
            print "r.y", self.conc1[tstep][:]
        self.view_sol(self.times, self.conc1)
        raw_input("view solution")
  
    def view_sol(self, times, conc):
        """
        Show the solution in conc with times.
        conc[i][:] contains solution at time times[i]
        """
        self.solution_view = CellVariable(name = "yarn concentration", mesh = self.mesh_yarn, value = conc[0][:])
        self.viewer =  Viewer(vars = self.solution_view, datamin=0., datamax=conc.max()+0.2*conc.max())
        for time, con in zip(times[1:], conc[1:]):
            self.solution_view.setValue(con)
            self.viewer.plot()
            #if time == 200.0:
             #   dump.write({'space_position': self.grid, 'conc1': con},
              #              filename = utils.OUTPUTDIR + os.sep + 'ode_t2.gz', extension = '.gz')
    def run(self):        
        self.create_mesh()
        self.initial_yarn1d()
        self.solve_fiber()
        self.solve_ode()
        
