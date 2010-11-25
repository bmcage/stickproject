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
from scipy.integrate import odeint, ode
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

#-------------------------------------------------------------------------
#
#Fipy Imports
#
#-------------------------------------------------------------------------
from fipy import *

#-------------------------------------------------------------------------
#
# DiffusionModel class 
#
#-------------------------------------------------------------------------
class FiberModel(object):
    """
    fiberfipy.FiberModel is a special diffusion model for a single radial 
    symmetric fiber which is composed of a specific material (cotton, polyester,...)
    and has a number of coatings.
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        
        boundary_left: the boundary condition for left side of 1D domain;
        boundary_right: the boundary condition for right side of 1D domain;
        diffusion_co_l1: the diffusion coefficient of DEET in the first layer of binder;
        diffusion_co_l2: the diffusion coefficient of DEET in the second layer of binder;
        initial_c: the initial concentration of DEET in the whole domain
        """
        self.datatime = []
        self.cfg = config
        self.comps = self.cfg.get('general.components')
        self.method = self.cfg.get('general.method')
        if not (self.method in ['FVM']):
            print 'unkown solution method'
            sys.exit(0)
        self.submethod = self.cfg.get('general.submethod')
        if not (self.submethod in ['ode', 'odeint', 'fipy']):
            print 'unkown solution submethod'
            sys.exit(0)
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = self.time_period / self.delta_t
        #read the initial and boundary information for fiber
        self.n_point = self.cfg.get('fiber.n_point') #discretize the fiber radius
        self.boundary_fib_left = self.cfg.get('boundary.boundary_fib_left')
        self.boundary_fib_right = self.cfg.get('boundary.boundary_fib_right')
        self.diffusion_co_l1 =  self.cfg.get('diffusion.diffusion_co_l1')
        self.diffusion_co_l2 = self.cfg.get('diffusion.diffusion_co_l2')
        self.init_conc1_fiber = eval(self.cfg.get('initial.init_conc1_fiber'))
        self.transfer_conc1 = self.cfg.get('transfer.transfer_conc1')
        
        self.verbose = self.cfg.get('general.verbose')

    def create_mesh(self, scaleL=1.):
        """
        Create a mesh for use in the model
        scaleL : scale factor to use in the length scale
        
        grid: the space position of each central point for every cell element;
        """
        self.scaleL = scaleL
        self.beginning_point = self.cfg.get('fiber.beginning_point')
        self.end_point = self.cfg.get('fiber.end_point')
        scale_beginning = self.beginning_point * self.scaleL
        print 'this is the scale beginning point:',scale_beginning
        scale_end = self.end_point * self.scaleL
        print 'this is the scale end point', scale_end
        self.grid = sp.linspace(scale_beginning, scale_end, self.n_point)
        if self.submethod == 'fipy':
            self.delta_r = self.grid[1] - self.grid[0]
            self.mesh_fiber = CylindricalGrid1D(nr = self.n_point, dr = self.delta_r)

        self.diffusion_coeff = sp.empty(self.n_point, float)
    
    def initial_fiber(self):
        """ initial concentration over the domain"""
        self.initial_c1 = sp.empty(self.n_point, float)
        for i in sp.arange(0, self.n_point, 1):
            ##TODO, this is not good, change is in middle
            if i <= (self.n_point - 1) / 2:
                self.initial_c1[i] = self.init_conc1_fiber(i)[0] 
                self.diffusion_coeff[i] = self.diffusion_co_l1
            elif i > (self.n_point - 1) /2:
                self.initial_c1[i] = self.init_conc1_fiber(i)[1] 
                self.diffusion_coeff[i] = self.diffusion_co_l2
        if self.submethod == 'fipy':
            #radial symmetry, solve for w with w = r*C
            self.initial_c2 = self.initial_c1 * self.grid
        print 'initial condition is:', self.initial_c1
        
        #self.conc1 = CellVariable(name = "solution concentration1", 
        #            mesh = self.mesh2d, value = self.init_conc1)
        #self.viewer = None
        #self.viewer = Viewer(vars = self.conc1, datamin = 0., datamax =0.005)

    def f_conc1_ode(self, t, w_rep):
        return self.f_conc1(w_rep, t)

    def f_conc1(self, w_rep, t):
        grid = self.grid
        n_cell = len(grid)
        delta_r=grid[1]-grid[0]
        #Initialize the left side of ODE equations
        diff_w_t = sp.zeros(n_cell, float)
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = sp.zeros(n_cell+1, float)
        flux_edge[0] = self.boundary_fib_left
        flux_edge[-1] = self.boundary_fib_right
        #Diffusion coefficient changes with the concentration changing
        position_diffusion = sp.empty(n_cell, float)
        position_diffusion[:] = self.diffusion_co_l1 + (self.diffusion_co_l2 - self.diffusion_co_l1)/(1 + \
                 sp.exp(-100*(grid[:]-grid[(n_cell-1)/2])))
        diffusion_co = sp.empty(n_cell, float)
        ## TODO, this is wrong !!
        diffusion_co[:] = position_diffusion[:] * sp.exp(self.initial_c1[:]) * grid[:]
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = ((diffusion_co[:-1] + diffusion_co[1:]) / 2.)*\
                          (w_rep[1:] - w_rep[:-1])/delta_r
        diff_w_t[:]=(flux_edge[1:]-flux_edge[:-1])/delta_r
        return diff_w_t


    def f_conc2_ode(self, t, w_rep):
        return self.f_conc2(w_rep, t)

    def f_conc2(self, t, w_rep):
        grid = self.grid
        n_cell = len(grid)
        delta_r=grid[1]-grid[0]
        #Initialize the left side of ODE equations
        diff_w_t = sp.zeros(n_cell, float)
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = sp.zeros(n_cell+1, float)
        flux_edge[0] = self.boundary_fib_left
        flux_edge[-1] = self.boundary_fib_right
        #calculate flux rate in each edge of the domain
        ## TODO, this is wrong, should be r * ? !!
        flux_edge[1:-1] = self.diffusion_co_l1*\
                          (w_rep[1:]/grid[1:] - w_rep[:-1]/grid[:-1])/delta_r
        diff_w_t[:]=(flux_edge[1:]-flux_edge[:-1])/delta_r
        return diff_w_t
    
    def solve_odeint(self):
        initial_w = self. initial_c1 * self.grid
        self.solv=odeint(self.f_conc1, initial_w, times)
        print 'This is solution', self.solv
        self.conc1=self.solv/ self.grid
        print 'is this true', self.conc1[-1] == self.solv[-1]/self.grid
    
    def solve_ode(self):
        self.delta_t = self.times[1]-self.times[0]
        self.initial_t = self.times[0]
        endT = self.times[-1]
        self.conc1 = np.empty((len(times), len(self.initial_c1)), float)
        r = ode(self.f_conc1_ode).set_integrator('vode', method = 'bdf')
        initial_w1 = self.initial_c1 * self.grid
        r.set_initial_value(initial_w1, self.initial_t)#.set_f_params(2.0)
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        while r.successful() and r.t < endT - self.delta_t /10.:
            r.integrate(r.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] = r.y / self.grid

    def solve_fipy(self):
        #using fipy to solve 1D problem in fiber
        self.solution_fiber = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber,
                            value = self.initial_c1, hasOld = 1)
        self.viewer =  Viewer(vars = self.solution_fiber, datamin=0., datamax=1.1)
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        print 'the length of solution:', len(self.solution_fiber)
        print 'length of the diffusion coefficient', len(self.diffusion_coeff)
        print len(self.grid)
        
        self.BCs_fiber = (FixedFlux(faces = self.mesh_fiber.getFacesRight(), 
                                    value = self.boundary_fib_right),
                          FixedFlux(faces = self.mesh_fiber.getFacesLeft(), 
                                    value = 0.0))
        self.eqX_fiber = TransientTerm() == DiffusionTerm(coeff = 
                        self.diffusion_coeff * sp.exp(-self.solution_fiber))
        """
        self.eqX_fiber = TransientTerm() == DiffusionTerm(coeff = (self.diffusion_co_l1 + (self.diffusion_co_l2 - \
                                        self.diffusion_co_l1)/(1 + sp.exp(-100 * (self.grid - (self.beginning_point * self.grid.scaleL + \
                                        self.grid[(self.n_point - 1)]))))) * sp.exp(-solution_fiber))
        """
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        for time in self.times[1:]:
            self.solve_fipy_step()
            if self.viewer is not None:
                self.viewer.plot()
            tstep += 1
            print time
            self.conc1[tstep][:] = self.solution_fiber.getValue()

    def solve_fipy_step(self):
        res = 1e+1
        while res > 1e-8:
            res = self.eqX_fiber.sweep(var = self.solution_fiber, 
                                        boundaryConditions = self.BCs_fiber,
                                        dt = self.delta_t)
        self.solution_fiber.updateOld()

    def solve(self):
        """
        Solve the diffusion process in the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        discretization_t = self.steps + 1
        self.times = sp.linspace(0, self.time_period, discretization_t)
        if self.submethod == 'fipy':
            self.solve_fipy()
        else:            
            if self.submethod == 'odeint':
                self.solve_odeint()
            elif  self.submethod == 'ode':
                self.solve_ode()
        print self.conc1
        self.fiber_surface = sp.empty(len(self.times), float)
        for i in sp.arange(1,len(self.times) + 1,1):
            self.fiber_surface[i - 1] = self.conc1[i - 1][-1]
        print self.fiber_surface[:]
    
    def run(self):        
        self.create_mesh()
        self.initial_fiber()
        self.solve()
        
