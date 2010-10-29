#
# Copyright (C) 2010  P.Li, B. Malengier
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

#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
import os,sys
import math,cmath
import scipy as sp
import pylab as pl
from scipy.integrate import odeint, ode

class Solving1DFiber(object):
    """
    Solving the 1D convection-diffusion equation on the single fiber with conser
    -vative forum-Finite Volume Scheme:
     c_t=flux_right-flux_right
    In current model, all the cells are full cells and the boundary condition is
    Neumann Condition (Constant Flux rate on the boundary)
    """
    def __init__(self, grid, initial_c1, boundary_fib_left, boundary_fib_right, diffusion_co_l1, 
                diffusion_co_l2):
        """
        grid: the space position of each central point for every cell element;
        boundary_left: the boundary condition for left side of 1D domain;
        boundary_right: the boundary condition for right side of 1D domain;
        diffusion_co_l1: the diffusion coefficient of DEET in the first layer of binder;
        diffusion_co_l2: the diffusion coefficient of DEET in the second layer of binder;
        initial_c: the initial concentration of DEET in the whole domain
        """
        self.grid = grid
        self.initial_c1 = initial_c1
        self.boundary_fib_left = boundary_fib_left
        self.boundary_fib_right = boundary_fib_right
        self.diffusion_co_l1 = diffusion_co_l1
        self.diffusion_co_l2 = diffusion_co_l2
        
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
        diffusion_co[:] = position_diffusion[:] * sp.exp(self.initial_c1[:]) * grid[:]
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = ((diffusion_co[:-1] + diffusion_co[1:]) / 2.)*\
                          (w_rep[1:] - w_rep[:-1])/delta_r
        diff_w_t[:]=(flux_edge[1:]-flux_edge[:-1])/delta_r
        return diff_w_t
    
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
        flux_edge[1:-1] = self.diffusion_co_l1*\
                          (w_rep[1:]/grid[1:] - w_rep[:-1]/grid[:-1])/delta_r
        diff_w_t[:]=(flux_edge[1:]-flux_edge[:-1])/delta_r
        return diff_w_t
    
    def solver_c(self, times):
        self.times = times
        initial_w = self. initial_c1 * self.grid
        self.solv=odeint(self.f_conc1, initial_w, times)
        print 'This is solution', self.solv
        self.conc1=self.solv/ self.grid
        print 'is this true', self.conc1[-1] == self.solv[-1]/self.grid
    
    def solver_c_2(self, delta_t, initial_t):
        r = ode(self.f_conc2).set_integrator('vode', method = 'bdf')
        self.delta_t = delta_t
        self.initial_t = initial_t
        initial_w1 = self.initial_c1 * self.grid
        r.set_initial_value(initial_w1, self.initial_t)#.set_f_params(2.0)
        while r.successful() and r.t < self.initial_t + self.delta_t:
            r.integrate(r.t + self.delta_t)
            self.conc1 = r.y / self.grid    


        