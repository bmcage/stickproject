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
from scipy.integrate import odeint, ode, trapz
import matplotlib.pyplot as plt
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
from fiber1d.config import METHOD, FLUX, TRANSFER, BOUND_TYPE

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
    fiber1d.FiberModel is a special diffusion model for a single radial 
    symmetric fiber which is composed of a specific material (cotton, polyester,...)
    and has a number of coatings.
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        
        boundary_left: the Neumann boundary condition for left side of 1D domain;
        boundary_right: the Neumann boundary condition for right side of 1D domain;
        transfer: the Robin BC for the right side of 1D domain
        initial_c: the initial concentration of DEET in the whole domain
        """
        self.datatime = []
        self.cfg = config
        self.method = self.cfg.get('general.method')
        if not (self.method in METHOD):
            print 'ERROR: unkown solution method %s' % self.method
            sys.exit(0)
        self.submethod = self.cfg.get('general.submethod')
        if not (self.submethod in METHOD[self.method][1]):
            print 'ERROR: unkown solution submethod %s' % self.submethod
            sys.exit(0)
        self.fiber_diff = self.cfg.get('fiber.internal_diffusion')
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = (self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.delta_t = self.times[1] - self.times[0]
        #read the initial and boundary information for fiber
        self.n_edge = self.cfg.get('fiber.n_edge') #discretize the fiber radius
        self.bound_left = BOUND_TYPE[self.cfg.get('boundary.type_left')]
        self.bound_right = BOUND_TYPE[self.cfg.get('boundary.type_right')]
        self.boundary_fib_left = self.cfg.get('boundary.boundary_fib_left')
        self.boundary_fib_right = self.cfg.get('boundary.boundary_fib_right')
        self.boundary_transf_right = self.cfg.get('boundary.transfer_right')
        
        self.verbose = self.cfg.get('general.verbose')

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh
        grid: the space position of each central point for every cell element;
        """
        if not self.fiber_diff:
            n_edge = [0]
            self.diff_coef = [0.]
            self.diff_exp_fact = [0.]
            self.init_conc = [lambda x: 0.0 ]
            self.ind_first_zone = 1
        else:
            n_edge = [self.cfg.get('fiber.n_edge')]
            self.diff_coef = [self.cfg.get('fiber.diffusion_coef')]
            self.diff_exp_fact = [self.cfg.get('fiber.diffusion_polymer_exp_factor')]
            self.init_conc = [eval(self.cfg.get('fiber.init_conc')) ]
            self.ind_first_zone = 0
        self.nrlayers = self.cfg.get('fiber.nrlayers')
        self.surf_begin = [0.]
        self.surf = [self.cfg.get('fiber.radius_pure_fiber')]
        for i in range(self.nrlayers):
            section = 'fiberlayer_%i' % i
            n_edge += [self.cfg.get(section + '.n_edge')]
            self.surf_begin += [self.surf[-1]]
            self.surf += [self.surf[-1] + self.cfg.get(section + '.thickness')]
            self.diff_coef += [self.cfg.get(section + '.diffusion_coef')]
            self.diff_exp_fact += [self.cfg.get(section + '.diffusion_polymer_exp_factor')]
            self.init_conc += [eval(self.cfg.get(section + '.init_conc')) ]

        #we now construct the full edge grid
        self.tot_edges = 0
        first = True
        for nr in n_edge:
            if nr == 0 and self.tot_edges == 0 and not first:
                print 'ERROR, no discretization points given'
                sys.exit(0)
            if not (nr == 0):
                first = False
                if self.tot_edges:
                    self.tot_edges += nr - 1
                else:
                    self.tot_edges = nr
        self.grid_edge = sp.empty(self.tot_edges , float)
        left = 0.
        totnr = 0
        first = True
        for nr, right in zip(n_edge, self.surf):
            if nr:
                if first:
                    self.grid_edge[totnr:totnr+nr] = sp.linspace(left, right, nr)
                    totnr += nr
                else:
                    self.grid_edge[totnr:totnr+nr-1] = sp.linspace(left, right, nr)[1:]
                    totnr += nr-1
                first = False
            left = right
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]
        #create a fipy mesh for visualization and fipy computation
        self.mesh_fiber = CylindricalGrid1D(dx=tuple(self.delta_r))
        self.mesh_fiber.periodicBC = False
        self.mesh_fiber = self.mesh_fiber + \
                        ((self.surf_begin[self.ind_first_zone],),)

    def initial_fiber(self):
        """ initial concentration over the domain"""
        if self.fiber_diff:
            #porosity value of cotton fiber
            self.porosity_in = self.cfg.get('fiber.porosity_in')
            #percentage of active component
            self.percentage_active = self.cfg.get('fiber.percentage_active')
            
        self.initial_c1 = sp.empty(self.tot_edges-1, float)
        self.diffusion_coeff = sp.empty(self.tot_edges-1, float)
        self.diffusion_exp_fact = sp.empty(self.tot_edges-1, float)
        st = 0
        surf = self.surf[st]
        for i, pos in enumerate(self.grid):
            while pos > surf:
                st += 1
                surf = self.surf[st]
            if st == 0:
                self.initial_c1[i] = self.init_conc[st](pos) * \
                                            (self.porosity_in / 100.) * \
                                            (self.percentage_active / 100.)
            else:
                self.initial_c1[i] = self.init_conc[st](pos)
            self.diffusion_coeff[i] = self.diff_coef[st]
            self.diffusion_exp_fact[i] = self.diff_exp_fact[st]

        print 'initial mass = ', self.calc_mass(self.initial_c1)

    def calc_mass(self, conc_r):
        """calculate the mass of component present given value in cell center
        This is given by 2 \pi int_r1^r2 C(r)r dr
        
        conc_r: concentration in self.grid
        """
        return sp.sum(conc_r *  self.grid* self.delta_r) * 2. * sp.pi 
        
    def _set_bound_flux(self, flux_edge, w_rep):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, w_rep contains solution in the cell centers
        """
        if self.bound_left == FLUX:
            flux_edge[0] = -self.boundary_fib_left *  self.grid_edge[0]
        else:
            print 'ERROR: boundary type left not implemented'
            sys.exit(0)
        if self.bound_right == FLUX:
            flux_edge[-1] = self.boundary_fib_right *  self.grid_edge[-1]
        else:
            # a transfer coeff to the right
            flux_edge[-1] = -self.boundary_transf_right * w_rep[-1] + self.diffusion_coeff[-1] * \
                            sp.exp(-self.diffusion_exp_fact[-1] * w_rep[-1]/self.grid[-1])*w_rep[-1]/self.grid[-1]

    def _set_bound_fluxu(self, flux_edge, conc_r):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, conc_r contains solution in the cell centers
        """
        if self.bound_left == FLUX:
            flux_edge[0] = -self.boundary_fib_left
        else:
            print 'ERROR: boundary type left not implemented'
            sys.exit(0)
        if self.bound_right == FLUX:
            flux_edge[-1] = self.boundary_fib_right
        else:
            # a transfer coeff to the right
            flux_edge[-1] = -self.boundary_transf_right * conc_r[-1] 
            #print 'flux', flux_edge[-1]

    def f_conc1_ode(self, t, w_rep):
        return self.f_conc1(w_rep, t)

    def f_conc1(self, w_rep, t):
        grid = self.grid
        n_cell = len(grid)
        #Initialize the left side of ODE equations
        diff_w_t = sp.empty(n_cell, float)
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = sp.empty(n_cell+1, float)
        self._set_bound_flux(flux_edge, w_rep)
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = (self.diffusion_coeff[:-1] * \
                            sp.exp(-self.diffusion_exp_fact[:-1] * w_rep[:-1]/self.grid[:-1]) \
                         + self.diffusion_coeff[1:] * \
                            sp.exp(-self.diffusion_exp_fact[1:] * w_rep[1:]/self.grid[1:]))/2.\
                    * self.grid_edge[1:-1] \
                    * (w_rep[1:]/self.grid[1:] - w_rep[:-1]/self.grid[:-1])\
                    / ((self.delta_r[:-1] + self.delta_r[1:])/2.)
        diff_w_t[:]=(flux_edge[1:]-flux_edge[:-1])/self.delta_r[:]
        return diff_w_t
    
    def f_conc1_odeu(self, t, conc_r):
        return self.f_conc1u(conc_r, t)
    
    def f_conc1u(self, conc_r, t):
        grid = self.grid
        n_cell = len(grid)
        #Initialize the left side of ODE equations
        diff_u_t = sp.empty(n_cell, float)
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = sp.empty(n_cell+1, float)
        self._set_bound_fluxu(flux_edge, conc_r)
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = (self.diffusion_coeff[:-1] * \
                            sp.exp(-self.diffusion_exp_fact[:-1] * conc_r[:-1]) \
                         + self.diffusion_coeff[1:] * \
                            sp.exp(-self.diffusion_exp_fact[1:] * conc_r[1:]))/2.\
                    * self.grid_edge[1:-1] \
                    * (conc_r[1:] - conc_r[:-1])\
                    / ((self.delta_r[:-1] + self.delta_r[1:])/2.)
        #diff_u_t[:]=2*(flux_edge[1:]-flux_edge[:-1])/(2*self.grid_edge[1:]*self.delta_r[:]+self.delta_r[:]**2)
        diff_u_t[:]=(flux_edge[1:]-flux_edge[:-1])/(self.grid_edge[1:]**22-self.grid_edge[:-1]**2/2)
        return diff_u_t
    
    def solve_odeint(self):
        initial_w = self.initial_c1 * self.grid
        self.solv=odeint(self.f_conc1, initial_w, self.times)
        self.conc1=self.solv/ self.grid
        self.view_sol(self.times, self.conc1)
        
    def solve_odeintu(self):
        self.solv=odeint(self.f_conc1u, initial_c1, self.times)
        self.conc1=self.solv
        self.view_sol(self.times, self.conc1)
    
    def solve_ode(self):
        self.delta_t = self.times[1]-self.times[0]
        self.initial_t = self.times[0]
        endT = self.times[-1]
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        r = ode(self.f_conc1_ode).set_integrator('vode', method = 'bdf')
        initial_w1 = self.initial_c1 * self.grid
        r.set_initial_value(initial_w1, self.initial_t)#.set_f_params(2.0)
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        while r.successful() and r.t < endT - self.delta_t /10.:
            r.integrate(r.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] = r.y / self.grid
            #print 'mass = ', self.calc_mass(self.conc1[tstep])
            print self.conc1[tstep][-1]
        self.view_sol(self.times, self.conc1)
        
    def solve_odeu(self):
        self.delta_t = self.times[1]-self.times[0]
        self.initial_t = self.times[0]
        endT = self.times[-1]
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        r = ode(self.f_conc1_odeu).set_integrator('vode', method = 'bdf')
        r.set_initial_value(self.initial_c1, self.initial_t)#.set_f_params(2.0)
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        while r.successful() and r.t < endT - self.delta_t /10.:
            r.integrate(r.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] = r.y 
            #print 'mass = ', self.calc_mass(self.conc1[tstep])
            print self.conc1[tstep][-1]
        self.view_sol(self.times, self.conc1)
       

    def solve_fipy(self):
        #using fipy to solve 1D problem in fiber
        self.solution_fiber = CellVariable(name = "fiber concentration", 
                                mesh = self.mesh_fiber,
                                value = self.initial_c1, hasOld = 1)
        self.viewer =  Viewer(vars = self.solution_fiber, datamin=0., datamax= 1.1)
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)

        if self.bound_left == FLUX and self.bound_right == FLUX:
            self.BCs_fiber = (FixedFlux(faces = self.mesh_fiber.getFacesRight(), 
                                    value = self.boundary_fib_right),
                              FixedFlux(faces = self.mesh_fiber.getFacesLeft(), 
                                    value = -self.boundary_fib_left))
        elif self.bound_left == FLUX and self.bound_right == TRANSFER:
            self.BCs_fiber = (FixedFlux(faces = self.mesh_fiber.getFacesRight(), 
                                       value = self.boundary_transf_right \
                                            * self.solution_fiber.getFaceValue()),
                              FixedFlux(faces = self.mesh_fiber.getFacesLeft(), 
                                    value = -self.boundary_fib_left))
        else:
            print 'ERROR: boundary type left not implemented'
            sys.exit(0)
        ## TODO: in following diffusion is given as cellvariable, would a 
        ##      FaceVariable already not be better?
        self.eqX_fiber = TransientTerm() == DiffusionTerm(coeff = 
                            self.diffusion_coeff * 
                            sp.exp(-self.diffusion_exp_fact * self.solution_fiber))
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1[:]
        for time in self.times[1:]:
            self.solve_fipy_step()
            if self.viewer is not None:
                self.viewer.plot()
                #raw_input("please<return>.....")
            tstep += 1
            self.conc1[tstep][:] = self.solution_fiber.getValue()
            #if time == 200.0:
            #    dump.write({'space_position': self.grid, 'conc1': self.conc1[tstep][:]},
            #            filename = utils.OUTPUTDIR + os.sep + 'fipy_t1.gz', extension = '.gz')
            #    print 'finish file'
            #print 'mass = ', self.calc_mass(self.conc1[tstep])
            print self.conc1[tstep][-1]

    def solve_fipy_step(self):
        res = 1e+1
        while res > 1e-8:
            res = self.eqX_fiber.sweep(var = self.solution_fiber,
                                        dt = self.delta_t,
                                        boundaryConditions = self.BCs_fiber)
        self.solution_fiber.updateOld()

    def solve(self):
        """
        Solve the diffusion process in the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        if self.submethod == 'fipy':
            self.solve_fipy()
        else:            
            if self.submethod == 'odeintw':
                self.solve_odeint()
            elif  self.submethod == 'odew':
                self.solve_ode()
            elif self.submethod == 'odeintu':
                self.solve_odeintu()
            elif self.submethod == 'odeu':
                self.solve_odeu()
        self.fiber_surface = sp.empty(len(self.times), float)
        for i in sp.arange(0,len(self.times),1):
            self.fiber_surface[i] = self.conc1[i][-1]

    def view_sol(self, times, conc):
        """
        Show the solution in conc with times.
        conc[i][:] contains solution at time times[i]
        """
        self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber,
                            value = conc[0])
        self.viewer =  Viewer(vars = self.solution_view, datamin=0., datamax=conc[0].max())
        self.viewer.plot()
        for time, con in zip(times[1:], conc[1:]):
            self.solution_view.setValue(con)
            self.viewer.plot()
            #if time == 200.0:
            #    dump.write({'space_position': self.grid, 'conc1': con},
            #                filename = utils.OUTPUTDIR + os.sep + 'ode_t1.gz', extension = '.gz')

    def dump_solution(self): 
        """write out the solution to disk for future use"""
        dump.write({'space_position': self.grid, 'conc': self.conc1},
            filename = utils.OUTPUTDIR + os.sep + 'sol_%s.gz' % self.submethod,
            extension = '.gz')

    def run(self, wait=False, output=False):
        self.create_mesh()
        self.initial_fiber()
        self.solve()

        print 'end mass = ', self.calc_mass(self.conc1[-1])
        if output:
            self.dump_solution()
        if wait:
            raw_input("Finished fiber1d run")
        