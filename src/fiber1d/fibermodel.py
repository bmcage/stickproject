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
from fiber1d.config import (METHOD, FLUX, TRANSFER, BOUND_TYPE, FIBER_FORM,
                    CIRCLE, ELLIPSE)

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
        print 'the time period', self.time_period
        self.delta_t = self.cfg.get('time.dt')
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        #set correct delta_t
        self.delta_t = self.times[1]-self.times[0]
        print "Timestep used in fiber model:", self.delta_t
        #storage for output
        self.fiber_surface = sp.empty(len(self.times), float)
        self.transfer_boundary = sp.empty(len(self.times), float)
        
        #print 'the times', self.times
        #self.delta_t = 0.1#self.times[1] - self.times[0]
        #read the initial and boundary information for fiber
        self.n_edge = self.cfg.get('fiber.n_edge') #discretize the fiber radius
        self.bound_left = BOUND_TYPE[self.cfg.get('boundary.type_left')]
        self.bound_right = BOUND_TYPE[self.cfg.get('boundary.type_right')]
        self.boundary_fib_left = self.cfg.get('boundary.boundary_fib_left')
        self.boundary_fib_right = self.cfg.get('boundary.boundary_fib_right')
        self.boundary_transf_right = self.cfg.get('boundary.transfer_right')
        
        #data for stepwise operation
        self.initialized = False
        
        self.__Rf_pure = None
        self.__Rf = None

        self.plotevery = self.cfg.get("plot.plotevery")
        
        self.verbose = self.cfg.get('general.verbose')
        
    def radius_pure(self):
        """method that returns the radius of the fiber seen as a circle,
           without coatings
        """
        if self.__Rf_pure:
            return self.__Rf_pure
        rad = self.cfg.get('fiber.radius_pure_fiber')
        form = FIBER_FORM[self.cfg.get('fiber.form')]
        if form == CIRCLE:
            pass
        elif form == ELLIPSE:
            #radius_pure_fiber is the long axis
            ecc = self.cfg.get('fiber.eccentricity')
            if not (ecc>0. and ecc <= 1.):
                raise Exception, 'Eccentricity must be between 0 and 1'
            #we determine radius of the fiber with same surface area
            #pi a b = pi r**2
            a = rad
            b = sqrt(a**2 * (1-ecc**2))
            rad = sqrt(a*b)
        else:
            raise Exception, 'Fiber form is not supported for a 1D fiber model'
        self.__Rf_pure = rad
        return self.__Rf_pure
    
    def radius(self):
        """ method that returns the total radius of the fiber
        """
        if self.__Rf:
            return self.__Rf
        rad = self.radius_pure()
        for i in range(self.cfg.get('fiber.nrlayers')):
            section = 'fiberlayer_%i' % i
            rad += self.cfg.get(section + '.thickness')
        self.__Rf = rad
        return self.__Rf

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
            self.porosity = [0.]
        else:
            n_edge = [self.cfg.get('fiber.n_edge')]
            self.diff_coef = [self.cfg.get('fiber.diffusion_coef')]
            self.diff_exp_fact = [self.cfg.get('fiber.diffusion_polymer_exp_factor')]
            self.init_conc = [eval(self.cfg.get('fiber.init_conc')) ]
            self.ind_first_zone = 0
            self.porosity = [self.cfg.get('fiber.porosity_in')]
        self.nrlayers = self.cfg.get('fiber.nrlayers')
        self.surf_begin = [0.]
        self.surf = [self.radius_pure()]
        for i in range(self.nrlayers):
            section = 'fiberlayer_%i' % i
            n_edge += [self.cfg.get(section + '.n_edge')]
            self.surf_begin += [self.surf[-1]]
            self.surf += [self.surf[-1] + self.cfg.get(section + '.thickness')]
            self.diff_coef += [self.cfg.get(section + '.diffusion_coef')]
            self.diff_exp_fact += [self.cfg.get(section + '.diffusion_polymer_exp_factor')]
            self.init_conc += [eval(self.cfg.get(section + '.init_conc')) ]
            self.porosity += [self.cfg.get(section + '.porosity_layer')]

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
                    #left = 0.
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
        print 'mesh', self.grid, self.grid_edge, self.delta_r
        #create a fipy mesh for visualization and fipy computation
        self.mesh_fiber = CylindricalGrid1D(dx=tuple(self.delta_r))
        self.mesh_fiber.periodicBC = False
        self.mesh_fiber = self.mesh_fiber + \
                        ((self.surf_begin[self.ind_first_zone],),)

    def initial_fiber(self):
        """ initial concentration over the domain"""
        if self.fiber_diff:
            #percentage of active component
            self.percentage_active = self.cfg.get('fiber.percentage_active')
        print 'TOT EDGES', self.tot_edges
        self.initial_c1 = sp.empty(self.tot_edges-1, float)
        self.diffusion_coeff = sp.empty(self.tot_edges-1, float)
        self.diffusion_exp_fact = sp.empty(self.tot_edges-1, float)
        self.porosity_domain = sp.empty(self.tot_edges - 1, float)
        st = 0
        surf = self.surf[st]
        for i, pos in enumerate(self.grid):
            while pos > surf:
                st += 1
                surf = self.surf[st]
            if st == 0:
                self.initial_c1[i] = self.init_conc[st](pos) * \
                                            (self.porosity[st]) * \
                                            (self.percentage_active / 100.)
            else:
                self.initial_c1[i] = self.init_conc[st](pos) * self.porosity[st]
            self.diffusion_coeff[i] = self.diff_coef[st]
            self.diffusion_exp_fact[i] = self.diff_exp_fact[st]
            self.porosity_domain[i] = self.porosity[st]
        self.initial_w1 = self.initial_c1 * self.grid
        print 'initial mass = ', self.calc_mass(self.initial_c1)

    def calc_mass(self, conc_r):
        """calculate the mass of component present given value in cell center
        This is given by 2 \pi int_r1^r2 C(r)r dr
        
        conc_r: concentration in self.grid
        """
        return sp.sum(conc_r * self.grid* self.delta_r * self.porosity_domain) * 2. * sp.pi 
        
    def _set_bound_flux(self, flux_edge, w_rep):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, w_rep contains solution in the cell centers
        """
        if self.bound_left == FLUX:
            flux_edge[0] = -self.boundary_fib_left * self.grid_edge[0]
        else:
            print 'ERROR: boundary type left not implemented'
            sys.exit(0)
        if self.bound_right == FLUX:
            flux_edge[-1] = self.boundary_fib_right * self.grid_edge[-1]
        else:
            flux_edge[-1] = self.boundary_transf_right * w_rep[-1]  \
                            *self.porosity_domain[-1]

    def _set_bound_fluxu(self, flux_edge, conc_r):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, conc_r contains solution in the cell centers
        """
        if self.bound_left == FLUX:
            flux_edge[0] = -self.boundary_fib_left * self.grid_edge[0]
        else:
            print 'ERROR: boundary type left not implemented'
            sys.exit(0)
        if self.bound_right == FLUX:
            flux_edge[-1] = self.boundary_fib_right * self.grid_edge[-1]
        else:
            # a transfer coeff to the right
            flux_edge[-1] = self.boundary_transf_right * conc_r[-1] \
                             * self.grid_edge[-1]  \
                             * self.porosity_domain[-1]
            #print 'flux', flux_edge[-1]

    def f_conc1_ode(self, t, w_rep):
        return self.f_conc1(w_rep, t)

    def f_conc1(self, w_rep, t):
        #print 'w_rep', w_rep
        grid = self.grid
        n_cell = len(grid)
        #Initialize the left side of ODE equations
        diff_w_t = sp.empty(n_cell, float)
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = sp.empty(n_cell+1, float)
        self._set_bound_flux(flux_edge, w_rep)
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
##        print (len(self.porosity_domain[:-1]), 
##            len(self.diffusion_coeff[:-1]),
##            len(self.diffusion_exp_fact[:-1]), 
##            len(w_rep[:-1]), 
##            len(self.grid[:-1]), 
##            len(self.grid_edge[1:-1]), 
##            len(self.grid[:-1]), 
##            )
        flux_edge[1:-1] = (self.porosity_domain[:-1] * self.diffusion_coeff[:-1] * 
                            sp.exp(-self.diffusion_exp_fact[:-1] * w_rep[:-1]/self.grid[:-1]) \
                         + self.porosity_domain[1:] * self.diffusion_coeff[1:] * 
                            sp.exp(-self.diffusion_exp_fact[1:] * w_rep[1:]/self.grid[1:]))/2.\
                        * self.grid_edge[1:-1] \
                        * (w_rep[1:]/self.grid[1:] - w_rep[:-1]/self.grid[:-1])\
                        / ((self.delta_r[:-1] + self.delta_r[1:])/2.)
        diff_w_t[:]=(flux_edge[1:]-flux_edge[:-1])/self.delta_r[:] / self.porosity_domain[:]
        print flux_edge[-2],flux_edge[-1], diff_w_t[-1]
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
        flux_edge[1:-1] = (self.porosity_domain[:-1] * self.diffusion_coeff[:-1] * 
                            sp.exp(-self.diffusion_exp_fact[:-1] * conc_r[:-1]) \
                         + self.porosity_domain[1:] * self.diffusion_coeff[1:] * 
                            sp.exp(-self.diffusion_exp_fact[1:] * conc_r[1:]))/2.\
                        * self.grid_edge[1:-1] \
                        * (conc_r[1:] - conc_r[:-1])\
                        / ((self.delta_r[:-1] + self.delta_r[1:])/2.)
        diff_u_t[:] = 2.*((flux_edge[1:]-flux_edge[:-1])
                        /(self.grid_edge[1:]**2-self.grid_edge[:-1]**2)
                        / self.porosity_domain[:])
        return diff_u_t
    
    def solve_odeint(self):
        initial_w = self.initial_c1 * self.grid
        self.solv=odeint(self.f_conc1, initial_w, self.times, rtol=1.0e-9, atol = 1.0e-10)
        self.conc1=self.solv/ self.grid
        self.view_sol(self.times, self.conc1)
        
    def solve_odeintu(self):
        self.solv=odeint(self.f_conc1u, initial_c1, self.times)
        self.conc1=self.solv
        self.view_sol(self.times, self.conc1)

    def solve_ode_init(self):
        """
        Initialize the ode solver
        """
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        self.step_old_sol = self.initial_w1
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        self.conc1[0][:] = self.initial_c1
        self.tstep = 0
        self.solve_ode_reinit()
        self.initialized = True

    def solve_ode_reinit(self):
        """
        Reinitialize the ode solver to start again
        """
        self.initial_t = self.times[0]
        self.solver = ode(self.f_conc1_ode).set_integrator('vode', 
                            method = 'bdf',
                            nsteps=5000)
        self.solver.set_initial_value(self.step_old_sol, self.step_old_time)

    def solve_ode_step(self, step, needreinit=True):
        """Solve the fibermodel for one step, continuing from the present
           state, return the concentration after step
        """
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'
        if needreinit:
            self.solve_ode_reinit()
        curt = self.solver.t
        #print 'curt value', curt
        while self.solver.successful() and self.solver.t < curt + step - self.delta_t /10. and \
        self.solver.t < self.step_old_time + step - self.delta_t /10.:
            #print 'length of solution', len(self.solver.y)
            #print 'length of grid', len(self.grid)
            self.solver.integrate(self.solver.t + self.delta_t)
            self.tstep += 1
            self.conc1[self.tstep][:] = self.solver.y / self.grid
            #print 'conc1', self.conc1[self.tstep][:]
            self.fiber_surface[self.tstep] = self.conc1[self.tstep][-1]
            self.transfer_boundary[self.tstep] = self.boundary_transf_right * self.fiber_surface[self.tstep]
            #print 'mass = ', self.calc_mass(self.conc1[tstep])
        #return the concentration after step
        #self.initial_w1 = self.initial_c1 * self.grid
        self.solver.integrate(curt + step)
        self.step_old_time += step
        self.step_old_sol = self.solver.y
        assert self.solver.t == self.step_old_time, "%f %f" % (self.solver.t, self.step_old_time)
        return self.solver.y
        
    def solve_ode(self):
        self.solve_ode_init()
        endT = self.times[-1]
        self.initial_w1 = self.initial_c1 * self.grid
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        while self.solver.successful() and self.solver.t < endT - self.delta_t /10.:
            self.solver.integrate(self.solver.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] = self.solver.y / self.grid
            self.fiber_surface[tstep] = self.conc1[tstep][-1]
            self.transfer_boundary[tstep] = self.boundary_transf_right * self.fiber_surface[tstep]
            #print 'mass = ', self.calc_mass(self.conc1[tstep])

        self.view_sol(self.times, self.conc1)
        
    def solve_odeu(self):
        self.initial_t = self.times[0]
        endT = self.times[-1]
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        r = ode(self.f_conc1_odeu).set_integrator('vode', method = 'bdf',nsteps = 10000)
        r.set_initial_value(self.initial_c1, self.initial_t)#.set_f_params(2.0)
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        while r.successful() and r.t < endT - self.delta_t /10.:
            r.integrate(r.t + self.delta_t)
            tstep += 1
            self.conc1[tstep][:] = r.y 
            self.fiber_surface[tstep] = self.conc1[tstep][-1]
            self.transfer_boundary[tstep] = self.boundary_transf_right * self.fiber_surface[tstep]
            #print 'mass = ', self.calc_mass(self.conc1[tstep])
        self.view_sol(self.times, self.conc1)

    def solve_fipy(self):
        #using fipy to solve 1D problem in fiber
        self.solution_fiber = CellVariable(name = "fiber concentration", 
                                mesh = self.mesh_fiber,
                                value = self.initial_c1 * self.porosity_domain, hasOld = 1)
        self.viewer =  Viewer(vars = self.solution_fiber / self.porosity_domain, datamin=0., datamax= 1.1)
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
                            sp.exp(-self.diffusion_exp_fact * self.solution_fiber / self.porosity_domain)*
                            self.porosity_domain) 
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1[:]
        print 'solving for times', self.times[1:]
        for time in self.times[1:]:
            self.solve_fipy_sweep()
##            if self.viewer is not None:
##                self.viewer.plot()
                #raw_input("please<return>.....")
            tstep += 1
            self.conc1[tstep][:] = self.solution_fiber.getValue()
            self.fiber_surface[tstep] = self.conc1[tstep][-1]
            self.transfer_boundary[tstep] = self.boundary_transf_right * self.fiber_surface[tstep]
            #if time == 200.0:
            #    dump.write({'space_position': self.grid, 'conc1': self.conc1[tstep][:]},
            #            filename = utils.OUTPUTDIR + os.sep + 'fipy_t1.gz', extension = '.gz')
            #    print 'finish file'
            #print 'mass = ', self.calc_mass(self.conc1[tstep])

    def solve_fipy_sweep(self):
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
            elif  self.submethod == 'odew_step':
                self.solve_ode_init()
                self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber,
                            value = self.conc1[0][:])
                self.viewer =  Matplotlib1DViewer(vars = self.solution_view, 
                                        datamin=0., 
                                        datamax=1.2 * self.conc1[0].max())
                self.viewer.axes.set_title('time 0.0')
                self.viewer.plot()
                self.viewerplotcount = 1
                for nrt, time in enumerate(self.times[1:]):
                    res = self.solve_ode_step(self.delta_t, needreinit=False)
                    if self.viewerplotcount == 0:
                        self.solution_view.setValue(self.conc1[nrt+1])
                        self.viewer.axes.set_title('time %s' %str(time))
                        self.viewer.plot()
                    self.viewerplotcount += 1
                    self.viewerplotcount = self.viewerplotcount % self.plotevery
            elif self.submethod == 'odeintu':
                self.solve_odeintu()
            elif self.submethod == 'odeu':
                self.solve_odeu()
        print 'finished the fiber calculation'
        self.view_time(self.times, self.transfer_boundary)

    def solve_step(self, step):
        """
        Solve the diffusion process in the fiber over a step. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. 
        The resulting concentration is returned
        The mesh in this 1-D condition is 
        uniform
        """
        if self.submethod == 'fipy':
            res = self.solve_fipy_step()
        else:            
            if self.submethod == 'odeintw':
                raise Exception, 'Not supported'
            elif  self.submethod == 'odew':
                res = self.solve_ode_step(step)
            elif self.submethod == 'odeintu':
                raise Exception, 'Not supported'
            elif self.submethod == 'odeu':
                raise Exception, 'Not supported'
        return res

    def solve_init(self):
        """
        Initialize the solvers so they can be solved stepwize
        """
        if self.submethod == 'fipy':
            self.solve_fipy()
        else:            
            if self.submethod == 'odeintw':
                raise Exception, 'Not supported to step for odeint'
                self.solve()
            elif  self.submethod == 'odew':
                self.solve_ode_init()
            elif self.submethod == 'odeintu':
                raise Exception, 'Not supported to step for odeint'
                self.solve()
            elif self.submethod == 'odeu':
                self.solve_ode_init()

    def view_sol(self, times, conc):
        """
        Show the solution in conc with times.
        conc[i][:] contains solution at time times[i]
        """
        self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber,
                            value = conc[0][:])
        self.viewer =  Matplotlib1DViewer(vars = self.solution_view, datamin=0., datamax=conc.max()+0.20*conc.max())
        self.viewerplotcount = 1
        for time, con in zip(times[1:], conc[1:][:]):
            if self.viewerplotcount == 0:
                self.solution_view.setValue(con)
                self.viewer.plot()
                self.viewer.axes.set_title('time %s' %str(time))
            self.viewerplotcount += 1
            self.viewerplotcount = self.viewerplotcount % self.plotevery
    
    def view_time(self, times, conc):
        draw_time = times/(3600.*24.*30.) # convert seconds to months
        draw_conc = conc *1.0e4
        plt.figure()
        plt.plot(draw_time, draw_conc, '-', color = 'red')
        #plt.xlim(0.0, 3.0)
        plt.xlabel('Time (month)')
        plt.ylabel('Flux of DEET ($\mathrm{mg\cdot cm/s}$)')
        plt.draw()

    def dump_solution(self): 
        """write out the solution to disk for future use"""
        dump.write({'space_position': self.grid, 'conc': self.conc1},
            filename = utils.OUTPUTDIR + os.sep + 'sol_%s.gz' % self.submethod,
            extension = '.gz')
        dump.write({'time_step':self.times, 'flux':self.transfer_boundary},
            filename = utils.OUTPUTDIR + os.sep + 'flux_boundary', 
            extension = '.gz')

    def run_init(self):
        self.create_mesh()
        self.initial_fiber()
        
    def run(self, wait=False, output=False):
        self.run_init()
        if not self.initialized:
            self.solve_init()
        self.solve()

        print 'end mass = ', self.calc_mass(self.conc1[-1])
        if output:
            self.dump_solution()
        if wait:
            raw_input("Finished fiber1d run")

    def run_step(self, step):
        if not self.initialized:
            self.solve_init()
        return self.solve_step(step)
