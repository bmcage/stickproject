#
# Copyright (C) 2012  B. Malengier
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
import shutil
import sys
import numpy as np
import scipy as sp
from scipy.integrate import ode, trapz
import matplotlib.pyplot as plt
import time

HAVE_ODES = False
try:
    from scikits.odes import ode as sc_ode
    HAVE_ODES = True
except:
    print 'Could not load scikits.odes, odes solver not available'

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import const
import lib.utils.utils as utils

DIR = const.DATA_DIR + os.sep + 'sandbox_domainsplit'
#create outputdir if not existing, delete if present
if not os.path.isdir(DIR):
    os.mkdir(DIR)
else:
    shutil.rmtree(DIR)
    os.mkdir(DIR)
utils.set_outputdir(DIR)

#-------------------------------------------------------------------------
#
#Fipy Imports
#
#-------------------------------------------------------------------------
from fipy import *

#-------------------------------------------------------------------------
#
# DomainslitTestModel class 
#
#-------------------------------------------------------------------------
class DomainsplitTestModel(object):
    """
    Test if we can split a domain in two pieces
    """
    def __init__(self):
        """ 
        Domain is x in (0, L). Source S(x) in (0,l), l<L
        Diffusion D = 1, L =1
        Split in (0,s), s = l + epsilon < L
        and (l, L), with adapted BC and source terms
        """
        self.solver = None
        self.verbose = True
        self.time_period = 30.
        self.delta_t = 0.1
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        #set correct delta_t
        self.delta_t = self.times[1]-self.times[0]
        if self.verbose:
            print "Timestep used :", self.delta_t
        #storage for output
        self.fiber_surface = sp.empty(len(self.times), float)
        #the radial flux at surface. Do times 2 pi radius() to obtain all flux
        self.flux_at_surface = sp.empty(len(self.times), float)
        
        self.L = 1
        self.l = 0.4
        self.s = 0.5
        self.diffusion_coeff = 1e-2
        self.tot_edges = 41
        self.tot_edges_sp1 = 41
        self.tot_edges_sp2 = 41
            
        #data for stepwise operation
        self.initialized = False
        self.__userdata = None

        self.plotevery = 5

    def set_userdata(self, data):
        """
        Read the data of the concentration in the void space and overwrite the 
        default value 
        """
        self.__userdata = data

    def get_userdata(self):
        return self.__userdata

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh
        grid: the space position of each central point for every cell element;
        """
        self.grid_edge = sp.linspace(0., self.L, self.tot_edges)  
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]
        
        ##grid for split1 and split2
        self.grid_edge_sp1 = sp.linspace(0., self.s, self.tot_edges_sp1)  
        #construct cell centers from this
        self.grid_sp1 = (self.grid_edge_sp1[:-1] + self.grid_edge_sp1[1:])/2.
        #obtain cell sizes
        self.delta_r_sp1 = self.grid_edge_sp1[1:] - self.grid_edge_sp1[:-1]
        self.grid_edge_sp2 = sp.linspace(self.l, self.L, self.tot_edges_sp2)  
        #construct cell centers from this
        self.grid_sp2 = (self.grid_edge_sp2[:-1] + self.grid_edge_sp2[1:])/2.
        #obtain cell sizes
        self.delta_r_sp2 = self.grid_edge_sp2[1:] - self.grid_edge_sp2[:-1]

        #create a fipy mesh for visualization and fipy computation
        self.mesh_fiber = CylindricalGrid1D(dx=tuple(self.delta_r))
        self.mesh_fiber.periodicBC = False

    def initial_fiber(self):
        """ initial concentration over the domain"""
        self.initial_c1 = sp.zeros(self.tot_edges-1, float)
        for i, grid in enumerate(self.grid):
            if grid < self.l:
                self.initial_c1[i] = 1.
        self.initial_w1 = self.initial_c1 * self.grid
        self.volume = self.calc_volume()
        if self.verbose:
            print 'initial mass = ', self.calc_mass(self.initial_c1)
        ## the init for the split domains
        self.initial_c1_sp2 = sp.zeros(self.tot_edges_sp2-1, float)
        self.initial_c1_sp1 = sp.zeros(self.tot_edges_sp1-1, float)
        for i, grid in enumerate(self.grid_sp1):
            if grid < self.l:
                self.initial_c1_sp1[i] = 1.
        self.initial_w1_sp1 = self.initial_c1_sp1 * self.grid_sp1
        self.initial_w1_sp2 = self.initial_c1_sp2 * self.grid_sp2

    def calc_mass(self, conc_r):
        """calculate the mass of component present given value in cell center
        This is given by 2 \pi int_r1^r2 C(r)r dr
        
        conc_r: concentration in self.grid
        """
        return sp.sum(conc_r * (sp.power(self.grid_edge[1:], 2) - 
                                sp.power(self.grid_edge[:-1], 2)) 
                        ) * sp.pi 

    def calc_volume(self):
        """calculate the volume over which the compound can move. We have
           Cavg = mass/volume
        """
        return sp.sum((sp.power(self.grid_edge[1:],2) - 
                            sp.power(self.grid_edge[:-1],2)) 
                        ) * sp.pi 

    def _set_bound_flux(self, flux_edge, w_rep, t):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, w_rep contains solution in the cell centers
        This is the flux of w over radius! So times 2 pi R for full flux w, or
        for C times 2 pi.
        """
        # no flux
        flux_edge[0] = 0.
        #and correct, flux needed is r times the x coord flux
        flux_edge[0] = flux_edge[0] * self.grid_edge[0]
        flux_edge[-1] = 0.
        #and correct, flux needed is r times the x coord flux
        flux_edge[-1] = flux_edge[-1] * self.grid_edge[-1]
    
    def _set_bound_flux_sp1(self, flux_edge, w_rep, t):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, w_rep contains solution in the cell centers
        This is the flux of w over radius! So times 2 pi R for full flux w, or
        for C times 2 pi.
        """
        # no flux
        flux_edge[0] = 0.
        #and correct, flux needed is r times the x coord flux
        flux_edge[0] = flux_edge[0] * self.grid_edge_sp1[0]
        flux_edge[-1] = 0.
        #and correct, flux needed is r times the x coord flux
        flux_edge[-1] = flux_edge[-1] * self.grid_edge_sp1[-1]

    
    def _set_bound_flux_sp2(self, flux_edge, w_rep, t):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, w_rep contains solution in the cell centers
        This is the flux of w over radius! So times 2 pi R for full flux w, or
        for C times 2 pi.
        """
        # no flux
        flux_edge[0] = 0.
        #and correct, flux needed is r times the x coord flux
        flux_edge[0] = flux_edge[0] * self.grid_edge_sp2[0]
        flux_edge[-1] = 0.
        #and correct, flux needed is r times the x coord flux
        flux_edge[-1] = flux_edge[-1] * self.grid_edge_sp2[-1]

    def f_conc1_odes(self, t, w_rep, diff_w_t):
        grid = self.grid
        n_cell = len(grid)
        #Initialize the left side of ODE equations
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = self.__tmp_flux_edge
        self._set_bound_flux(flux_edge, w_rep, t)
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = - self.diffusion_coeff \
                        * self.grid_edge[1:-1] \
                        * (w_rep[1:]/self.grid[1:] - w_rep[:-1]/self.grid[:-1])\
                        / ((self.delta_r[:-1] + self.delta_r[1:])/2.)
        diff_w_t[:] = (flux_edge[:-1] - flux_edge[1:]) / self.delta_r[:] 
    
    def f_conc1_odes_sp1(self, t, w_rep, diff_w_t):
        grid = self.grid_sp1
        n_cell = len(grid)
        #Initialize the left side of ODE equations
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = self.__tmp_flux_edge_sp1
        self._set_bound_flux_sp1(flux_edge, w_rep, t)
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = - self.diffusion_coeff \
                * self.grid_edge_sp1[1:-1] \
                * (w_rep[1:]/self.grid_sp1[1:] - w_rep[:-1]/self.grid_sp1[:-1])\
                / ((self.delta_r_sp1[:-1] + self.delta_r_sp1[1:])/2.)
        diff_w_t[:] = (flux_edge[:-1] - flux_edge[1:]) / self.delta_r_sp1[:] 

    def f_conc1_odes_sp2(self, t, w_rep, diff_w_t):
        grid = self.grid_sp2
        n_cell = len(grid)
        #Initialize the left side of ODE equations
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = self.__tmp_flux_edge_sp2
        self._set_bound_flux_sp2(flux_edge, w_rep, t)
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = - self.diffusion_coeff \
                * self.grid_edge_sp2[1:-1] \
                * (w_rep[1:]/self.grid_sp2[1:] - w_rep[:-1]/self.grid_sp2[:-1])\
                / ((self.delta_r_sp2[:-1] + self.delta_r_sp2[1:])/2.)
        diff_w_t[:] = (flux_edge[:-1] - flux_edge[1:]) / self.delta_r_sp2[:] 

    def solve_odes_init(self):
        """
        Initialize the cvode solver
        """
        if not HAVE_ODES:
            raise Exception, 'Not possible to solve with given method, scikits.odes not available'
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        self.step_old_sol = self.initial_w1
        self.step_old_sol_sp1 = self.initial_w1_sp1
        self.step_old_sol_sp2 = self.initial_w1_sp2
        #data storage
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        self.ret_y = sp.empty(len(self.initial_c1), float)
        
        self.conc1_sp1 = np.empty((len(self.times), len(self.initial_c1_sp1)), float)
        self.ret_y_sp1 = sp.empty(len(self.initial_c1_sp1), float)
        self.conc1_sp2 = np.empty((len(self.times), len(self.initial_c1_sp2)), float)
        self.ret_y_sp2 = sp.empty(len(self.initial_c1_sp2), float)

        self.conc1[0][:] = self.initial_c1
        self.conc1_sp1[0][:] = self.initial_c1_sp1
        self.conc1_sp2[0][:] = self.initial_c1_sp2
        n_cell = len(self.grid)
        self.__tmp_diff_w_t = sp.empty(n_cell, float)
        self.__tmp_flux_edge = sp.empty(n_cell+1, float)
        n_cell_sp1 = len(self.grid_sp1)
        self.__tmp_diff_w_t_sp1 = sp.empty(n_cell_sp1, float)
        self.__tmp_flux_edge_sp1 = sp.empty(n_cell_sp1+1, float)
        n_cell_sp2 = len(self.grid_sp2)
        self.__tmp_diff_w_t_sp2 = sp.empty(n_cell_sp2, float)
        self.__tmp_flux_edge_sp2 = sp.empty(n_cell_sp2+1, float)
        
        self.tstep = 0
        self.solve_odes_reinit()
        self.initialized = True

    def solve_odes_reinit(self):
        """
        Reinitialize the cvode solver to start again
        """
        self.initial_t = self.times[0]
        self.solver = sc_ode('cvode', self.f_conc1_odes,
                             max_steps=50000, lband=1, uband=1)
        self.solver.init_step(self.step_old_time, self.step_old_sol)
        self.solver_sp1 = sc_ode('cvode', self.f_conc1_odes_sp1,
                             max_steps=50000, lband=1, uband=1)
        self.solver_sp1.init_step(self.step_old_time, self.step_old_sol_sp1)
        self.solver_sp2 = sc_ode('cvode', self.f_conc1_odes_sp2,
                             max_steps=50000, lband=1, uband=1)
        self.solver_sp2.init_step(self.step_old_time, self.step_old_sol_sp2)

    def do_solve_step(self, stoptime):
        """
        Solve up to time t. This does:
           1. solve the global model up to t
           2. solve split part 1 up to t
           3. set correct flux/source term for part 2
           4. solve split part 2 up to t
           3. set correct flux/source term for part 1
        """
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            #solve of global model
            time, conc = self.do_step_odes(t)
            ## TODO Now the split problem, and compare it

        if self.plotevery:
            if self.viewerplotcount == 0:
                self.solution_view.setValue(conc)
                self.viewer.axes.set_title('time %s' %str(time))
                self.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'conc%s.png' % str(int(10*time)))
            self.viewerplotcount += 1
            self.viewerplotcount = self.viewerplotcount % self.plotevery
        return time, conc

    def solve_odes(self, run_per_step = None, viewend = True):
        #self.solve_odes_init()
        endT = self.times[-1]
        self.initial_w1 = self.initial_c1 * self.grid
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        for time in self.times[1:]:
            flag, realtime = self.solver.step(time, self.conc1[tstep+1])
            if flag != 0:
                print 'ERROR: unable to compute solution, flag', flag
                break
            if self.verbose:
                print 'INFO: fibermodel at t = ', realtime
            tstep += 1
            self.conc1[tstep][:] = self.conc1[tstep][:] / self.grid[:]
            self.fiber_surface[tstep] = self.conc1[tstep][-1]
            self.flux_at_surface[tstep] = 0.
            if run_per_step:
                run_per_step(self, time, self.conc1[tstep])

            #print 'mass = ', self.calc_mass(self.conc1[tstep])
        if viewend:
            self.view_sol(self.times, self.conc1)

    def do_step_odes(self, stoptime, needreinit=True):
        """
        Solve the fibermodel up to stoptime, continuing from the present
        state, return time and r * concentration after step
        It is needed that run_init and solve_init method have been called
        before calling this method.
        
        if needreinit = True, the solver is initialized first with the
            data present in step_old_time amd step_old_sol
            
        Return: concentration over the grid
        """
        assert stoptime > self.step_old_time, "%f > %f" % (stoptime, self.step_old_time)
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'
        if needreinit:
            self.solve_odes_reinit()
        else:
            self.solver.set_options(tstop=stoptime)
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            flag, realtime = self.solver.step(t, self.ret_y)
            if flag < 0:
                raise Exception, 'could not find solution'
        
        self.step_old_time = realtime
        self.step_old_sol = self.ret_y
        assert np.allclose(realtime, stoptime, atol=1e-6, rtol=1e-6)
        return realtime, self.ret_y / self.grid

    def solve(self):
        """
        Solve the diffusion process in the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        def run_every_step(object, time, conc):
            if self.plotevery:
                if object.viewerplotcount == 0:
                    object.solution_view.setValue(conc)
                    object.viewer.axes.set_title('time %s' %str(time))
                    object.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'conc%s.png' % str(int(10*time)))
                object.viewerplotcount += 1
                object.viewerplotcount = self.viewerplotcount % self.plotevery

        if self.plotevery:
            self.solution_view = CellVariable(name = "fiber concentration", 
                    mesh = self.mesh_fiber,
                    value = self.conc1[0][:])
            self.viewer =  Matplotlib1DViewer(vars = self.solution_view, 
                                datamin=0., 
                                datamax=1.2 * self.conc1[0].max())
            self.viewer.axes.set_title('time 0.0')
            self.viewer.plot()
            self.viewerplotcount = 1
        self.solve_odes(run_per_step =run_every_step, viewend=False)

        if self.verbose:
            print 'Finished the fiber calculation'

    def do_step(self, stoptime, needreinit=True):
        """
        Solve the diffusion process in the fiber up to stoptime 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. 
        The resulting time and r*concentration is returned
        """
        return self.do_step_odes(stoptime, needreinit)

    def solve_init(self):
        """
        Initialize the solvers so they can be solved stepwize
        """
        self.solve_odes_init()

    def view_sol(self, times, conc):
        """
        Show the solution in conc with times.
        conc[i][:] contains solution at time times[i]
        """
        self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber,
                            value = conc[0][:])
        name = self.solution_view.name                    
        if self.plotevery:
            self.viewer =  Matplotlib1DViewer(vars = self.solution_view, datamin=0., datamax=conc.max()+0.20*conc.max())
        self.viewerplotcount = 0
        for time, con in zip(times[1:], conc[1:][:]):
            self.solution_view.setValue(con)
            if self.plotevery:
                self.viewer.axes.set_title('Fiber Conc vs radius at time %s' %str(time))
                self.viewer.axes.set_xlabel('Radius')
                self.viewer.axes.set_ylabel('Conc')
                if self.viewerplotcount == 0:
                   self.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'fiber%sconc%08.4f.png' % (name,time))
                #else:
                #    self.viewer.plot()
                    
                self.viewerplotcount += 1
                self.viewerplotcount = self.viewerplotcount % self.plotevery   

    def view_time(self, times, conc, title=None):
        draw_time = times/(3600.*24.*30.) # convert seconds to months
        draw_conc = conc *1.0e4
        plt.figure(num=None)
        plt.plot(draw_time, draw_conc, '-', color = 'red')
        plt.xlabel('Time (month)')
        plt.ylabel(title)
        plt.show()

    def dump_solution(self): 
        """write out the solution to disk for future use"""
        if self.method == 'FVM':
            dump.write({'space_position': self.grid, 'conc': self.conc1},
                filename = utils.OUTPUTDIR + os.sep + 'sol.gz',
                extension = '.gz')
        dump.write({'time_step':self.times, 'flux': self.flux_at_surface},
            filename = utils.OUTPUTDIR + os.sep + 'flux_boundary', 
            extension = '.gz')

    def run_init(self):
        self.create_mesh()
        self.initial_fiber()
        
    def run(self, wait=False, output=False):
        self.run_init()
        if not self.initialized:
            self.solve_init()
        
        self.solution_view = CellVariable(name = "fiber concentration", 
                mesh = self.mesh_fiber,
                value = self.conc1[0][:])
        self.viewer =  Matplotlib1DViewer(vars = self.solution_view, 
                            datamin=0., 
                            datamax=1.2 * self.conc1[0].max())
        self.viewer.axes.set_title('time 0.0')
        self.viewer.plot()
        self.viewerplotcount = 1
        
        for ind, t in enumerate(self.times[1:]):
            #print 'solving t', t
            rt, conc = self.do_solve_step(t)
            self.conc1[ind+1][:] = conc[:]
            self.fiber_surface[ind+1] = conc[-1]
            self.flux_at_surface[ind+1] = 0.

        if output:
            self.dump_solution()
        if wait:
            raw_input("Finished domainsplit run")

    def __del__(self):
        #remove the memory
        if self.solver:
            print 'del self.solver'
            del self.solver
        self.solver = None

if __name__ == '__main__':
    dsm = DomainsplitTestModel()
    dsm.run()