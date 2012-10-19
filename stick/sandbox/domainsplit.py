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
from math import pi

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
        Domain is x in (0, L). Source S(x) in (l,s), l<L
        Diffusion D = 1, L =1
        Split in (0,s), s = l + epsilon < L
        and (l, L), with adapted BC and source terms
        """
        self.solver = None
        self.verbose = True
        self.plotting = False
        self.time_period = 20
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
        self.volume_overlap = pi * (self.s**2 - self.l**2)
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
        self.grid_edge_sp = []
        self.grid_edge_sp1 = sp.linspace(0., self.s, self.tot_edges_sp1)
        #we add self.l here
        for ind, grid in enumerate(self.grid_edge_sp1):
            if grid > self.l: 
                break
            self.grid_edge_sp.append(grid)
        left = self.grid_edge_sp1[ind-1]
        right = self.grid_edge_sp1[ind]
        if self.l - left < right - self.l:
            self.grid_edge_sp1[ind-1] = self.l
            self.grid_edge_sp[-1] = self.l
        else:
            assert ind != len(self.grid_edge_sp1) - 1
            self.grid_edge_sp1[ind] = self.l
            self.grid_edge_sp.append(self.l)
        #construct cell centers from this
        self.grid_sp1 = (self.grid_edge_sp1[:-1] + self.grid_edge_sp1[1:])/2.
        #obtain cell sizes
        self.delta_r_sp1 = self.grid_edge_sp1[1:] - self.grid_edge_sp1[:-1]
        self.grid_edge_sp2 = sp.linspace(self.l, self.L, self.tot_edges_sp2)
        #we add self.s here
        for ind, grid in enumerate(self.grid_edge_sp2):
            if grid > self.s: 
                break
        if ind == 1:
            tmp = sp.linspace(self.l, self.L, self.tot_edges_sp2-1)
            self.grid_edge_sp2[1] = self.s
            self.grid_edge_sp2[2:] = tmp[1:]
        else:
            left = self.grid_edge_sp2[ind-1]
            right = self.grid_edge_sp2[ind]
            if self.s - left < right - self.s:
                self.grid_edge_sp2[ind-1] = self.s
            else:
                self.grid_edge_sp2[ind] = self.s
        for grid in self.grid_edge_sp2[1:]:
            self.grid_edge_sp.append(grid)
        self.grid_edge_sp = np.asarray(self.grid_edge_sp)
        self.tot_edges_sp = len(self.grid_edge_sp)
        #construct cell centers from this
        self.grid_sp2 = (self.grid_edge_sp2[:-1] + self.grid_edge_sp2[1:])/2.
        #obtain cell sizes
        self.delta_r_sp2 = self.grid_edge_sp2[1:] - self.grid_edge_sp2[:-1]
        #construct cell centers from this
        self.grid_sp = (self.grid_edge_sp[:-1] + self.grid_edge_sp[1:])/2.
        #obtain cell sizes
        self.delta_r_sp = self.grid_edge_sp[1:] - self.grid_edge_sp[:-1]
        print 'test grids'
        print 'sp1', self.grid_edge_sp1
        print 'sp2', self.grid_edge_sp2

        #create a fipy mesh for visualization and fipy computation
        self.mesh_fiber = CylindricalGrid1D(dx=tuple(self.delta_r))
        self.mesh_fiber.periodicBC = False
        self.mesh_fiber_sp1 = CylindricalGrid1D(dx=tuple(self.delta_r_sp1))
        self.mesh_fiber_sp1.periodicBC = False
        self.mesh_fiber_sp2 = CylindricalGrid1D(dx=tuple(self.delta_r_sp2))
        self.mesh_fiber_sp2.periodicBC = False
        self.mesh_fiber_sp2 = self.mesh_fiber_sp2 + ((self.l,),)
        self.mesh_fiber_sp = CylindricalGrid1D(dx=tuple(self.delta_r_sp))
        self.mesh_fiber_sp.periodicBC = False

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
        self.initial_c1_sp = sp.empty(self.tot_edges_sp-1, float)
        self.__tmp_sp = sp.empty(self.tot_edges_sp-1, float)
        self.ind_l_sp1 = None
        for i, grid in enumerate(self.grid_sp1):
            if grid < self.l:
                self.initial_c1_sp1[i] = 1.
            else:
                self.ind_l_sp1 = i
                break
        if self.ind_l_sp1 is None:
            raise ValueError, 'grid domain 1 must have grid cell > l'
        self.ind_s_sp2 = 1
        for i, grid in enumerate(self.grid_sp2):
            if grid < self.s:
                pass
            else:
                self.ind_s_sp2 = i
                break
        self.initial_c1_sp[:self.ind_l_sp1] = self.initial_c1_sp1[:self.ind_l_sp1]
        self.initial_c1_sp[self.ind_l_sp1:] = self.initial_c1_sp2[:]
        self.initial_w1_sp1 = self.initial_c1_sp1 * self.grid_sp1
        self.initial_w1_sp2 = self.initial_c1_sp2 * self.grid_sp2
        if self.verbose:
            m1 = self.calc_mass(self.initial_c1_sp1, split=1)
            m2 = self.calc_mass(self.initial_c1_sp2, split=2)
            print 'initial mass split = ', m1 , ' +' , m2, '=', m1+m2

    def calc_mass(self, conc_r, split=0):
        """calculate the mass of component present given value in cell center
        This is given by 2 \pi int_r1^r2 C(r)r dr
        
        conc_r: concentration in self.grid
        """
        if split == 1:
            grid = self.grid_edge_sp1
        elif split == 2:
            grid = self.grid_edge_sp2
        else:
            grid = self.grid_edge
        return sp.sum(conc_r * (sp.power(grid[1:], 2) - 
                                sp.power(grid[:-1], 2)) 
                        ) * sp.pi 

    def calc_mass_overlap_sp1(self, conc_r):
        mass_overlap = sp.sum(conc_r[self.ind_l_sp1:] \
                 * (sp.power(self.grid_edge_sp1[self.ind_l_sp1+1:], 2) - 
                    sp.power(self.grid_edge_sp1[self.ind_l_sp1:-1], 2)) 
                    ) * sp.pi
        return mass_overlap

    def calc_mass_overlap_sp2(self, conc_r):
        mass_overlap = sp.sum(conc_r[:self.ind_s_sp2] \
                 * (sp.power(self.grid_edge_sp2[1:self.ind_s_sp2+1], 2) - 
                    sp.power(self.grid_edge_sp2[:self.ind_s_sp2], 2)) 
                    ) * sp.pi
        return mass_overlap

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
        (This is the flux of w over radius! So times 2 pi R for full flux w, or)
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
        #add the source term in the overlap part
        diff_w_t[self.ind_l_sp1:] += self.grid_sp1[self.ind_l_sp1:] * \
                        self.source_sp1

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
        #add the source term in the overlap part
        #print 'w_t bef', diff_w_t[:self.ind_s_sp2]
        diff_w_t[:self.ind_s_sp2] += self.grid_sp2[:self.ind_s_sp2] * \
                self.source_sp2 
        #print 'w_t aft', diff_w_t[:self.ind_s_sp2]
        #raw_input('')

    def solve_odes_init(self):
        """
        Initialize the cvode solver
        """
        if not HAVE_ODES:
            raise Exception, 'Not possible to solve with given method, scikits.odes not available'
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        self.step_old_time_sp1 = self.initial_t
        self.step_old_time_sp2 = self.initial_t
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
        self.solver_sp1.init_step(self.step_old_time_sp1, self.step_old_sol_sp1)
        self.solver_sp2 = sc_ode('cvode', self.f_conc1_odes_sp2,
                             max_steps=50000, lband=1, uband=1)
        self.solver_sp2.init_step(self.step_old_time_sp2, self.step_old_sol_sp2)

    def do_solve_step(self, stoptime):
        """
        Solve up to time t. This does:
           1. solve the global model up to t
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
        self.plotting = False
        if self.plotevery:
            if self.viewerplotcount == 0:
                self.solution_view.setValue(conc)
                self.viewer.axes.set_title('time %s' %str(time))
                self.plotting = True
                self.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'conc%s.png' % str(int(10*time)))
            self.viewerplotcount += 1
            self.viewerplotcount = self.viewerplotcount % self.plotevery
        return time, conc

    def do_solve_step_sp1(self, stoptime):
        """
        Solve up to time t. This does:
           2. solve split part 1 up to t
        """
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time_sp1
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            #solve of global model
            time, conc = self.do_step_odes_sp1(t)

        if self.plotting:
            self.solution_view_sp1.setValue(conc)
            self.viewer_sp1.axes.set_title('time %s' %str(time))
            self.viewer_sp1.plot(filename=utils.OUTPUTDIR + os.sep + 'conc_sp1_%s.png' % str(int(10*time)))
        return time, conc

    def do_solve_step_sp2(self, stoptime):
        """
        Solve up to time t. This does:
           4. solve split part 2 up to t
        """
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time_sp2
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            #solve of global model
            time, conc = self.do_step_odes_sp2(t)

        if self.plotting:
            self.solution_view_sp2.setValue(conc)
            self.viewer_sp2.axes.set_title('time %s' %str(time))
            self.viewer_sp2.plot(filename=utils.OUTPUTDIR + os.sep + 'conc_sp2_%s.png' % str(int(10*time)))
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
            self.solver_sp1.set_options(tstop=stoptime)
            self.solver_sp2.set_options(tstop=stoptime)
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

    def do_step_odes_sp1(self, stoptime):
        """
        Solve the fibermodel part 1 up to stoptime, continuing from the present
        state, return time and r * concentration after step
        It is needed that run_init and solve_init method have been called
        before calling this method.
        
        A needreinit has to be set in the global model

        Return: concentration over the grid
        """
        assert stoptime > self.step_old_time_sp1, "%f > %f" % (stoptime, self.step_old_time_sp1)
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time_sp1
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            flag, realtime = self.solver_sp1.step(t, self.ret_y_sp1)
            if flag < 0:
                raise Exception, 'could not find solution'
        
        self.step_old_time_sp1 = realtime
        self.step_old_sol_sp1 = self.ret_y_sp1
        assert np.allclose(realtime, stoptime, atol=1e-6, rtol=1e-6)
        return realtime, self.ret_y_sp1 / self.grid_sp1

    def do_step_odes_sp2(self, stoptime):
        """
        Solve the fibermodel part 2 up to stoptime, continuing from the present
        state, return time and r * concentration after step
        It is needed that run_init and solve_init method have been called
        before calling this method.
        
        A needreinit has to be set in the global model

        Return: concentration over the grid
        """
        assert stoptime > self.step_old_time_sp2, "%f > %f" % (stoptime, self.step_old_time_sp2)
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time_sp2
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            flag, realtime = self.solver_sp2.step(t, self.ret_y_sp2)
            if flag < 0:
                raise Exception, 'could not find solution'
        
        self.step_old_time_sp2 = realtime
        self.step_old_sol_sp2 = self.ret_y_sp2
        assert np.allclose(realtime, stoptime, atol=1e-6, rtol=1e-6)
        return realtime, self.ret_y_sp2 / self.grid_sp2

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
    
    def plot_sp(self, time, conc_sp1, conc_sp2):
        """
        Combines plot of splitted parts to create a single plot
        """
        if self.plotting:
            self.__tmp_sp[:self.ind_l_sp1] = conc_sp1[:self.ind_l_sp1]
            self.__tmp_sp[self.ind_l_sp1:] = conc_sp2[:]
            self.solution_view_sp.setValue(self.__tmp_sp)
            self.viewer_sp.axes.set_title('time %s' %str(time))
            self.viewer_sp.plot(filename=utils.OUTPUTDIR + os.sep + 
                                        'conc_sp_%s.png' % str(int(10*time)))

    def run(self, wait=False, output=False):
        self.run_init()
        if not self.initialized:
            self.solve_init()
        
        self.solution_view = CellVariable(name = "fiber conc single", 
                mesh = self.mesh_fiber,
                value = self.conc1[0][:])
        self.viewer =  Matplotlib1DViewer(vars = self.solution_view, 
                            datamin=0., 
                            datamax=1.2 * self.conc1[0].max())
        self.viewer.axes.set_title('time 0.0')
        self.viewer.plot()
        self.viewerplotcount = 1
        self.solution_view_sp1 = CellVariable(name = "fiber conc part 1", 
                mesh = self.mesh_fiber_sp1,
                value = self.conc1_sp1[0][:])
        self.viewer_sp1 =  Matplotlib1DViewer(vars = self.solution_view_sp1, 
                            datamin=0., 
                            datamax=1.2 * self.conc1[0].max())
        self.viewer_sp1.axes.set_title('time 0.0')
        self.viewer_sp1.plot()
        self.solution_view_sp2 = CellVariable(name = "fiber conc part 2", 
                mesh = self.mesh_fiber_sp2,
                value = self.conc1_sp2[0][:])
        self.viewer_sp2 =  Matplotlib1DViewer(vars = self.solution_view_sp2, 
                            datamin=0., 
                            datamax=1.2 * self.conc1[0].max())
        self.viewer_sp2.axes.set_title('time 0.0')
        self.viewer_sp2.plot()
        self.solution_view_sp = CellVariable(name = "fiber conc split", 
                mesh = self.mesh_fiber_sp,
                value = self.initial_c1_sp)
        self.viewer_sp =  Matplotlib1DViewer(vars = self.solution_view_sp, 
                            datamin=0., 
                            datamax=1.2 * self.conc1[0].max())
        self.plot_sp(0., self.conc1_sp1[0], self.conc1_sp2[0])
        
        """
        Solve up to time t. This does:
           1. solve global model up to time t
           2. solve split part 1 up to t
           3. set correct flux/source term for part 2
           4. solve split part 2 up to t
           5. set correct flux/source term for part 1
        """
        self.source_sp1 = 0
        self.source_sp2 = 0
        mass_overlap_sp1 = self.calc_mass_overlap_sp1(self.initial_c1_sp1)
        mass_overlap_sp2 = self.calc_mass_overlap_sp2(self.initial_c1_sp2)
        newmass_sp1 = 0.
        newmass_sp2 = 0.
        for ind, t in enumerate(self.times[1:]):
            #step 1.
            rt, conc = self.do_solve_step(t)
            self.conc1[ind+1][:] = conc[:]
            self.fiber_surface[ind+1] = conc[-1]
            self.flux_at_surface[ind+1] = 0.
            #no we do same for the split function
            #step 2. solve in domain 1
            mass_overlap_sp1 = self.calc_mass_overlap_sp1(self.conc1_sp1[ind])
            rt, conc_sp1 = self.do_solve_step_sp1(t)
            self.conc1_sp1[ind+1][:] = conc_sp1[:]
            #step 3. mass that went to domain 2
            mass_overlap_sp1new = self.calc_mass_overlap_sp1(conc_sp1)
            newmass = mass_overlap_sp1new - mass_overlap_sp2
            mass_overlap_sp1 = mass_overlap_sp1new
            # source term per area per second
            self.source_sp2 = newmass / self.volume_overlap / self.delta_t
            #step 4. solve in domain 2
            rt, conc_sp2 = self.do_solve_step_sp2(t)
            self.conc1_sp2[ind+1][:] = conc_sp2[:]
            #step 5. mass that was extracted from domain 1
            mass_overlap_sp2new = self.calc_mass_overlap_sp2(conc_sp2)
            newmass_sp2 = mass_overlap_sp2new - mass_overlap_sp1
            mass_overlap_sp2 = mass_overlap_sp2new
            # source term per area per second
            self.source_sp1 = (newmass_sp2) / self.volume_overlap / self.delta_t
            self.plot_sp(rt, self.conc1_sp1[ind+1], self.conc1_sp2[ind+1])

        if output:
            self.dump_solution()
        if self.verbose:
            print 'end mass = ', self.calc_mass(conc)
            m1 = self.calc_mass(conc_sp1, split=1)
            mover = self.calc_mass_overlap_sp1(conc_sp1)
            m2 = self.calc_mass(conc_sp2, split=2)
            print 'end mass split = ', m1-mover , ' +' , m2, '=', m1-mover+m2
        if wait:
            raw_input("Finished domainsplit run")

    def __del__(self):
        #remove the memory
        if self.solver:
            print 'del self.solver'
            del self.solver
        self.solver = None

""" function to start this sandbox test
Typical calling sequence 
PYTHONPATH=/home/benny/git/odes/build/lib.linux-x86_64-2.7/:/home/benny/git/stickproject/stick python domainsplit.py
"""
if __name__ == '__main__':
    dsm = DomainsplitTestModel()
    dsm.run(wait=True)