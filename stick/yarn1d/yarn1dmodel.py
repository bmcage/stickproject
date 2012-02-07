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

HAVE_ODES = False
try:
    from scikits.odes import ode as sc_ode
    HAVE_ODES = True
except:
    print 'Could not load scikits.odes, odes solver not available'


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
import yarn.config as conf
from fiber.config import FiberConfigManager
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
        self.verbose = self.cfg.get('general.verbose')
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0., self.time_period, num=self.steps+1)
        self.delta_t = self.times[1] - self.times[0]
        if self.verbose:
            print "Timestep used in yarn1d model:", self.delta_t
        
        self.diff_coef = self.cfg.get('diffusion.diffusion_coeff')
        self.init_conc_func = eval(self.cfg.get('initial.init_conc'))
        
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        self.nr_models = self.cfg.get('fiber.number_type')
        assert self.nr_models == len(self.blend) == len(self.cfg.get('fiber.fiber_config'))

        #construct the config for the fibers
        self.cfg_fiber = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            self.cfg_fiber.append(FiberConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.time_period)
            if self.cfg_fiber[-1].get("time.dt") > self.cfg.get("time.dt"):
                self.cfg_fiber[-1].set("time.dt", self.cfg.get("time.dt"))
            #we need stepwize solution, we select cvode
            self.cfg_fiber[-1].set("general.method", 'FVM')
            self.cfg_fiber[-1].set("general.submethod", 'cvode_step')
            #we check that boundary is transfer or evaporation
            bty = self.cfg_fiber[-1].get("boundary.type_right")
            if bty not in ['evaporation', 'transfer']:
                raise ValueError, 'Boundary type for a fiber should be evaporation or transfer'
            if self.verbose:
                print 'NOTE: Fiber has boundary out of type %s' %  bty

        #some memory
        self.cache_index_t_yarn = 0
        self.cache_index_t_fiber = [0] * self.nr_models
        
        #Initialize the tortuosity
        self.tortuosity= self.cfg.get('yarn.tortuosity')
        self.boundary_transf_right = self.cfg.get('boundary.transfer_conc1')
        self.evap_equilibrium = self.cfg.get('boundary.evap_equilibrium')
        self.nr_fibers = self.cfg.get('fiber.number_fiber')
        
        self.plotevery = self.cfg.get("plot.plotevery")
        self.writeevery = self.cfg.get("plot.writeevery")

        self.initialized = False

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh!
        
        grid: the space position of each central point for every cell element (r-coordinate);
        """
        self.beginning_point = 0 #center of the yarn, r=0, with r the distance from center yarn.
        self.end_point = self.cfg.get('domain.yarnradius')
        self.nr_edge = self.cfg.get('domain.n_edge')
        #we now construct the full edge grid
        self.grid_edge = sp.linspace(self.beginning_point, self.end_point, self.nr_edge)
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]
        grid_square = np.power(self.grid_edge, 2)
        self.delta_rsquare = grid_square[1:] - grid_square[:-1]
        
        #create fiber models as needed: one per fibertype and per cell in the yarn model
        self.fiber_models = [0] * (self.nr_edge - 1)
        self.fiber_mass = np.empty((self.nr_edge - 1, self.nr_models), float)
        self.source_mass = np.empty((self.nr_edge - 1, self.nr_models), float)
        self.source = np.empty(self.nr_edge - 1, float)
        for ind in range(self.nr_edge-1):
            self.fiber_models[ind] = []
            for cfg in self.cfg_fiber:
                self.fiber_models[ind].append(FiberModel(cfg))

        #create cylindrical 1D grid over domain for using fipy to view.
        if self.plotevery:
            self.mesh_yarn = CylindricalGrid1D(dr=tuple(self.delta_r))
            self.mesh_yarn.periodicBC = False
            self.mesh_yarn = self.mesh_yarn + (self.beginning_point,)

        print 'mesh', self.grid_edge, ', delta_r', self.delta_r

    def initial_yarn1d(self):
        """ initial concentration over the domain"""
        self.init_conc = sp.ones(self.nr_edge-1, float)
        for ind, r in enumerate(self.grid):
            self.init_conc[ind] = self.init_conc_func(r)

    def out_conc(self, cellnr, t):
        """
        return the concentration of compound in the void zone of cell cellnr at
        time t
        """
        timenowyarn = self.times[self.tstep]
        if t >= timenowyarn:
            return self.conc1[self.tstep, cellnr]
        raise ValueError, 'out concentration should only be requested at a later time'

    def solve_fiber_init(self):
        """
        Solve the diffusion process for a repellent on the fiber at radial position r in the yarn.
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform.
        """        
##        self.nr_timesteps = np.empty((self.nr_models),int)
##        self.timesteps = [0]*self.nr_models
##        self.fiber_surface = [0] * self.nr_models
        
        for ind, models in enumerate(self.fiber_models):
            for type, model in enumerate(models):
                model.run_init()
                model.solve_init()
                #rebind the out_conc method to a call to yarn1d
                model.yarndata = ind
                model.out_conc = lambda t, data: self.out_conc(data, t)
                self.fiber_mass[ind, type] = model.calc_mass(model.initial_c1)

    def do_fiber_step(self, stoptime):
        """
        Solve the diffusion process on the fiber up to stoptime, starting
        from where we where last. 
        The flux is the BC: S*h(C_equi - C_yarn(t))*H(C-C_b,C_equi-C_yarn(t))
        """
        for ind, models in enumerate(self.fiber_models):
            for type, model in enumerate(models):
                time, result = model.do_step(stoptime, needreinit=False)
                tmp = model.calc_mass(result)
                self.source_mass[ind, type] = self.fiber_mass[ind, type] - tmp
                self.fiber_mass[ind, type] = tmp

    def _set_bound_flux(self, flux_edge, conc_r):
        """
        Method that takes BC into account to set flux on edge
        flux here is the flux per radial
        Data is written to flux_edge, conc_r contains solution in the cell centers
        """
        flux_edge[0] = 0.
        # tranfer flux, in x: flux_x  = tf * C, so radially per radial a flux
        #  flux_radial = 2 Pi * radius * flux_x / 2 * Pi
        flux_edge[-1] = self.boundary_transf_right * conc_r[-1] * self.grid_edge[-1]

    def set_source(self, timestep):
        """
        Method to calculate the radial source term
        Per radial we have the global equation 
           \partial_t (r C) = \partial_r (D/tau) r \partial_r C + r Source
        where Source is amount per time per volume released/absorbed
        This equation is integrated over a shell, and we determine 
            d_t w, with w = rC, where the sourceterm is \int_{r_i}^{r_{i+1}} r Source
        and is the term here calculated and stored in self.source
        As we assume Source constant over a shell, we have
            sourceterm = Source * \Delta r_i^2 / 2
        self.source_mass contains per shell how much mass was released in 
        previous step by a fiber. Suppose this mass is A. 
        We determine how many fibers there are radially, multiply this with A
        and divide by volume and timestep to obtain Source
        """
        
        #nrf is number of fibers in the shell at that grid position
        # per radial
        nrf_shell = (self.delta_rsquare/self.end_point**2 * self.nr_fibers 
                     / 2. / np.pi)
        for ind, pos in enumerate(self.grid):
            self.source[ind] = 0.
            for type, blend in enumerate(self.blend):
                #nrf is number of fibers of blend in the shell at that grid position
                # per radial
                self.source[ind] += self.source_mass[ind, type] * nrf_shell[ind] * blend
            self.source[ind] /= timestep
        self.source *= self.delta_rsquare / 2.

##    def get_source(self, t):
##        #find right index of interval for t in model.times and in self.times
##        if self.times[self.cache_index_t_yarn] <= t and \
##                t< self.times[self.cache_index_t_yarn+1] :
##            self.index_t_yarn = self.cache_index_t_yarn
##        else:
##            self.index_t_yarn = None
##            i = max([self.cache_index_t_yarn - 1,0])
##            while i < self.steps-1:
##                if self.times[i] <= t and t< self.times[i+1] :
##                    self.index_t_yarn = i
##                    break
##                i += 1
##            if self.index_t_yarn is None:
##                #backward in time, so reducing timestep it seems
##                i = self.cache_index_t_yarn-1
##                while i > 0:
##                    if self.times[i] <= t and t< self.times[i+1] :
##                        self.index_t_yarn = i
##                        break
##                    i -= 1
##            if self.index_t_yarn is None:
##                #no interval found
##                if t > self.times[-1]:
##                    self.index_t_yarn = self.steps-1
##                    print 'time over endtime', t, '>', self.times[-1], ", set index t to max", self.index_t_yarn
##                else:
##                    self.index_t_yarn = self.steps-1
##                    print "endtime,", t, self.times, ", set index t to max", self.index_t_yarn
##                    
##                    #raise exception, 'something wrong'
##        self.cache_index_t_yarn = self.index_t_yarn
##        
##        #the same now for the time of the fiber models
##        nr = 0
##        self.index_t_fiber = [0]*self.nr_models
##        
##        while nr < self.nr_models:
##            if self.timesteps[nr][self.cache_index_t_fiber[nr]] <= t and \
##                t< self.timesteps[nr][self.cache_index_t_fiber[nr]+1] :
##                self.index_t_fiber[nr] = self.cache_index_t_fiber[nr]
##                #print "interval found in loop 1"
##            else:
##                #self.index_t_fiber[nr] = None
##                i = max([self.cache_index_t_fiber[nr] - 1,0])
##                while i < self.nr_timesteps[nr] - 1:
##                        if self.timesteps[nr][i] <= t and t< self.timesteps[nr][i+1] :
##                            self.index_t_fiber[nr] = i
##                            break
##                        #print 'i', i, "index_t", self.index_t_fiber, "interval found in loop 2"
##                        i += 1                                    
##                if self.index_t_fiber[nr] is None:
##                    #backward in time, so reducing timestep it seems
##                    i = self.cache_index_t_fiber[nr]-1
##                    while i > 0:
##                        if self.timesteps[nr][i] <= t and t< self.timesteps[nr][i+1] :
##                            self.index_t_fiber[nr] = i
##                            break
##                        #print "interval found in loop 3"
##                        i -= 1                    
##                if self.index_t_fiber[nr] is None:
##                    #no interval found
##                    if t > self.timesteps[nr][-1]:
##                        print 'ERROR: time over endtime', t, '>', self.timesteps[nr][-1]
##                        self.index_t_fiber[nr] = self.nr_timesteps[nr] - 1
##                        break
##                    else:
##                        print nr, t, self.timesteps
##                        raise Exception, 'something wrong'
##            self.cache_index_t_fiber[nr] = self.index_t_fiber[nr]
##            nr+=1
##
##        #source term is n*Cf(R,r_i+,t)/2pi=(m*delta(r**2)_i/Ry**2)*Cf(R,r_i+,t)/2pi with n the number of fibers in a shell,
##        #m the number of fibers per yarn.
##        grid_square = np.power(self.grid_edge, 2)
##        self.delta_rsquare = grid_square[1:] - grid_square[:-1]
##        n = self.nr_fibers*self.delta_rsquare/(self.end_point**2)
##
##        fibersurf = 0.
##        for ind, blend in enumerate(self.blend):
##            #print "size fiber_surf", size(self.fiber_surface), "size index_t", size(self.index_t_fiber), "ind", ind, "t", self.index_t_fiber[ind]
##            fiber_surf_t = self.fiber_surface[ind][self.index_t_fiber[ind]]  +\
##                    (self.fiber_surface[ind][self.index_t_fiber[ind]+1] 
##                      - self.fiber_surface[ind][self.index_t_fiber[ind]])\
##                    /(self.timesteps[ind][self.index_t_fiber[ind]+1]
##                       -self.timesteps[ind][self.index_t_fiber[ind]])\
##                    *(t-self.timesteps[ind][self.index_t_fiber[ind]])
##            fibersurf = fibersurf + fiber_surf_t * blend/100
##        self.source[self.index_t_yarn,:]=n*fibersurf/(2*np.pi)

    def f_conc1_ode(self, t, conc_r, diff_u_t):
        ##print 't, concr', t, conc_r
        grid = self.grid
        n_cellcenters = len(grid)
        #Initialize the flux rate on the edges
        flux_edge = self.__tmp_flux_edge
        self._set_bound_flux(flux_edge, conc_r)

        #calculate flux rate in each edge of the domain
        flux_edge[1:-1] = (2 * (self.diff_coef/self.tortuosity) *
                           self.grid_edge[1:-1]*(conc_r[1:]-conc_r[:-1])
                           /(self.delta_r[:-1]+self.delta_r[1:]))
        diff_u_t[:] = ( ((flux_edge[1:]-flux_edge[:-1]) + self.source[:])
                        / self.delta_r[:] 
                      )

        return diff_u_t

    def solve_ode_init(self):
        """
        Initialize the ode solver
        """
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        n_cells = len(self.init_conc)
        self.conc1 = np.empty((len(self.times), n_cells), float)
        self.ret_y = np.empty(n_cells, float)
        self.__tmp_flux_edge = sp.empty(n_cells+1, float)
        self.tstep = 0
        self.conc1[0][:] = self.init_conc[:]
        
        self.solver = sc_ode('cvode', self.f_conc1_ode,
                             max_steps=50000, lband=1, uband=1)
        self.solver.init_step(self.step_old_time, self.init_conc)
        self.initialized = True

    def do_ode_step(self, stoptime):
        """Solve the yarnmodel up to stoptime, continuing from the present
           state, return the time, concentration after step
        """
        self.solver.set_tcrit(tcrit=stoptime)
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'

        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        while compute:
            t = self.step_old_time + self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            flag, realtime = self.solver.step(t, self.ret_y)
            if flag < 0:
                raise Exception, 'could not find solution, flag %d' % flag
        assert np.allclose(realtime, stoptime), "%f %f" % (realtime, stoptime)
        self.tstep += 1
        self.conc1[self.tstep][:] = self.ret_y[:]/self.grid[:]
        self.step_old_time = stoptime
        return stoptime, self.ret_y

    def view_sol(self, times, conc):
        """
        Show the solution in conc with times.
        conc[i][:] contains solution at time times[i]
        """
        if self.plotevery:
            self.solution_view = CellVariable(name = "Yarn radial concentration", mesh = self.mesh_yarn, value = conc[0][:])
            self.viewer =  Matplotlib1DViewer(vars = self.solution_view, datamin=0., datamax=conc.max()+0.20*conc.max())
            self.viewerplotcount = 0
            self.viewerwritecount = 0
            for time, con in zip(times, conc):
                self.solution_view.setValue(con)
                if self.viewerplotcount == 0:
                   self.viewer.plot()
                   self.viewer.axes.set_title('time %s' %str(time))
                self.viewerplotcount += 1
                self.viewerplotcount = self.viewerplotcount % self.plotevery
                if self.writeevery:
                    if self.viewerwritecount == 0:
                        self.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'yarnconc%08.4f.png' % time)
                    self.viewerwritecount += 1
                    self.viewerwritecount = self.viewerwritecount % self.writeevery

    def run(self, wait=False): 
        self.create_mesh()
        self.initial_yarn1d()
        self.solve_fiber_init()
        if not self.initialized:
            self.solve_ode_init()
        
        print 'Start mass of DEET per grid cell per fiber type'
        for ind, masses in enumerate(self.fiber_mass):
            print 'cell', ind,
            for mass in masses:
                print mass, ' - ',
            print ' '

        for t in self.times[1:]:
            #print 'solving t', t
            self.do_fiber_step(t)
            self.set_source(t-self.step_old_time)
            self.do_ode_step(t)

        print 'Final mass of DEET per grid cell per fiber type'
        for ind, masses in enumerate(self.fiber_mass):
            print 'cell', ind,
            for mass in masses:
                print mass, ' - ',
            print ' '

        self.view_sol(self.times, self.conc1)

        if wait:
            raw_input("Finished fiber1d run")
