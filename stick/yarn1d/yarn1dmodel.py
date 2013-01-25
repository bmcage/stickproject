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
import numpy as np

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
import stick.const as const
import stick.lib.utils.utils as utils
import stick.yarn.config as conf
from stick.fiber.config import FiberConfigManager
from stick.fiber1d.fibermodel import FiberModel

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
    An overlap region outside of the yarn can be added for multiscale simulations
    Fipy solve the transient diffusion problem in the whole domain
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        """
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = np.linspace(0., self.time_period, num=self.steps+1)
        self.delta_t = self.times[1] - self.times[0]
        if self.verbose:
            print "Timestep used in yarn1d model:", self.delta_t
        
        self.diff_coef = self.cfg.get('diffusion.diffusion_coeff')
        self.init_conc_func = eval(self.cfg.get('initial.init_conc1d'))
        
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        self.blend = [x/100. for x in self.blend]
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
        self.step_old_time = None
        self.step_old_sol = None
        
        #Initialize the tortuosity
        self.tortuosity= self.cfg.get('yarn.tortuosity')
        #use the area function for calculating porosity
        self.prob_area = eval(self.cfg.get('fiber.prob_area'))
        # boundary data
        self.bound_type = conf.BOUND_TYPE[self.cfg.get('boundary.type_right')]
        self.boundary_conc_out = self.cfg.get('boundary.conc_out')
        self.boundary_D_out = self.cfg.get('boundary.D_out')
        self.boundary_dist = self.cfg.get('boundary.dist_conc_out')
        self.boundary_transf_right = self.cfg.get('boundary.transfer_coef')
        self.nr_fibers = self.cfg.get('fiber.number_fiber')
        
        self.plotevery = self.cfg.get("plot.plotevery")
        self.writeevery = self.cfg.get("plot.writeevery")
        
        #allow a multiscale model to work with a source in overlap zone
        self.source_overlap = 0.

        self.initialized = False

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh!
        
        grid: the space position of each central point for every cell element (r-coordinate);
        """
        self.end_point = self.cfg.get('domain.yarnradius')
        self.nr_edge = self.cfg.get('domain.n_edge')
        self.nr_cell = self.nr_edge - 1
        self.use_extend = self.cfg.get("domain.useextension")
        self.fiberlayout_method = self.cfg.get('domain.fiberlayout_method')
        self.areaextend = 0.
        if self.use_extend:
            self.end_extend = self.end_point + \
                self.cfg.get('domain.extensionfraction') * self.end_point
            self.nr_edge_extend = max(2, 
                    int(self.nr_edge*self.cfg.get('domain.extensionfraction')))
            self.areaextend = np.pi * (self.end_extend**2 - self.end_point**2)
        else:
            self.end_extend = self.end_point
            self.nr_edge_extend = 1
        #we now construct the full edge grid
        self.nr_edge_tot = self.nr_edge + self.nr_edge_extend - 1
        self.nr_cell_tot = self.nr_edge_tot - 1
        self.grid_edge = np.empty(self.nr_edge_tot, float)
        self.grid_edge[:self.nr_edge] = np.linspace(0., 
                                    self.end_point, self.nr_edge)
        self.grid_edge[self.nr_edge:] = np.linspace(self.end_point,
                                    self.end_extend, self.nr_edge_extend)[1:]
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]
        grid_square = np.power(self.grid_edge, 2)
        self.delta_rsquare = grid_square[1:] - grid_square[:-1]
        
        #nrf is number of fibers in the shell at that grid position
        # per radial
        self.nrf_shell = (self.delta_rsquare[:self.nr_cell]\
                                / (self.end_point**2) * self.nr_fibers)
        #create fiber models as needed: one per fibertype and per cell in the yarn model
        self.fiber_models = [0] * (self.nr_edge - 1)
        self.fiber_mass = np.empty((self.nr_edge - 1, self.nr_models), float)
        self.source_mass = np.empty((self.nr_edge - 1, self.nr_models), float)
        self.source = np.zeros(self.nr_edge_tot - 1, float)
        for ind in range(self.nr_edge-1):
            self.fiber_models[ind] = []
            for cfg in self.cfg_fiber:
                self.fiber_models[ind].append(FiberModel(cfg))
        
        #calculate the porosity as n=(pi Ry^2-nr_fibers pi Rf^2) / pi Ry^2
        #porosity in the yarn
        self.porosity = np.ones(self.nr_cell_tot, float)
        self.volfracfib = []  # volume fraction of the fiber types
        if self.fiberlayout_method == 'virtlocoverlap':
            value_from_areafunction = np.zeros(self.nr_cell, float)
            for i_porosity in range(len(self.prob_area)):
                function_area = self.prob_area[i_porosity]
                value_from_areafunction += function_area(self.grid[:])
                #plot the value from the function and the porosity value to check
                #the calculation for porosity
                plt.figure()
                plt.plot(self.grid[:], 1 - value_from_areafunction, '-', color = 'red')
                plt.plot(self.grid[:], 1 - function_area(self.grid[:]), '*')
                plt.xlabel('Yarn domain')
                plt.ylabel('value')
                plt.ylim(0., 1.0)
                plt.show()
            self.porosity[:self.nr_cell] = 1. - value_from_areafunction[:self.nr_cell]
            
        else:        
            for blend, model in zip(self.blend, self.fiber_models[0]):
                print 'fiberradius', model.radius(), 'yarnradius', self.end_point   
                self.volfracfib.append(
                        blend * self.nr_fibers *  np.power(model.radius(),2)
                                / np.power(self.end_point,2) )
                if np.sum(self.volfracfib)>1:
                    raise ValueError, 'porosity  negative, unrealistic number of fibers in yarn cross section, %f fibers per yarn * Rf^2/Ry^2 = %f' % (self.nr_fibers,np.sum(self.volfracfib))
                    raw_input()
            self.porosity[:self.nr_cell] = 1- np.sum(self.volfracfib)
            print 'porosity in yarn', self.porosity[:self.nr_cell], 
            #if self.porosity[:self.nr_cell]<0:
                # raise ValueError, 'porosity  negative'
        #create cylindrical 1D grid over domain for using fipy to view.
        if self.plotevery:
            self.mesh_yarn = CylindricalGrid1D(dr=tuple(self.delta_r[:self.nr_cell]))
            self.mesh_yarn.periodicBC = False
            self.mesh_yarn = self.mesh_yarn

        print 'mesh yarn', self.grid_edge, ', delta_r yarn', self.delta_r

    def initial_yarn1d(self):
        """ initial concentration over the domain"""
        self.init_conc = np.empty(self.nr_cell_tot, float)
        #zero to the outside
        for ind, r in enumerate(self.grid[:self.nr_cell]):
            self.init_conc[ind] = self.init_conc_func(r)
        self.init_conc[self.nr_cell:] = self.cfg.get('boundary.conc_out')

    def get_data(self, cellnr):
        index = cellnr
        return index

    def out_conc(self, data, t):
        """
        return the concentration of compound in the void zone of cell cellnr at
        time t
        """
        timenowyarn = self.step_old_time
        if t >= timenowyarn:
            #return data
            return self.step_old_sol[data]
        raise ValueError, 'out concentration should only be requested at a later time'

    def solve_fiber_init(self):
        """
        Solve the diffusion process for a repellent on the fiber at radial position r in the yarn.
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform.
        """
        for ind, models in enumerate(self.fiber_models):
            for type, model in enumerate(models):
                model.run_init()
                model.solve_init()
                #rebind the out_conc method to a call to yarn1d
                model.set_userdata(self.get_data(ind))
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
        #avergage porosity on edge:
        if self.use_extend:
            porright = (self.porosity[self.nr_cell-1] + self.porosity[self.nr_cell]) / 2
        else:
            porright = (self.porosity[self.nr_cell-1] + 1.) / 2
        diffright = (self.diff_coef/self.tortuosity + self.boundary_D_out)/ 2
        if self.bound_type == conf.TRANSFER:
            # tranfer flux, in x: flux_x  = tf * C, so radially per radial a flux
            #  flux_radial = 2 Pi * radius * flux_x / 2 * Pi
            flux_edge[self.nr_edge-1] = self.boundary_transf_right * porright \
                    * conc_r[self.nr_cell-1] * self.grid_edge[self.nr_edge-1]
        elif self.bound_type == conf.DIFF_FLUX:
            # diffusive flux with the outside
            # flux radial = - D_out * (conc_out - yarn_edge_conc)/dist_conc_out * radius
            if self.use_extend:
                conright = conc_r[self.nr_cell]
                bcdist = (self.delta_r[self.nr_cell-1]+self.delta_r[self.nr_cell])/2.
            else:
                conright = self.boundary_conc_out
                bcdist = self.boundary_dist
            ##print 'fluxdiff', conright - conc_r[self.nr_cell-1], self.boundary_conc_out,  conc_r[self.nr_cell-1]
            flux_edge[self.nr_edge-1] = -(diffright * porright * 
                    (conright - conc_r[self.nr_cell-1]) 
                    / bcdist 
                    * self.grid_edge[self.nr_edge-1])
        if self.use_extend:
            #zero flux at right boundary
            flux_edge[-1] = 0.

    def calc_mass(self, conc):
        """
        calculate current amount of mass of volatile based on data currently
        stored
        """
        #first we calculate the mass in the void space:
        mass = np.sum(conc[:self.nr_cell] * 
                      (np.power(self.grid_edge[1:self.nr_edge], 2) -
                        np.power(self.grid_edge[:self.nr_edge-1], 2)) *
                      self.porosity[:self.nr_cell]
                     ) * np.pi
        #now we add the mass in the fibers
        for ind, pos in enumerate(self.grid[:self.nr_cell]):
            for type, blend in enumerate(self.blend):
                #nrf is number of fibers of blend in the shell at that grid position
                #print 'fiber mass', self.fiber_mass
                massfib = (self.fiber_mass[ind, type]
                            * self.nrf_shell[ind] * blend)
                #print 'mass fiber', self.fiber_mass[ind,type], 'nr fibers per shell', self.nrf_shell[ind]
                mass += massfib
        #print 'yarn conc', conc
        #print 'yarn totalmass',  mass, 'microgram'
        return mass

    def calc_mass_overlap(self, conc):
        """
        calculate current amount of mass of volatile in the overlap region
        based on data currently stored
        """
        #we calculate the mass in the void space:
        mass = np.sum(conc[self.nr_cell:] * 
                      (np.power(self.grid_edge[self.nr_edge:], 2) -
                        np.power(self.grid_edge[self.nr_edge-1:-1], 2)) *
                      self.porosity[self.nr_cell:]
                     ) * np.pi
        #print "porosity",self.porosity
        ##print 'yarn mass overlap', conc[self.nr_cell:], mass
        return mass

    def set_source(self, timestep):
        """
        Method to calculate the radial source term
        Per radial we have the global equation 
           \partial_t (n r C) = \partial_r (D/tau) r \partial_r (n C) + r Source
        where Source is amount per time per volume released/absorbed
        So Source = masssourceyarns / (V \Delta t)
        This equation is integrated over a shell and devided by n (the porosity), and we determine 
           n d_t w, with w = rC, where the sourceterm is \int_{r_i}^{r_{i+1}} r Source
        and is the term here calculated and stored in self.source
        As we assume Source constant over a shell by averaging out the mass over the area of a shell (nV), we have
            Source = masssourceyarns/(V \Delta t) * \Delta r_i^2 / 2
        self.source_mass contains per shell how much mass was released in 
        previous step by a fiber. Suppose this mass is M. 
        We determine how many fibers there are radially, multiply this with M
        and divide by volume V \delta t to obtain Source-concentration, since concentration is mass/volume time. 
        Afterwards we multiply this Source with \Delta r_i^2 / (2 n \Delta r)
        coming from the integration (int n d_t w gives the term n \Delta r).
        """
        for ind, pos in enumerate(self.grid_edge[:self.nr_cell]):
            self.source[ind] = 0.
            #V is the area of the shell
            V = np.pi*((pos+self.delta_r[ind])**2-pos**2)
            for type, blend in enumerate(self.blend):
                #nrf is number of fibers of blend in the shell at that grid position
                # per radial
                self.source[ind] += (self.source_mass[ind, type] 
                                        * self.nrf_shell[ind] * blend)
            #self.source[ind] /= timestep
            self.source[ind] /= V * timestep
        ##print 'source', self.source
        ## TODO Tine, what formula is this? Comment above says just multiply...
        ##self.source *= self.delta_rsquare / (2.*V*self.porosity)
        ## I think:
        ##    nrfibershell*blend* mass diff / Volume shell
        ## so only remains to divide by volume shell. HOWEVER, this must be 
        ## in for loop at level of cell, so pos. In essence, you have a 
        ## factor self.delta_rsquare / 2   more that I don't see
        ## Also: source must be per second, so divided by the timestep

    def f_conc1_ode(self, t, conc_r, diff_u_t):
        """
        Solving the radial yarn 1D diffusion equation: 
        
          \partial_t (rC) =  \partial_r (D/tau r \partial_r C) + Source * r
        
        with Source the conc amount per time unit added at r. 
        Solution is obtained by integration over a cell, so
        
           \delta r d_t (r C) = flux_right - flux_left + Source (\delta r^2 /2)
        
        so 
        
          d_t C = 1 / (r \delta r) * (flux_right - flux_left + Source (\delta r^2 /2) )
        """
        grid = self.grid
        n_cellcenters = len(grid)
        #Initialize the flux rate on the edges
        flux_edge = self.__tmp_flux_edge
        #set flux on edge 0, self.nr_edge-1 and -1
        self._set_bound_flux(flux_edge, conc_r)

        #calculate flux rate in each edge of the domain
        flux_edge[1:self.nr_edge-1] = -(2 * (self.diff_coef/self.tortuosity) *
            self.grid_edge[1:self.nr_edge-1] *
            (conc_r[1:self.nr_cell]-conc_r[:self.nr_cell-1])
            /(self.delta_r[:self.nr_cell-1]+self.delta_r[1:self.nr_cell])
            * (self.porosity[:self.nr_cell-1] + self.porosity[1:self.nr_cell])/2
            )
        if self.use_extend:
            # diffusion in the outside region
            flux_edge[self.nr_edge:-1] = -(2 * self.boundary_D_out *
            self.grid_edge[self.nr_edge:-1] *
            (conc_r[self.nr_cell+1:]-conc_r[self.nr_cell:-1])
            /(self.delta_r[self.nr_cell:-1]+self.delta_r[self.nr_cell+1:])
            * (self.porosity[self.nr_cell:-1] + self.porosity[self.nr_cell+1:])/2
            )
            ##print 'conc overlap', conc_r[self.nr_cell-1:]
            ##print 'flux overlap', flux_edge[self.nr_edge-1:]
        
        diff_u_t[:] = ((flux_edge[:-1]-flux_edge[1:])
                            / self.delta_r[:]/ self.porosity[:]
                    + self.source[:] / self.porosity[:] * self.delta_rsquare / 2
                      / self.delta_r
                      )
##        print 'source part', self.source[:] / self.porosity[:] * self.delta_rsquare / 2\
##                      / self.delta_r / self.grid[:]
##        print 'flux edges', flux_edge
        if self.use_extend and self.source_overlap:
            #porosity assumed 1 in extend!
            diff_u_t[self.nr_cell:] += (self.source_overlap
                            * self.delta_rsquare[self.nr_cell:] / 2
                            / self.delta_r[self.nr_cell:])
        diff_u_t[:] = diff_u_t[:] / self.grid[:]  # still division by r to move from w to C
        ##print 'u_t overlap', diff_u_t[self.nr_cell-1:]

    def solve_ode_init(self):
        """
        Initialize the ode solver
        """
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        n_cells = len(self.init_conc)
        self.conc1 = np.empty((len(self.times), n_cells), float)
        self.ret_y = np.empty(n_cells, float)
        self.__tmp_flux_edge = np.empty(n_cells+1, float)
        self.tstep = 0
        self.conc1[0][:] = self.init_conc[:]
        self.step_old_sol = self.conc1[0]
        
        self.solver = sc_ode('cvode', self.f_conc1_ode,
                             min_step_size=1e-8, rtol=1e-6, atol=1e-6, 
                             max_steps=50000, lband=1, uband=1)
        self.solver.init_step(self.step_old_time, self.init_conc)
        self.initialized = True

    def do_ode_step(self, stoptime):
        """Solve the yarnmodel up to stoptime, continuing from the present
           state, return the time, concentration after step
        """
        self.solver.set_options(tstop=stoptime)
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'

        flag, realtime = self.solver.step(stoptime, self.ret_y)
        if flag < 0:
            raise Exception, 'could not find solution, flag %d' % flag
        assert np.allclose(realtime, stoptime), "%f %f" % (realtime, stoptime)
        return stoptime, self.ret_y

    def do_yarn_init(self):
        """
        generic initialization needed before yarn can be solved
        """
        self.create_mesh()
        self.initial_yarn1d()
        self.solve_fiber_init()
        if not self.initialized:
            self.solve_ode_init()

    def do_yarn_step(self, stoptime):
        """
        Solve yarn up to time t. This does:
           1. solve the fiber up to t
           2. set correct source term for the yarn
           3. solve the yarn up to t
        """
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time
##        cm1 = self.calc_mass(self.step_old_sol)
##        cm2 = self.calc_mass_overlap(self.step_old_sol)
##        print 'mass in yarn domain', cm1, cm2, cm1 + cm2, 'src', \
##                self.source_overlap * self.areaextend * 0.1
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            self.do_fiber_step(t)
            self.set_source(t-self.step_old_time)
            realtime, self.step_old_sol = self.do_ode_step(t)
            cm1 = self.calc_mass(self.step_old_sol)
            cm2 = self.calc_mass_overlap(self.step_old_sol)
##            print 'mass in yarn domain', cm1, cm2, cm1 + cm2
##            print 'new sol', self.step_old_sol
##            raw_input('')
            self.tstep += 1
            self.step_old_time = t

        return realtime, self.step_old_sol

    def view_sol(self, times, conc):
        """
        Show the solution in conc with times.
        conc[i][:] contains solution at time times[i]
        """
        if self.plotevery:
            self.solution_view = CellVariable(name="Yarn radial concentration",
                        mesh=self.mesh_yarn, value=conc[0][:self.nr_cell])
            self.viewer =  Matplotlib1DViewer(vars=self.solution_view,
                                datamin=0., datamax=conc.max()+0.20*conc.max())
            self.viewerplotcount = 0
            self.viewerwritecount = 0
            for time, con in zip(times, conc):
                self.solution_view.setValue(con[:self.nr_cell])
                if self.viewerplotcount == 0 or (self.writeevery and 
                                            self.viewerwritecount == 0):
                    self.viewer.axes.set_title('time %s' %str(time))
                    if self.writeevery and self.viewerwritecount == 0:
                        #plot and savefig
                        self.viewer.plot(filename=utils.OUTPUTDIR + os.sep \
                                                + 'yarnconc%08.4f.png' % time)
                    else:
                        #only plot
                        self.viewer.plot()
                self.viewerplotcount += 1
                self.viewerplotcount = self.viewerplotcount % self.plotevery
                if self.writeevery:
                    self.viewerwritecount += 1
                    self.viewerwritecount = self.viewerwritecount % self.writeevery

    def run(self, wait=False):
        self.do_yarn_init()
        
        print 'Start mass of DEET per grid cell per fiber type'
        for ind, masses in enumerate(self.fiber_mass):
            print 'cell', ind,
            for mass in masses:
                print mass, ' - ',
            print ' '
        mc1 = self.calc_mass(self.step_old_sol)
        mc2 = self.calc_mass_overlap(self.step_old_sol)
        print 'Total mass in yarn', mc1, ', mass in overlap zone:', mc2, \
                'Sum', mc1 + mc2

        for t in self.times[1:]:
            #print 'solving t', t
            rt, rety = self.do_yarn_step(t)
            self.conc1[self.tstep][:] = self.ret_y[:]

        print 'Final mass of DEET per grid cell per fiber type'
        for ind, masses in enumerate(self.fiber_mass):
            print 'cell', ind,
            for mass in masses:
                print mass, ' - ',
            print ' '
        mc1 = self.calc_mass(self.step_old_sol)
        mc2 = self.calc_mass_overlap(self.step_old_sol)
        print 'Total mass in yarn', mc1, ', mass in overlap zone:', mc2, \
                'Sum', mc1 + mc2

        self.view_sol(self.times, self.conc1)

        if wait:
            raw_input("Finished yarn1d run")
