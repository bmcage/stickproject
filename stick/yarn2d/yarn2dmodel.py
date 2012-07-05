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
from copy import *
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

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
import yarn.config as conf
from mycorrection import MyDiffusionTermNoCorrection
from yarn2dgrid import Yarn2dGrid
from fiber.config import FiberConfigManager, CIRCLE
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
class Yarn2DModel(object):
    """
    Yarn2DModel is a special diffusion model for a single yarn which is composed 
    by a certain amount of fibers. Firstly, one cross-section of fiber is 
    generated. Then the uniform distribution function is used to generated the 
    postions of fibers in the yarn until the number reaches our requirement.
    After that a homogeneous mesh is generated by using Gmsh.
    Only diffusion processes in a single fiber and yarn are considered. 
    ODE of scipy solve the diffusion process in the layers of DEET and permithrine
    which are on the fiber
    Fipy solve the transient diffusion problem in the whole domain
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        """
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')
        #time data
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.delta_t = self.times[1] - self.times[0]
        if self.verbose:
            print "Timestep used in yarn1d model:", self.delta_t
        self.diffusion_DEET = self.cfg.get('diffusion.diffusion_conc')
        self.init_conc = eval(self.cfg.get('initial.init_conc'))
        
        self.Ry = self.cfg.get('domain.yarnradius')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        self.eps_value = self.cfg.get('fiber.eps_value')
        
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        self.nrtypefiber = self.cfg.get('fiber.number_type')
        self.fiber_edge_result = [0] * self.nrtypefiber
        self.boundary_diff_out = self.cfg.get('boundary.conc_out')
        assert self.nrtypefiber == len(self.blend) == len(self.cfg.get('fiber.fiber_config'))
        
        #construct cfg for fiber        
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
        #create fiber models
        self.fiber_models = []
        for cfg in self.cfg_fiber:
            self.fiber_models.append(FiberModel(cfg))

        #obtain some fiber data for quick reference
        self.Rf = []
        for fmodel in self.fiber_models:
            self.Rf.append(fmodel.radius())
        self.radius_fiber =  [self.scaleL * rad for rad in self.Rf]
        # boundary data
        self.bound_type = conf.BOUND_TYPE[self.cfg.get('boundary.type_right')]
        #self.boundary_conc_out = self.cfg.get('boundary.conc_out')
        self.boundary_D_out = self.cfg.get('boundary.D_out')
        self.boundary_dist = self.cfg.get('boundary.dist_conc_out')
        self.boundary_transf_right = self.cfg.get('boundary.transfer_coef')
        self.nr_fibers = self.cfg.get('fiber.number_fiber')
        self.plotevery = self.cfg.get("plot.plotevery")
        self.writeevery = self.cfg.get("plot.writeevery")
        self.writeoutcount = 0
        #take the value for the parameters of the fiber scale
        self.nr_edge = self.cfg.get("domain.n_edge")
        #some memory
        self.step_old_time = None
        self.step_old_sol = None
        
    def create_mesh(self):
        """
        Create a mesh for use in the model
        """
        #initiate the paramters for the fiber scale
        #fiber_models = [0] * (self.nr_edge - 1)
        self.fiber_mass = np.empty((self.nr_edge - 1, self.nrtypefiber), float)
        self.source_mass = np.empty((self.nr_edge - 1, self.nrtypefiber), float)
        #print "begin to create mesh"
        self.source = np.empty(self.nr_edge - 1, float)
        self.grid = Yarn2dGrid(self.cfg)
        self.mesh2d = self.grid.mesh_2d_generate(filename='yarn.geo',
                                regenerate=not self.cfg.get('general.read'))
        if self.cfg.onlymesh and self.verbose:
            print "Finished Mesh generation"
            
    def out_conc(self, cellnr, t):
        """
        return the concentration of compound in the void zone of cell cellnr at
        time t
        """
       #self.solve_fiber_init()
        timenowyarn = self.step_old_time
        if t >= timenowyarn:
            return self.step_old_sol[cellnr]
        raise ValueError, 'out concentration should only be requested at a later time'
    
    def initial_yarn2d(self):
        self.solve_fiber_init()
        datamax = self.cfg.get('plot.maxval')
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        cellCenter_x, cellCenter_y = self.mesh2d.getCellCenters()
        initialConc = []
        for i_x, i_y in zip(cellCenter_x, cellCenter_y):
            initialConc.append(self.init_conc(i_x, i_y))
        initialConc = sp.array(initialConc)
        self.conc = CellVariable(name = "Conc. Active Component", 
                    mesh = self.mesh2d, value = initialConc)
        self.viewer = None
        print 'the length for draw', len(self.conc)
        print 'the initial value', self.init_conc
        
        self.viewer = Viewer(vars = self.conc, title = 'Concentration of DEET', 
                            datamin = 0., datamax = 0.0005)
        self.viewer.plot()
        raw_input("take the example of the initial condition")
        self.viewerplotcount = 1

    def solve_fiber_init(self):
        """
        Initialize the solvers that do the fiber simulation
        """
        print 'the length of the self.fiber_models', len((self.fiber_models))
        
        for ind, model in enumerate(self.fiber_models):
            print 'the type of the models', ind
            model.run_init()
            model.solve_init()
            #rebind the out_conc to call yarn2D
            model.yarndata = ind
            model.out_conc = lambda t, data: self.out_conc(data, t)
            initial_concentration = model.init_conc[ind](1)
            print 'the concentration value', initial_concentration
            self.fiber_mass[ind] = model.calc_mass(initial_concentration)
            print 'the mass in the fiber', self.fiber_mass[ind]
            self.fiber_mass[ind] = model.calc_mass(model.initial_c1)

    def solve_fiber_step(self, stoptime):
        """
        Solve the diffusion process on the fiber up to time, starting
        from where we where last. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
##        for nyfib, model in enumerate(self.fiber_models):
##            #determine step neede to reach this time
##            step = time - model.step_old_time
##            self.fiber_edge_result[nyfib] = model.run_step(step)[-1]
        for nyfib, model in enumerate(self.fiber_models):
            #for type, model in enumerate(models):
            print 'begin to clculate step by step'
            print 'the stoptime value', stoptime
            raw_input('Enter')
            time_simulate, result = model.do_step(stoptime, needreinit = False)
            print 'begin to calculate in fiber', result
            tmp = model.calc_mass(result)
            print 'the mass in the void space', tmp
            self.source_mass[nyfib] = self.fiber_mass[nyfib] - tmp
            self.fiber_mass[nyfib] = tmp
            #self.fiber_edge_result[nyfib] = result[-1]

    def solve_single_component(self):
        """
        The DEET diffusion process is divided into two parts:
        (1) DEET diffuses through the layer containing permithrine on the fiber 
        and reaches the surface;
        (2) DEET begins to diffuse in the void space of  yarn
        So it means that the boundary condition of fiber has two steps:
        (1) When the DEET does not reach the surface of fiber, the inner and out
        boundaries are no flux boundary condition;
        (2) When the DEET reaches surface, the evaporation happens. So the boundaries 
        of fiber and yarn are changed to constant flux (Neumann boundary condition)
        """
        compute = True
        #input the trsient equation of diffusion        
        self.eq = TransientTerm() == MyDiffusionTermNoCorrection(
                                                coeff = self.diffusion_DEET)
        #get the position of the boundary faces
        xfc, yfc = self.mesh2d.getFaceCenters()
        xcc, ycc = self.mesh2d.getCellCenters()
        self.cell_volume = self.mesh2d.getCellVolumes()
        filepath3 = utils.OUTPUTDIR + os.sep + 'index_fiber.dat'
        filepath4 = utils.OUTPUTDIR + os.sep + 'yarn_out.dat'
        filepath5 = utils.OUTPUTDIR + os.sep + 'conc_fib.gz'
        self.fib_x = self.grid.x_position
        self.fib_y = self.grid.y_position
        self.all_fib_radius = self.grid.all_radius_fibers
        self.fiber_kind = self.grid.fiber_kind
        #we need to determine which fibers are of a specific kind
        self.fibers = []
        for nyfib in sp.arange(self.nrtypefiber):
            self.fibers += [self.fiber_kind == nyfib]
        #we need to determine which indexes of the fiber are of specific kind
        all_indexes = sp.arange(self.number_fiber)
        self.index_fiber = []
        for nyfib in sp.arange(self.nrtypefiber):
            self.index_fiber += [all_indexes[self.fibers[nyfib]]]

        #now determine boundaries
        face_in = ((self.mesh2d.getExteriorFaces()) & 
                    (sp.power(xfc,2) + sp.power(yfc,2) \
                     < (self.grid.radius_domain - self.grid.radius_boundlayer)**2))
        #we need to determine which nodes in the mesh are surface of 
        #of a certain fiber kind and create the inner boundary
        for form in self.grid.Rf_form:
            if form != CIRCLE:
                raise Exception, 'ERROR, not supported to define boundary conditions on ellipse'
        eps_fib = [val * self.eps_value for val in self.radius_fiber]
        self.int_bound = []
        for nyfib in sp.arange(self.nrtypefiber):
            print 'begin to find the yarn boundary'
            tmp = np.zeros(len(face_in), bool)
            print 'the length of each sub-loop', len(self.index_fiber[nyfib])
            for fib in self.index_fiber[nyfib]:
                #print 'define the boundary from the yarn'
                tmp = (((self.mesh2d.getExteriorFaces()) &
                    (sp.power(xfc - self.fib_x[fib], 2) + sp.power(yfc - self.fib_y[fib], 2)\
                    <= (self.radius_fiber[nyfib] + eps_fib[nyfib])**2)) | (tmp))
            self.int_bound.append(tmp)
        self.ext_bound = (~face_in) & (self.mesh2d.getExteriorFaces())
        
        #data structures to hold data
        conc_on_fib = sp.empty(self.nrtypefiber)
        flux_in_fib = sp.empty(self.nrtypefiber)
        conc_fib_out = [0] *  self.steps
        #self.initial_t = 0.
        value_face_out = np.zeros(len(self.ext_bound), float)
        determine_out = np.zeros(len(self.ext_bound), bool)
        self.record_conc = open(filepath4, "w")
        ## TODO: IMPROVE THIS BC
        extBC = FixedFlux(self.ext_bound, value = 0.0)
        t = self.step_old_time
        stop_time = self.time_period
        boundary_ex = self.boundary_diff_out
        each_step = 0
        print 'begin to calculate the DEET concentration'
        i = 0
        while compute:
            
            t += self.delta_t
            if t >= stop_time -self.delta_t / 100:
                t = stop_time
                compute = False
            print 'the time for the fiber_step', t
            self.solve_fiber_step(t)
            extBC = FixedFlux(self.ext_bound, value = boundary_ex)
            BCs = [extBC]
            for nyfib in sp.arange(self.nrtypefiber):
                conc_on_fib[nyfib] = self.fiber_edge_result[nyfib]
                flux_in_fib[nyfib] = (
                                    self.fiber_models[nyfib].boundary_transf_right * 
                                    conc_on_fib[nyfib]
                                    )
                BCs.append(FixedFlux(self.int_bound[nyfib], value = -flux_in_fib[nyfib]))
            conc_fib_out[each_step] = copy(conc_on_fib)
            self.initial_t = self.times[i+1]
            self.eq.solve(var = self.conc, boundaryConditions = tuple(BCs), 
                        dt = self.delta_t)
                        
##        for i in sp.arange(0, self.steps, 1, dtype=int):
##            #advance fiber solution one step
##            self.solve_fiber_step(self.times[i+1])
##            #update BC with new fiber solution
##            
##            extBC = FixedFlux(self.ext_bound, value = )
##            BCs = [extBC]
##            for nyfib in sp.arange(self.nrtypefiber):
##                conc_on_fib[nyfib] = self.fiber_edge_result[nyfib]
##                flux_in_fib[nyfib] = (
##                    self.fiber_models[nyfib].boundary_transf_right 
##                     * conc_on_fib[nyfib])
##                BCs.append(FixedFlux(self.int_bound[nyfib], value = -flux_in_fib[nyfib]))#[nyfib + 1] = FixedFlux(self.int_bound[nyfib], value = -flux_in_fib[nyfib])
##            conc_fib_out[i] = copy(conc_on_fib)
##            
##            self.initial_t = self.times[i+1]
##            self.eq.solve(var = self.conc, boundaryConditions = tuple(BCs), 
##                          dt = self.times[i+1] - self.times[i])
##            print 'Solution obtained at time = ', self.times[i+1]
            #self.tstep += 1
            self.step_old_time = t
            each_step += 1
            if self.writeoutcount == 0:
                conc_tot_each = self.conc.getValue()
                conc_face_ex = self.conc.getArithmeticFaceValue()
                for i_out in sp.arange(len(self.ext_bound)):
                    value_face_out[i_out] = float(conc_face_ex[i_out])
                    determine_out[i_out] = self.ext_bound[i_out]
                value_out_record = value_face_out[determine_out]
                conc1_average_out = np.sum(value_out_record) / len(value_out_record)
                if self.verbose:
                    print 'average concentration out', conc1_average_out
                boundary_ex = self._set_bound_flux(conc_face_ex, self.ext_bound)

                self.record_conc.write("%g, %g \n" %(self.times[i], conc1_average_out))
                print 'mass conservative with two fibers', self.cal_mass_void(conc_tot_each,
                                                    self.cell_volume) / self.scaleL
            self.writeoutcount += 1
            self.writeoutcount = self.writeoutcount % self.writeevery
                
            if self.viewer is not None:
                if self.viewerplotcount == 0:
                    self.viewer.plot()
                self.viewerplotcount += 1
                self.viewerplotcount = self.viewerplotcount % self.plotevery
            
        self.record_conc.close()
        dump.write({'time_step': self.times, 'conc_fib': conc_fib_out}, 
                    filename = filepath5, extension = '.gz')
        raw_input("Finished <press return>.....")
    
    def cal_mass_void(self, conc_void, cell_volume):
        """
        calculate the mass of component in the void space 
        conc_void: the concentration of materials in the void space
        This is given by: total_mass = sum(conc_void[] * cell_volume[])
        """
        return sp.sum(conc_void * cell_volume)
    
    def _set_bound_flux(self, conc_face_ex, ext_bound):
        """
        the method that takes BC into account to set the flux on edge
         
        """
        conc_face_ex = float(conc_face_ex)
        conc_face_edge = conc_face_ex(ext_bound)
        if self.bound_type == conf.TRANSFER:
            flux_out = conc_face_edge * self.boundary_transf_right
        elif self.bound_type == conf.DIFF_FLUX:
            flux_out = self.boundary_Diff * (self.boundary_diff_out - 
                        conc_face_edge) / self.distance_out
        return flux_out
               
    def run(self):
        """
        Method that is called to do a full model run
        """
        self.create_mesh()
        self.initial_yarn2d()
        if not self.cfg.onlymesh:
            self.initial_yarn2d()
            self.solve_single_component()
