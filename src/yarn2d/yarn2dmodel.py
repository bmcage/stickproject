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
import yarn2d.config as conf
from mycorrection import MyDiffusionTermNoCorrection
from yarn2dgrid import Yarn2dGrid
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
        self.datatime = []
        self.cfg = config
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.step = self.cfg.get('time.step')
        self.steps = (self.time_period*(1.+self.step*1e-6)) // self.step #self.delta_t
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.delta_t = self.times[1] - self.times[0]
        self.Ry = self.cfg.get('domain.yarnradius')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        self.Rf = self.cfg.get('fiber.radius_fiber')
        self.radius_fiber =  [self.scaleL * rad for rad in self.Rf]
        self.eps_value = self.cfg.get('fiber.eps_value')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        self.nrtypefiber = self.cfg.get('fiber.number_type')
        self.fiber_edge_result = [0] * self.nrtypefiber
        assert self.nrtypefiber == len(self.blend) == len(self.Rf)
        #computational radius
        self.cfg_fiber = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            self.cfg_fiber.append(Fiber1dConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.cfg.get("time.time_period"))
            if self.cfg_fiber[-1].get("time.dt") > self.cfg.get("time.time_period"):
                self.cfg_fiber[-1].set("time.dt", self.cfg.get("time.time_period"))
        
        #create fiber models
        self.fiber_models = []
        for cfg in self.cfg_fiber:
            self.fiber_models.append(FiberModel(cfg))

        self.verbose = self.cfg.get('general.verbose')
    
    def create_mesh(self):
        """
        Create a mesh for use in the model
        """
        self.grid = Yarn2dGrid(self.cfg)
        self.mesh2d = self.grid.mesh_2d_generate(filename='yarn.geo',
                                regenerate=not self.cfg.get('general.read'))
        if self.cfg.onlymesh:
            raw_input("Finished Mesh generation, press Enter to continue")
    
    def initial_yarn2d(self):
        self.init_conc = self.cfg.get('initial.init_conc')
        self.conc = CellVariable(name = "Conc. Active Component", 
                    mesh = self.mesh2d, value = self.init_conc)
        ##self.conc1 = CellVariable(name = "solution concentration1",
        ##            mesh = self.mesh2d, value = self.init_conc)
        ##self.conc2 = CellVariable(name = "solution concentration2",
        ##            mesh = self.mesh2d, value = self.init_conc)
        self.viewer = None
        self.viewer = Viewer(vars = self.conc, datamin = 0., datamax =0.0005)
    def solve_fiber_init(self):
        for model in self.fiber_models:
            model.run_init()

    def solve_fiber_step(self, step):
        """
        Solve the diffusion process on the fiber for a single step, starting
        from where we where last. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        for nyfib, model in enumerate(self.fiber_models):
            #model.run_init()
            print 'test', model.tot_edges
            self.fiber_edge_result[nyfib] = model.run_step(step)[-1]
            #times of the models must coincide at the moment as later on we do
            #conc_on_fib[nyfib] = self.fiber_models[nyfib].fiber_surface[i+1]
            #we should interpolate to avoid that
            #assert len(self.times) == len(model.times), 'fiber and yarn model need same time steps'
            #test = (sp.array(self.times) == sp.array(model.times))
            #assert test.all(), 'fiber and yarn model need same time steps'
            
    def solve_fiber(self):
        """
        Solve the diffusion process on the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        for model in self.fiber_models:
            model.run()
            #self.fiber_edge_result[nyfib] = model.run
            #times of the models must coincide at the moment as later on we do
            #conc_on_fib[nyfib] = self.fiber_models[nyfib].fiber_surface[i+1]
            #we should interpolate to avoid that
            assert len(self.times) == len(model.times), 'fiber and yarn model need same time steps'
            test = (sp.array(self.times) == sp.array(model.times))
            assert test.all(), 'fiber and yarn model need same time steps'

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
        #self.solve_fiber()
        self.solve_fiber_init()
        #self.solve_fiber_step()
        self.diffusion_DEET = self.cfg.get('diffusion.diffusion_conc')
        #input the trsient equation of diffusion        
        self.eq = TransientTerm() == MyDiffusionTermNoCorrection(coeff = self.diffusion_DEET)
        #get the position of the boundary faces
        
        conc_fib_out1 = sp.empty(self.steps)
        conc_fib_out2 = sp.empty(self.steps)
        
        xfc, yfc = self.mesh2d.getFaceCenters()
        xcc, ycc = self.mesh2d.getCellCenters()
        self.cell_volume = self.mesh2d.getCellVolumes()
        face_in = ((self.mesh2d.getExteriorFaces()) & 
                    (sp.power(xfc,2) + sp.power(yfc,2) \
                        < (self.grid.radius_domain - self.grid.radius_boundlayer)**2))
        filepath3 = utils.OUTPUTDIR + os.sep + 'index_fiber.dat'
        filepath4 = utils.OUTPUTDIR + os.sep + 'yarn_out.dat'
        filepath5 = utils.OUTPUTDIR + os.sep + 'conc_fib1.gz'
        filepath6 = utils.OUTPUTDIR + os.sep + 'conc_fib2.gz'
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
        #outside
        self.ext_bound = ((self.mesh2d.getExteriorFaces()) & 
                    (sp.power(xfc,2) + sp.power(yfc,2) \
                        < (self.grid.radius_domain - self.grid.radius_boundlayer)**2))
        
        #we need to determine which nodes in the mesh are surface of 
        #of a certain fiber kind and create the inner boundary
        eps_fib = [val * self.eps_value for val in self.radius_fiber]
        self.int_bound = []
        for nyfib in sp.arange(self.nrtypefiber):
            tmp = np.zeros(len(face_in), bool)
            for fib in self.index_fiber[nyfib]:
                tmp = (((self.mesh2d.getExteriorFaces()) &
                    (sp.power(xfc - self.fib_x[fib], 2) + sp.power(yfc - self.fib_y[fib], 2)\
                    <= (self.radius_fiber[nyfib] + eps_fib[nyfib])**2)) | (tmp))
            self.int_bound.append(tmp)
        face_ex = (~face_in) & (self.mesh2d.getExteriorFaces())
        self.initial_t = 0.
        filename1 = 'concentration_out.gz'
        filepath1 = utils.OUTPUTDIR + os.sep + filename1
        conc1_out_yarn = sp.zeros(self.steps, float)
        conc_on_fib = sp.empty(self.nrtypefiber)
        flux_in_fib = sp.empty(self.nrtypefiber)

        ## TODO: IMPROVE THIS BC
        conc1_out_yarn = []
        value_face_out = np.empty(len(face_ex), float)
        determine_out = np.empty(len(face_ex), bool)
        self.record_conc = open(filepath4, "w")
        for i in sp.arange(0, self.steps, 1):
            self.solve_fiber_step(self.delta_t)
            BCs = []
            BCs.append(FixedFlux(face_ex, value = 0.0))
            for nyfib in sp.arange(self.nrtypefiber):
                ## TODO, fix so that self.times need not be == model.times
##                conc_on_fib[nyfib] = (self.fiber_models[nyfib].fiber_surface[i] +
##                                    self.fiber_models[nyfib].fiber_surface[i+1]) / 2.
                conc_on_fib[nyfib] = self.fiber_edge_result[nyfib]
                
                flux_in_fib[nyfib] = self.fiber_models[nyfib].boundary_transf_right * conc_on_fib[nyfib]
                BCs.append(FixedFlux(self.int_bound[nyfib], value = -flux_in_fib[nyfib]))#[nyfib + 1] = FixedFlux(self.int_bound[nyfib], value = -flux_in_fib[nyfib])
            conc_fib_out1[i] = conc_on_fib[0]
            conc_fib_out2[i] = conc_on_fib[-1]
            self.initial_t += self.delta_t
            self.eq.solve(var = self.conc, boundaryConditions = tuple(BCs), dt = self.delta_t, )
            print 'time = ', (i+1) * self.delta_t
            self.conc_tot_each = self.conc.getValue()
            self.conc_face_ex = self.conc.getArithmeticFaceValue()
            for i_out in sp.arange(len(face_ex)):
                value_face_out[i_out] = float(self.conc_face_ex[i_out])
                determine_out[i_out] = face_ex[i_out]
            value_out_record = value_face_out[determine_out]
            conc1_average_out = np.sum(value_out_record) / len(value_out_record)
            print 'average conccentration out', conc1_average_out
            conc1_out_yarn = np.append(conc1_out_yarn, conc1_average_out)
            self.record_conc.write("%g, %g \n" %(self.times[i], conc1_out_yarn[-1]))
            print 'mass conservative with two fibers', self.cal_mass_void(self.conc_tot_each,
                                                self.cell_volume) / self.scaleL
##            if i%10 == 0:                                    
            if self.viewer is not None:
                self.viewer.plot()
##        dump.write({'time_step': self.times, 'conc_out': conc1_out_yarn},
##                               filename = filepath4, extension = '.gz')
        self.record_conc.close()
           # raw_input("next time step <return>....")
        dump.write({'time_step': self.times, 'conc_fib1': conc_fib_out1}, 
                    filename = filepath5, extension = '.gz')
        dump.write({'time_step': self.times, 'conc_fib2': conc_fib_out2},
                    filename= filepath6, extension = '.gz')
        raw_input("Finished <press return>.....")
    
    def cal_mass_void(self, conc_void, cell_volume):
        """
        calculate the mass of component in the void space 
        conc_void: the concentration of materials in the void space
        This is given by: total_mass = sum(conc_void[] * cell_volume[])
        """
        return sp.sum(conc_void * cell_volume)
        
    
    def run(self):
        self.create_mesh()
        if not self.cfg.onlymesh:
            self.initial_yarn2d()
            self.solve_single_component()
