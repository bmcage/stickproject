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
from fiberfipy.config import FiberfipyConfigManager
from fiberfipy.fibermodel import FiberModel

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
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = self.time_period / self.delta_t
        self.cfg_fiber = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            self.cfg_fiber.append(FiberfipyConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.cfg.get("time.time_period"))
        
        #create fiber models
        self.fiber_models = []
        for cfg in self.cfg_fiber:
            self.fiber_models.append(FiberModel(cfg))

        self.verbose = self.cfg.get('general.verbose')

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh!
        
        grid: the space position of each central point for every cell element (r-coordinate);
        """
        self.beginning_point = 0 #center of the yarn, r=0, with r the distance from center yarn.
        self.end_point = self.cfg.get('domain.yarnradius')
        self.nr_edge = 10
        self.diff_coef = 0.
        self.init_conc = lambda x: 0.0
        #we now construct the full edge grid
        self.grid_edge = sp.linspace(self.beginning_point, self.end_point, self.nr_edge)
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]
        #create cylindrical 1D grid over domain.
        self.mesh_yarn = CylindricalGrid1D(dr=tuple(self.delta_r))
        self.mesh_yarn.periodicBC = False
        self.mesh_yarn = self.mesh_yarn + (self.beginning_point,)
                     
    def initial_yarn1d(self):
        """ initial concentration over the domain"""
        self.init_conc = self.cfg.get('initial.init_conc')
        self.conc = CellVariable(name = "", 
                    mesh = self.mesh_yarn, value = self.init_conc)
        #self.viewer = None
        self.viewer = Viewer(vars = self.conc, datamin = 0., datamax = none)
#Viewer(vars = self.conc, datamin = 0., datamax =0.005)

   def solve_fiber(self):
        """
        Solve the diffusion process on the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform.
        
        """
        for model in self.fiber_models:
            model.run()
        
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
        flux_edge[1:-1] = (self.diffusion_coeff[:-1] * sp.exp(-self.diff_exp_fact * w_rep[:-1]/self.grid[:-1]) \
                         + self.diffusion_coeff[1:] * sp.exp(-self.diff_exp_fact * w_rep[1:]/self.grid[1:]))/2.\
                    * self.grid_edge[1:-1] \
                    * (w_rep[1:]/self.grid[1:] - w_rep[:-1]/self.grid[:-1])\
                    / ((self.delta_r[:-1] + self.delta_r[1:])/2.)
        diff_w_t[:]=(flux_edge[1:]-flux_edge[:-1])/self.delta_r[:]
        return diff_w_t
    
    def solve_odeint(self):
        initial_w = self.initial_c1 * self.grid
        self.solv=odeint(self.f_conc1, initial_w, self.times)
        self.conc1=self.solv/ self.grid
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
        self.view_sol(self.times, self.conc1)

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
        self.diffusion_DEET = self.cfg.get('diffusion.diffusion_conc')
        #input the trsient equation of diffusion        
        self.eq = TransientTerm() == MyDiffusionTermNoCorrection(coeff = self.diffusion_DEET)
        #get the position of the boundary faces
        xfc, yfc = self.mesh2d.getFaceCenters()
        xcc, ycc = self.mesh2d.getCellCenters()
        face_in = ((self.mesh2d.getExteriorFaces()) & 
                    (sp.power(xfc,2) + sp.power(yfc,2) \
                        < (self.grid.radius_domain - self.grid.radius_boundlayer)**2))
        face_ex = (~face_in) & (self.mesh2d.getExteriorFaces())
        self.initial_t = 0.
        filename1 = 'concentration_out.gz'
        filepath1 = utils.OUTPUTDIR + os.sep + filename1
        conc1_out_yarn = sp.zeros(1, float)
##        #calculate the bed net part
##        n_point_net = int(self.yarn_length / self.net_width) + 1
##        delta_effect = self.domain_effect / self.dis_effect
##        self.distance_yarn = sp.empty(4 * n_point_net, float)
        for i in sp.arange(0, self.steps, 1):
            ## TODO: take the blend into account!!
            conc_on_fib = self.fiber_models[0].fiber_surface[i+1]
            flux_in_fib = self.fiber_models[0].boundary_transf_right * conc_on_fib
            #loss to outside of yarn is 0.01 conc at outside
            BCs = (FixedFlux(face_ex, value = 0.01 * self.conc.getArithmeticFaceValue()), 
                   FixedFlux(face_in, value = -flux_in_fib),)
            self.eq.solve(var = self.conc, boundaryConditions = BCs, dt = self.delta_t, )
            self.initial_t += self.delta_t
            print 'time = ', (i+1) * self.delta_t
##            value_face_out = np.empty(len(face_ex), float)#save the concentration at the face-out
##            determine_out = np.empty(len(face_ex), bool)#save the boolean value at the face-out
##            for i_out in sp.arange(0, len(face_ex), 1):
##                value_face_out[i_out] = float(self.conc_face_ex[i_out])
##                determine_out[i_out] = face_ex[i_out]
##            value_out_record = value_face_out[determine_out]#get the value at the face out
##            conc1_average_out = np.sum(value_out_record) / len(value_out_record)
##            conc1_out_yarn = np.append(conc1_out_yarn, conc1_average_out)
                    
            if self.viewer is not None:
                self.viewer.plot(title = "")
##        dump.write({'time_step': self.times, 'conc_out': conc1_out_yarn},
##                                filename = filepath1, extension = '.gz')
        raw_input("Finshed <return>.....")
    
    def run(self):        
        self.create_mesh()
        self.initial_yarn2d()
        self.solve_fiber()
        self.solve_single_component()
        
