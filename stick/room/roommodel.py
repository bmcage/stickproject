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
Module holding a generic diffusion model for a room used for textile exp. 
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
import scipy as sp
import matplotlib.pyplot as plt
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import stick.const as const
import stick.room.config as conf

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
class RoomModel(object):
    """
    RoomModel is a special diffusion model for a room in which we will do 
    textile experiments.
    Fipy solves the transient diffusion problem in the whole domain
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
            print "Timestep used in room model:", self.delta_t
        
        #construct cfg for fabric
        self.fabpos = self.cfg.get('fabric.fabricposition')
        if conf.PLACEMENT[self.fabpos] == conf.NONE:
            #no fabric, only model the room
            self.cfg_fabric = None
            self.fabric_model = None
        else:
            filename = self.cfg.get('fabric.fabric_config')
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            from stick.fiberfabric.config import FiberFabricConfigManager
            from stick.fiberfabric.fiberfabricmodel import FiberFabricModel
            self.cfg_fabric = FiberFabricConfigManager.get_instance(filename)
            #set values from the yarn on this inifile
            self.cfg_fabric.set("time.time_period", self.time_period)
            if self.cfg_fabric.get("time.dt") > self.cfg.get("time.dt"):
                self.cfg_fabric.set("time.dt", self.cfg.get("time.dt"))
            if conf.PLACEMENT[self.fabpos] == conf.BOTCENT:
                raise NotImplementedError
            else:
                raise NotImplementedError('Wrong placement given')
                
            #create fabric model
            self.fabric_model = FiberFabricModel(cfg)

        self.plotevery = self.cfg.get("plot.plotevery")
        self.writeevery = self.cfg.get("plot.writeevery")
        self.writeoutcount = 0
        
        self.initialized = False
        
    def create_mesh(self):
        """
        Create a mesh for use in the model
        """
        from fipy.meshes import Grid3D
        length = self.cfg.get('domain.length')
        width = self.cfg.get('domain.width')
        height = self.cfg.get('domain.height')
        el_length = self.cfg.get('discretization.el_length')
        el_width = self.cfg.get('discretization.el_width')
        el_height = self.cfg.get('discretization.el_height')
        
        dx = length / el_length # mesh size in x direction 
        dy = width / el_width # mesh size in x direction 
        dz = height / el_height # mesh size in x direction 

        if self.fabric_model:
            raise NotImplementedModel
        else:
            self.mesh = Grid3D(dx=dx, nx=el_length, dy=dy, ny=el_width, dz=dz, 
                   nz=el_height)

    def initial_room(self):
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        
        self.init_concvap = eval(self.cfg.get('initial.init_concvap'))
        self.init_concair = eval(self.cfg.get('initial.init_concair'))
        self.init_temp = eval(self.cfg.get('initial.init_temp'))
        
        cellCenter_x, cellCenter_y, cellCenter_z = self.mesh.getCellCenters()
        initialConcVap = []
        initialConcAir = []
        initialTemp = []
        for i_x, i_y, i_z in zip(cellCenter_x, cellCenter_y, cellCenter_z):
            initialConcVap.append(self.init_concvap(i_x, i_y, i_z))
            initialConcAir.append(self.init_concair(i_x, i_y, i_z))
            initialTemp.append(self.init_temp(i_x, i_y, i_z))
        initialConcVap = sp.array(initialConcVap)
        initialConcAir = sp.array(initialConcAir)
        initialTemp = sp.array(initialTemp)
        
        self.step_old_sol_vap = initialConcVap[0]
        self.step_old_sol_air = initialConcAir[0]
        self.step_old_sol_tmp = initialTemp[0]
        
        if self.fabric_model:
            self.solve_fabric_init()
        self.conVap = CellVariable(name = "Conc. water vapour", 
                                    mesh = self.mesh,
                                    value = initialConcVap)
        self.concAir = CellVariable(name = "Conc. air", 
                                    mesh = self.mesh,
                                    value = initialConcAir)
        self.concTmp = CellVariable(name = "Temperature", 
                                    mesh = self.mesh,
                                    value = initialTemp)

        # constans in the equations
        self.Da = self.cfg.get('roomcoeff.diff_coef')
        self.Ka = self.cfg.get('roomcoeff.therm_cond_K')
        self.ca = self.cfg.get('roomcoeff.spec_heat_c')

        self.valueDirTmp = self.cfg.get('boundary.dirichletval_T_BC')

        self.viewer = None
        self.viewer = Viewer(vars = self.concTmp, title = 'Temperature Distribution', 
                            datamin = 0., datamax = 30.)
        self.viewer.plot()
        #raw_input("take the example of the initial condition")
        self.viewerplotcount = 1
        self.viewerplotcount = self.viewerplotcount % self.plotevery

    def solve_fabric_init(self):
        """
        Initialize the solver that does the fabric simulation
        """
        self.fabric_model.run_init()
        self.fabric_model.solve_init()

    def solve_room(self):
        """
        Solve the unknowns in the room
        """
        #input the transient equation
        self.eqVap = TransientTerm() == DiffusionTerm(coeff=self.Da)
        self.eqAir = TransientTerm() == DiffusionTerm(coeff=self.Da)
        self.eqTmp = TransientTerm(coeff=self.ca) == DiffusionTerm(coeff=self.Ka)
                            
        #get the position of the boundary faces
        xfc, yfc, zfc = self.mesh.faceCenters

        #Dirichlet boundary conditions
        facesBound = (self.mesh.facesLeft | self.mesh.facesRight | 
                        self.mesh.facesTop | self.mesh.facesBottom |
                        self.mesh.facesBack | self.mesh.facesFront)
        
        self.concTmp.constrain(self.valueDirTmp, facesBound)
        
        #all other boundaries are automatically Neumann BC
        
        #now loop in time to solve
        t = self.initial_t
        stop_time = self.times[-1]
        compute = True
        while compute:
            t += self.delta_t
            print 'computing time', t
            if t >= stop_time -self.delta_t / 100:
                t = stop_time
                compute = False
            self.eqTmp.solve(var = self.concTmp,
                          dt = self.delta_t)
            
            if self.viewer is not None and self.viewerplotcount == 0:
                    self.viewer.plot()
            self.viewerplotcount += 1
            self.viewerplotcount = self.viewerplotcount % self.plotevery

            if self.writeoutcount == 0:
                pass
            self.writeoutcount += 1
            self.writeoutcount = self.writeoutcount % self.writeevery

        raw_input("Finished <press return>.....")

    def run(self):
        """
        Method that is called to do a full model run
        """
        self.create_mesh()
        self.initial_room()
        self.solve_room()
