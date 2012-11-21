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
import shutil
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
import stick.lib.utils.utils as utils

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
        length = self.cfg.get('domain.length')
        width = self.cfg.get('domain.width')
        height = self.cfg.get('domain.height')
        el_length = self.cfg.get('discretization.el_length')
        el_width = self.cfg.get('discretization.el_width')
        el_height = self.cfg.get('discretization.el_height')
        flength = self.cfg.get('fabric.length')
        fwidth = self.cfg.get('fabric.width')
        fheight = self.cfg.get('fabric.height')
        
        dx = length / el_length # mesh size in x direction 
        dy = width / el_width # mesh size in x direction 
        dz = height / el_height # mesh size in x direction 
            
        load_msh = self.cfg.get('domain.load_msh')
        if load_msh:
            filenamemsh = os.path.normpath(os.path.join(
                                os.path.dirname(self.cfg.filename), 
                                self.cfg.get('domain.msh_file')))
            if not os.path.isfile(filenamemsh):
                raise Exception('File mesh file does not exist: %s' % filenamemsh)
        else:
            self.meshsize = height/5
            self.meshsizesmall = self.meshsize
            meshgeo = """cl1 = %(ms)g;
cl2 = %(ms_small)g;
Point(1) = {%(L)g, %(W)g, 0, cl1};
Point(2) = {%(L)g, -%(W)g, 0, cl1};
Point(3) = {-%(L)g, -%(W)g, 0, cl1};
Point(4) = {-%(L)g, %(W)g, 0, cl1};
Point(5) = {-%(L)g, %(W)g, %(H)g, cl1};
Point(6) = {-%(L)g, -%(W)g, %(H)g, cl1};
Point(7) = {%(L)g, -%(W)g, %(H)g, cl1};
Point(8) = {%(L)g, %(W)g, %(H)g, cl1};
Point(9) = {%(l)g, %(w)g, 0, cl1};
Point(10) = {%(l)g, -%(w)g, 0, cl1};
Point(11) = {-%(l)g, -%(w)g, 0, cl1};
Point(12) = {-%(l)g, %(w)g, 0, cl1};
Point(13) = {%(l)g, %(w)g, %(h)g, cl2};
Point(14) = {%(l)g, -%(w)g, %(h)g, cl2};
Point(15) = {-%(l)g, -%(w)g, %(h)g, cl2};
Point(16) = {-%(l)g, %(w)g, %(h)g, cl2};
Point(17) = {%(l)g, %(w)g, %(overlaph)g, cl1};
Point(18) = {-%(l)g, %(w)g, %(overlaph)g, cl1};
Point(19) = {-%(l)g, -%(w)g, %(overlaph)g, cl1};
Point(20) = {%(l)g, -%(w)g, %(overlaph)g, cl1};
Line(1) = {3, 6};
Line(2) = {6, 5};
Line(3) = {5, 8};
Line(4) = {8, 7};
Line(5) = {7, 6};
Line(6) = {2, 7};
Line(7) = {1, 8};
Line(8) = {4, 5};
Line(9) = {3, 2};
Line(10) = {2, 1};
Line(11) = {1, 4};
Line(12) = {4, 3};
Line(13) = {11, 15};
Line(14) = {15, 19};
Line(15) = {12, 12};
Line(16) = {12, 16};
Line(17) = {16, 18};
Line(18) = {9, 13};
Line(19) = {13, 17};
Line(20) = {10, 14};
Line(21) = {14, 20};
Line(22) = {10, 11};
Line(23) = {11, 12};
Line(24) = {12, 9};
Line(25) = {9, 10};
Line(26) = {14, 15};
Line(27) = {15, 16};
Line(28) = {16, 13};
Line(29) = {13, 14};
Line(30) = {20, 19};
Line(31) = {19, 18};
Line(32) = {18, 17};
Line(33) = {17, 20};
Line(54) = {4, 12};
Line(55) = {11, 3};
Line(56) = {10, 2};
Line(57) = {9, 1};
Line Loop(35) = {30, -14, -26, 21};
Plane Surface(35) = {35};
Line Loop(37) = {14, 31, -17, -27};
Plane Surface(37) = {37};
Line Loop(39) = {17, 32, -19, -28};
Plane Surface(39) = {39};
Line Loop(41) = {19, 33, -21, -29};
Plane Surface(41) = {41};
Line Loop(43) = {30, 31, 32, 33};
Plane Surface(43) = {43};
Line Loop(45) = {28, 29, 26, 27};
Plane Surface(45) = {45};
Line Loop(47) = {27, -16, -23, 13};
Plane Surface(47) = {47};
Line Loop(49) = {24, 18, -28, -16};
Plane Surface(49) = {49};
Line Loop(51) = {29, -20, -25, 18};
Plane Surface(51) = {51};
Line Loop(53) = {26, -13, -22, 20};
Plane Surface(53) = {53};
Line Loop(59) = {54, 24, 57, 11};
Plane Surface(59) = {59};
Line Loop(61) = {57, -10, -56, -25};
Plane Surface(61) = {61};
Line Loop(63) = {56, -9, -55, -22};
Plane Surface(63) = {63};
Line Loop(65) = {55, -12, 54, -23};
Plane Surface(65) = {65};
Line Loop(67) = {11, 8, 3, -7};
Plane Surface(67) = {67};
Line Loop(69) = {4, -6, 10, 7};
Plane Surface(69) = {69};
Line Loop(71) = {2, -8, 12, 1};
Plane Surface(71) = {71};
Line Loop(73) = {5, -1, 9, 6};
Plane Surface(73) = {73};
Line Loop(75) = {5, 2, 3, 4};
Plane Surface(75) = {75};
Surface Loop(77) = {43, 35, 37, 39, 41, 45};
Volume(77) = {77};
Surface Loop(79) = {75, 73, 71, 67, 59, 65, 63, 61, 69, 51, 53, 47, 49, 43, 35, 37, 39, 41};
Volume(79) = {79};
""" % { 
        'ms':   self.meshsize,
        'ms_small': self.meshsizesmall,
        'L': length/2,
        'W': width/2,
        'H': height,
        'l': flength/2,
        'w': fwidth/2,
        'h': fheight,
        'overlaph': fheight + fheight,
        }
            filenamegeo = utils.OUTPUTDIR + os.sep + 'room.geo'
            filenamemsh = utils.OUTPUTDIR + os.sep + 'room.msh'
            filegeo = open(filenamegeo, 'wb')
            filegeo.write(meshgeo)
            filegeo.close()
            #new refine the file a couple of times
            refine = 0
            from subprocess import Popen, check_call
            check_call(['gmsh', '-3', filenamegeo, '-o', filenamemsh ])
            for ind in range(refine):
                print 'refining', ind, 'time'
                check_call(['gmsh', '-refine', filenamemsh, '-o', filenamemsh])

        # construct the fipy mesh with the filenamemsh
        self.mesh = Gmsh3D(filenamemsh)
##            self.mesh = Grid3D(dx=dx, nx=el_length, dy=dy, ny=el_width, dz=dz, 
##                   nz=el_height)

        #get the position of the boundary faces
        xfc, yfc, zfc = self.mesh.faceCenters
        # define different sizes for the boundary conditions
        self.facesLeft = (xfc < -length/2 + 1e-6)
        self.facesRight = (xfc > length/2 - 1e-6)
        self.facesTop = (zfc > height - 1e-6)
        self.facesBottom = (zfc < 1e-6)
        self.facesFront = (yfc < -width/2 + 1e-6)
        self.facesBack = (yfc > width/2 - 1e-6)

        self.facesBound = (self.facesLeft | self.facesRight | 
                        self.facesTop | self.facesBottom |
                        self.facesBack | self.facesFront)
        self.facesTextile = (xfc > -flength/2 - 1e-6) & (xfc < flength/2 + 1e-6) \
                & (yfc > -fwidth/2 - 1e-6) & (yfc < fwidth/2 + 1e-6) \
                & (zfc < fheight + 1e-6) & (zfc > fheight - 1e-6)
##        print 'test Left', list(xfc[self.facesLeft])
##        print 'test Right', list(xfc[self.facesRight])
##        print 'test Top', list(zfc[self.facesTop])
##        print 'test Front', list(yfc[self.facesFront])
##        print 'test Back', list(yfc[self.facesBack])
##        print 'test face Text', list(zfc[self.facesTextile])
##        print list(xfc[self.facesTextile])
##        print list(yfc[self.facesTextile])
##        #print 'bound', xfc[self.facesBound]
##        sys.exit()

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
                            datamin = 15., datamax = 45.)
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

        #Dirichlet boundary conditions
        self.concTmp.constrain(self.valueDirTmp, self.facesBound)
        self.concTmp.constrain(40., self.facesTextile)
        
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
                    outvtk =  utils.OUTPUTDIR + os.sep + \
                                        'roomTemp%08.3f.vtk' % t
                    invtk = self.viewer.vtkcellfname
                    shutil.copy(invtk, outvtk)
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
