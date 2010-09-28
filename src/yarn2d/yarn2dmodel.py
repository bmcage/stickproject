#
# Copyright (C) 2010  P.Li, B. Malengier
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

""" module holding a generic diffusion model. 
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
import lib.diff.diffusion as diffusion
from mycorrection import MyDiffusionTermNoCorrection
import fibersurface

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
class Yarn2DModel():
    """
    yar2dmodel is a special diffusion model for a single yarn which is composed 
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
        self.comps = self.cfg.get('general.components')
        self.read = self.cfg.get('general.read')
        self.Ry = self.cfg.get('domain.yarnradius')
        self.Rf = self.cfg.get('fiber.radius_fiber')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        #computational radius
        self.radius_yarn = self.scaleL * self.Ry
        self.radius_fiber =  self.scaleL * self.Rf
        self.radius_boundlayer = self.radius_fiber/2.
        self.radius_domain = self.radius_yarn + self.radius_boundlayer
        self.cellsize_centre = self.cfg.get('domain.cellsize_centre')
        self.cellSize = self.cfg.get('domain.cellsize_fiber')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = self.time_period / self.delta_t
        #read the initial and boundary information for fiber
        self.n_point = self.cfg.get('fiber.n_point') #discretize the fiber radius
        self.boundary_fib_left = self.cfg.get('boundary.boundary_fib_left')
        self.boundary_fib_right = self.cfg.get('boundary.boundary_fib_right')
        self.diffusion_co_l1 =  self.cfg.get('diffusion.diffusion_co_l1')
        self.diffusion_co_l2 = self.cfg.get('diffusion.diffusion_co_l2')
        self.init_conc1_fiber = eval(self.cfg.get('initial.init_conc1_fiber'))
        self.transfer_conc1 = self.cfg.get('transfer.transfer_conc1')
        
    def create_circle_domain_gmsh(self):
        """
        create gmsh file with the circle domain of yarn and the fibers in it
        returns string with defenition, path to file
        """
        filename = 'yarn.geo'
        filepath = utils.OUTPUTDIR + os.sep + filename
        start = time.clock()
        if self.read == 'False':            
            self.type = self.cfg.get('fiber.type')
            self.circle_file = open(filepath, "w")
            self.current_point = 0
            self.x_position = sp.empty(self.number_fiber, float)
            self.y_position = sp.empty(self.number_fiber, float)
            self.x_central = 0.
            self.y_central = 0.
            self.z = 0.
            index = 1
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                    self.x_central, self.y_central, self.z, 
                                    self.cellsize_centre))
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                    self.x_central - self.radius_domain, self.y_central, 
                                    self.z, self.cellsize_centre))
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                    self.x_central, self.y_central + self.radius_domain,
                                    self.z, self.cellsize_centre))
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                    self.x_central + self.radius_domain, self.y_central,
                                    self.z, self.cellsize_centre))                      
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                    self.x_central, self.y_central - self.radius_domain,
                                    self.z, self.cellsize_centre))
            index = index + 4
            print self.cellsize_centre, self.cellSize
            for i in sp.arange(1, self.number_fiber + 1, 1):
                if i == 1:
                    #generate the position of fiber
                    a = np.random.uniform(-0.5, 0.5)
                    b = np.random.uniform(-0.5, 0.5)
                    distance_center = sp.sqrt((a - self.x_central)**2 + (b - self.y_central)**2)
                    while distance_center + self.radius_fiber >= self.radius_yarn:
                        a = np.random.uniform(-1, 1)
                        b = np.random.uniform(-1, 1)
                        distance_center = sp.sqrt((a - self.x_central)**2 + (b - self.y_central)**2)
                    else:
                        self.x_position[i-1] = a
                        self.y_position[i-1] = b
                        self.current_point = self.current_point + 1
                        index = index + 1
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                            self.x_position[i-1], self.y_position[i-1], self.z,self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                            self.x_position [i-1] - self.radius_fiber,
                                            self.y_position[i-1], self.z, self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                            self.x_position[i-1], self.y_position[i-1] + self.radius_fiber,
                                            self.z, self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                            self.x_position[i-1] + self.radius_fiber,
                                            self.y_position[i-1], self.z, self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                            self.x_position[i-1], self.y_position[i-1] - self.radius_fiber,
                                            self.z, self.cellSize))
                        index = index + 4
                elif i > 1:
                    #generate the position of fiber
                    a = np.random.uniform(-1, 1)
                    b = np.random.uniform(-1, 1)
                    #distance between the current point and center
                    distance_center = sp.sqrt((a - self.x_central)**2 + (b - self.y_central)**2)
                    #distance between the current point and existing points
                    distance_each = sp.sqrt((a - self.x_position[:self.current_point])**2 + \
                                            (b - self.y_position[:self.current_point])**2)
                    while distance_center + self.radius_fiber >= self.radius_yarn or np.min(distance_each) <= (2*self.radius_fiber): 
                        a = np.random.uniform(-1, 1)
                        b = np.random.uniform(-1, 1)
                        distance_center = sp.sqrt((a - self.x_central)**2 + (b - self.y_central)**2)
                        #distance between the current point and existing points
                        distance_each = sp.sqrt((a - self.x_position[:self.current_point])**2 + \
                                            (b - self.y_position[:self.current_point])**2)
                    else:
                        self.x_position[i-1] = a
                        self.y_position[i-1] = b
                        self.current_point = self.current_point + 1
                        index = index + 1
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                            self.x_position[i-1], self.y_position[i-1], self.z,self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                            self.x_position [i-1] - self.radius_fiber,
                                            self.y_position[i-1], self.z, self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                            self.x_position[i-1], self.y_position[i-1] + self.radius_fiber,
                                            self.z, self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                            self.x_position[i-1] + self.radius_fiber,
                                            self.y_position[i-1], self.z, self.cellSize))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                            self.x_position[i-1], self.y_position[i-1] - self.radius_fiber,
                                            self.z, self.cellSize))
                        index = index + 4
            #above part is for generating the points of circle in the domain
            print "all the points has been generated"
            index_point_in_circle = 0 #the number of each point in Circle part
            for i1 in sp.arange(0, self.number_fiber + 1, 1):
                if i1 == 0:
                    index = index + 1
                    index_point_in_circle = index_point_in_circle + 1
                    t1 = index_point_in_circle + 1
                    t2 = index_point_in_circle + 2
                    t3 = index_point_in_circle + 3
                    t4 = index_point_in_circle + 4
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index,
                                                        t1, index_point_in_circle, t2))
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index + 1,
                                                        t2, index_point_in_circle, t3))
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index + 2,
                                                        t3, index_point_in_circle, t4))
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index + 3, 
                                                        t4, index_point_in_circle, t1))
                    index = index + 3
                elif i1 > 0:
                    index = index + 1
                    index_point_in_circle = index_point_in_circle + 5
                    t1 = index_point_in_circle + 1
                    t2 = index_point_in_circle + 2
                    t3 = index_point_in_circle + 3
                    t4 = index_point_in_circle + 4
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index,
                                                        t1, index_point_in_circle, t2))
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index + 1,
                                                        t2, index_point_in_circle, t3))
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index + 2,
                                                        t3, index_point_in_circle, t4))
                    self.circle_file.write("Circle(%d) = {%d,%d,%d};\n" %(index + 3, 
                                                        t4, index_point_in_circle, t1))
                    index = index + 3
            #above part is for generating the circle part of the domain
            index_circle_for_loop = index - 4* (self.number_fiber + 1) 
            index = index + 1
            self.circle_file.write("Line Loop(%d)= {" %(index))
            for i2 in sp.arange(0, 4*(self.number_fiber + 1), 1):
                if i2 < 4*(self.number_fiber + 1) -1:
                    index_circle_for_loop = index_circle_for_loop + 1
                    self.circle_file.write("%d," %(index_circle_for_loop))
                elif i2 == 4*(self.number_fiber + 1) -1:
                    index_circle_for_loop = index_circle_for_loop +1
                    self.circle_file.write("%d" %(index_circle_for_loop))
            self.circle_file.write("};\n")
            elapsed = (time.clock() - start)
            print "How much time is consumed is: %.3f" %(elapsed)
            #above part is for generating the surface loop in the yarn domain
            print "all the circles has been generated"
            
            index_loop_in_plane = index
            index = index + 1
            self.circle_file.write("Plane Surface(%d) = {%d};\n" %(index, index_loop_in_plane))
            self.circle_file.close()
            circledef = open(filepath, "r").readlines()
        elif self.read == 'True':
            circledef = open(filepath, "r").readlines()
        return ''.join(circledef)
    
    def gmsh_2d_generate(self):
        """
        Using Gmsh2D from the package Fipy and mesh file constructed before 
        generates the 2D mesh for convection-diffusion problem in 2D. 
        """
        self.mesh2d = Gmsh2D(self.create_circle_domain_gmsh())
    
    def initial_yarn2d(self):
        self.init_conc1 = self.cfg.get('initial.init_conc1')
        self.conc1 = CellVariable(name = "solution concentration1", 
                    mesh = self.mesh2d, value = self.init_conc1)
        self.viewer = None
        self.viewer = Viewer(vars = self.conc1, datamin = 0., datamax =0.00002)

    def solve_single_component_fiber(self):
        """
        Solve the diffusion process on the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        self.beginning_point = self.cfg.get('fiber.beginning_point')
        self.end_point = self.cfg.get('fiber.end_point')
        scale_beginning = self.beginning_point * self.scaleL
        print 'this is the scale beginning point:',scale_beginning
        scale_end = self.end_point * self.scaleL
        print 'this is the scale end point', scale_end
        self.grid = sp.linspace(scale_beginning, scale_end, self.n_point)
        #self.grid = grid
        initial_c1 = sp.empty(self.n_point, float)
        discretization_t = self.steps + 1
        times = sp.linspace(0, self.time_period, discretization_t)
        for i in sp.arange(0, self.n_point, 1):
            if i <= (self.n_point - 1) / 2:
                initial_c1[i] = self.init_conc1_fiber(i)[0]
            elif i > (self.n_point - 1) /2:
                initial_c1[i] = self.init_conc1_fiber(i)[1]
        print 'initial condition is:', initial_c1
        self.fiber_conc1 = fibersurface.Solving1DFiber(self.grid, initial_c1, 
                           self.boundary_fib_left, self.boundary_fib_right,
                           self.diffusion_co_l1, self.diffusion_co_l2)
        self.fiber_conc1.solver_c(times)
        self.fiber_surface = sp.empty(len(times), float)
        for i in sp.arange(1,len(times) + 1,1):
            self.fiber_surface[i - 1] = self.fiber_conc1.conc1[i - 1][-1]
        print self.fiber_surface[:]

    def solve_single_component(self):
        """
        The DEET diffusion process is divided into two parts:
        (1) DEET diffuses through the layer containing permithrine on the fiber 
        and reache the surface;
        (2) DEET begins to diffuse in the void space of  yarn
        So it means that the boundary condition of fiber has two steps:
        (1) When the DEET does not reach the surface of fiber, the inner and out
        boundaries are no flux boundary condition;
        (2) When the DEET reaches surface, the evaporation happens. So the boundaries 
        of fiber and yarn are changed to constant flux (Neumann boundary condition)
        """
        self.diffusion_DEET = self.cfg.get('diffusion.diffusion_conc1')
        #input the trsient equation of diffusion        
        self.eq = TransientTerm() == MyDiffusionTermNoCorrection(coeff = self.diffusion_DEET)
        #get the position of the boundary faces
        xfc, yfc = self.mesh2d.getFaceCenters()
        xcc, ycc = self.mesh2d.getCellCenters()
        face_in = ((self.mesh2d.getExteriorFaces()) & 
                    (sp.power(xfc,2) + sp.power(yfc,2) \
                        < (self.radius_domain - self.radius_boundlayer)**2))
        face_ex = (~face_in) & (self.mesh2d.getExteriorFaces())
        #BCs = (FixedFlux(face_ex, value = 0.), FixedValue(face_in, value = 1.),)
        """
        for i in sp.arange(0, self.steps, 1):
            surface_value = self.fiber_surface[i]
            BCs = (FixedFlux(face_ex, value = 0.), FixedValue(face_in, value = surface_value),)
            self.eq.solve(var = self.conc1, boundaryConditions = BCs, dt = self.delta_t, )
            print 'time = ', (i+1) * self.delta_t
            #raw_input("Finshed <return>.....")
            if self.viewer is not None:
                self.viewer.plot()
                #raw_input("continue to next step, please enter <return>.....")
        """
        i1 = 0
        self.initial_t = 0.
        initial_c2 = sp.empty(len(self.grid),float)
        filename = 'fiber_layer.data'
        filepath = utils.OUTPUTDIR + os.sep + filename
        self.fiber_file = open(filepath, "w")
        for i in sp.arange(0, self.steps, 1):
            determine_value = self.fiber_surface[i]
            if determine_value <= 0:
                surface_value = determine_value
                BCs = (FixedFlux(face_ex, value = sp.zeros(len(face_ex), float)), FixedValue(face_in, value = 0),)
                self.eq.solve(var = self.conc1, boundaryConditions = BCs, dt = self.delta_t, )
                print 'time = ', (i+1) * self.delta_t
                #raw_input("Finshed <return>.....")
            elif determine_value > 0:#the evaporation happens on the surface of the layer
                i1 = i1 + 1
                if i1 == 1:
                    surface_value = determine_value
                    self.boundary_fib_right = self.transfer_conc1 * surface_value
                    boundary_in_yarn = self.boundary_fib_right
                    BCs = (FixedFlux(face_ex, value = 0.), FixedFlux(face_in, value = boundary_in_yarn),)
                    self.eq.solve(var = self.conc1, boundaryConditions = BCs, dt = self.delta_t, )
                    try:
                        self.conc_face_ex = self.conc1.getArithmeticFaceValue()#[self.mesh2d.getExteriorFaces().getValue()]
                        print self.conc_face_ex
                    except ValueError:
                        print 'the method has problem'
                    for i2 in sp.arange(0, len(self.grid), 1):
                        initial_c2[i2] = self.fiber_conc1.conc1[i][i2]
                        #self.fiber_file.write("%f, %f \n" %(self.grid[i2], initial_c2[i2])) 
                    #self.fiber_file.close() 
                    self.fiber_conc1 = fibersurface.Solving1DFiber(self.grid, initial_c2, 
                                   self.boundary_fib_left, -self.boundary_fib_right,
                                   self.diffusion_co_l1, self.diffusion_co_l2)
                    self.fiber_conc1.solver_c_2(self.delta_t, self.initial_t)
                    self.domain_conc1 = self.fiber_conc1.conc1
                    surface_value = self.domain_conc1[-1]
                    self.initial_t = self.initial_t + self.delta_t
                    print 'time = ', (i+1) * self.delta_t
                    #raw_input("Finshed <return>.....")
                elif i1 > 1:
                    self.boundary_fib_right = self.transfer_conc1 * surface_value
                    boundary_in_yarn = self.boundary_fib_right
                    print 'the inner boundary of yarn', boundary_in_yarn
                    #boundary_yarn_out_estep = sp.empty(len(face_ex), float)
                    conc_face_ex = sp.zeros(len(face_ex), float)
                    conc_face_ex = 0.001 * self.conc_face_ex
                    """
                    for i3 in sp.arange(0, len(face_ex), 1):
                        if face_ex[i3] == True:
                            conc_face_ex[i3] = 0.001 * self.conc_face_ex[i3]
                    """
                    BCs = (FixedFlux(face_ex, value = conc_face_ex), FixedFlux(face_in, value = -boundary_in_yarn),)
                    self.eq.solve(var = self.conc1, boundaryConditions = BCs, dt = self.delta_t, ) 
                    initial_c2 = self.domain_conc1
                    self.conc_face_ex = self.conc1.getArithmeticFaceValue()
                    print 'initial in fiber', initial_c2
                    if (i+1) * self.delta_t == 10:
                        for i2 in sp.arange(0, len(self.grid), 1):
                           self.fiber_file.write("%f, %f \n" %(self.grid[i2], initial_c2[i2]))
                        self.fiber_file.close()
                    self.fiber_conc1 = fibersurface.Solving1DFiber(self.grid, initial_c2, 
                                   self.boundary_fib_left, -self.boundary_fib_right,
                                   self.diffusion_co_l1, self.diffusion_co_l2)
                    self.fiber_conc1.solver_c_2(self.delta_t, self.initial_t)
                    self.domain_conc1 = self.fiber_conc1.conc1
                    surface_value = self.domain_conc1[-1]            
                    self.initial_t = self.initial_t + self.delta_t
                    print 'time = ', (i+1) * self.delta_t
                    #raw_input("Finshed <return>.....")
            if self.viewer is not None:
                self.viewer.plot()
        raw_input("Finshed <return>.....")
        
    def run(self):        
        self.gmsh_2d_generate()
        self.initial_yarn2d()
        self.solve_single_component_fiber()
        self.solve_single_component()
        
