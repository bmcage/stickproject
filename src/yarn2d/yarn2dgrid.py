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
Module for functions of a yarn 2D grid. 
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
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
from fipy import Gmsh2D


#-------------------------------------------------------------------------
#
# Yarn2dGrid class
#
#-------------------------------------------------------------------------

class Yarn2dGrid(object):
    def __init__(self, cfg):
        self.cfg = cfg
        
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
        
        self.verbose = self.cfg.get('general.verbose')

    def create_circle_domain_gmsh(self, filename='yarn.geo', regenerate=True):
        """
        Create gmsh file with the circle domain of yarn and the fibers in it
        returns string with defenition, path to file
        """
        filepath = utils.OUTPUTDIR + os.sep + filename
        start = time.clock()
        if regenerate:            
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
                                            self.x_position[i-1], self.y_position[i-1], 
                                            self.z,self.cellSize))
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
                    while distance_center + self.radius_fiber >= self.radius_yarn or np.min(distance_each)\
                            <= (2*self.radius_fiber): 
                        a = np.random.uniform(-1, 1)
                        b = np.random.uniform(-1, 1)
                        distance_center = sp.sqrt((a - self.x_central)**2 + \
                        (b - self.y_central)**2)
                        #distance between the current point and existing points
                        distance_each = sp.sqrt((a - self.x_position[:self.current_point])**2 + \
                                            (b - self.y_position[:self.current_point])**2)
                    else:
                        self.x_position[i-1] = a
                        self.y_position[i-1] = b
                        self.current_point = self.current_point + 1
                        index = index + 1
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                            self.x_position[i-1], self.y_position[i-1], 
                                            self.z,self.cellSize))
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
            if self.verbose:
                print "all the points have been generated"
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
            if self.verbose:
                print "Time to generate the mesh file: %.3f" %(elapsed)
                #above part is for generating the surface loop in the yarn domain
                print "Mesh generation finished"
            
            index_loop_in_plane = index
            index = index + 1
            self.circle_file.write("Plane Surface(%d) = {%d};\n" %(index, index_loop_in_plane))
            self.circle_file.close()
            circledef = open(filepath, "r").readlines()
        else:
            circledef = open(filepath, "r").readlines()
        return ''.join(circledef)

    def mesh_2d_generate(self, filename='yarn.geo', regenerate=True):
        """
        Return a Gmsh2D object from Fipy for the yarn 2D grid
        The gmsh file is written to filename.
        If regenerate is True, it is looked if the file already exists, and if
        so, the file is reused instead of generated.
        """
        self.mesh = Gmsh2D(self.create_circle_domain_gmsh(filename, regenerate))
        return self.mesh
