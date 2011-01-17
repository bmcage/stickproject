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
# Yarn2dGrid new class
#
#-------------------------------------------------------------------------

class Yarn2dNewGrid(object):
    def __init__(self, cfg):
        self.cfg = cfg
        
        self.Ry = self.cfg.get('domain.yarnradius')
        self.Rf = self.cfg.get('fiber.radius_fiber')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        self.radius_yarn = self.scaleL * self.Ry
        self.radius_fiber =  [self.scaleL * rad for rad in self.Rf]
        self.radius_boundlayer = max(self.radius_fiber)/2.
        self.radius_domain = self.radius_yarn + self.radius_boundlayer
        self.cellsize_centre = self.cfg.get('domain.cellsize_centre')
        self.cellSize = self.cfg.get('domain.cellsize_fiber')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        
        self.verbose = self.cfg.get('general.verbose')
        
    def create_circle_domain_gmesh(self, filename='yarn_new.geo', filename1 = 'fib_centers_x_new.geo',
        filename2 = 'fib_centers_y_new.geo' ,regenerate=True):
        filepath = utils.OUTPUTDIR + os.sep + filename
        filepath1 = utils.OUTPUTDIR + os.sep + filename1
        filepath2 = utils.OUTPUTDIR + os.sep + filename2
        self.number_for_circles = int((self.radius_yarn - self.radius_fiber[0]) / (2 * \
                                    self.radius_fiber[0]))
        self.total_circles = self.number_for_circles + 1
        self.radius_circle = sp.empty(self.total_circles, float)
        self.theta_value = self.cfg.get('domain.theta_value')
        self.beta_value = self.cfg.get('domain.beta_value')
        self.x_position = sp.empty(self.number_fiber, float)
        self.y_position = sp.empty(self.number_fiber, float)
        for i_circle in sp.arange(self.total_circles):
            if i_circle == 0:
                self.radius_circle[i_circle] = (i_circle + 1) * self.radius_fiber[0]
            elif i_circle > 0:
                self.radius_circle[i_circle] = i_circle * 2. * self.radius_fiber[0]
        """
        Generate virtual location part
        """
        #calculate the number of virtual location
        self.number_circle = sp.empty(self.total_circles, int)
        for i_circle in sp.arange(self.total_circles):
            if i_circle == 0:
                self.number_circle[i_circle] = 1
            else:
                self.number_circle[i_circle] = int(sp.pi / (self.radius_fiber[0] \
                                                / self.radius_circle[i_circle]))
                print 'the value', self.radius_fiber[0] /self.radius_circle[i_circle]
                print 'the number in each circle', self.number_circle[i_circle]
        total_number_vl = sum(self.number_circle[:])
        if self.number_fiber > total_number_vl:
            print 'the number of fiber is more than the virtual locations'
            sys.exit(0)
        print 'the total number of virtual location is', total_number_vl
        #calculate the postion of each virtual location
        self.x_position_vl = []
        self.y_position_vl = []
        for i_circle in sp.arange(self.total_circles):
            if i_circle == 0:
                self.x_position_vl.append(0.)
                self.y_position_vl.append(0.)
            else:
                each_circle = self.number_circle[i_circle]
                for i_position in sp.arange(each_circle):
                    x_position = self.radius_circle[i_circle] * sp.cos(2 * i_position * sp.arcsin(self.radius_fiber[0] \
                                / self.radius_circle[i_circle]))
                    y_position = self.radius_circle[i_circle] * sp.sin(2 * i_position * sp.arcsin(self.radius_fiber[0] \
                                /self.radius_circle[i_circle]))
                    self.x_position_vl.append(x_position)
                    self.y_position_vl.append(y_position)
        #calculate the distribution value in each circle zone
        self.probability_value = sp.empty(self.total_circles, float)
        for i_num_fiber in sp.arange(self.total_circles):
            self.probability_value[i_num_fiber] = 0.8 #(1 - 2 * self.theta_value) * sp.power((sp.exp(1) - \
                                                #sp.exp(self.radius_circle[i_num_fiber] / self.radius_yarn))\
                                               #/(sp.exp(1) - 1), self.beta_value) + self.theta_value
        print 'the value of probability', self.probability_value
        #calculate how many fibers can be in each circle zone
        self.each_circle_zone_num = sp.empty(self.total_circles, int)
        for i_num_fiber in sp.arange(self.total_circles):
            self.each_circle_zone_num[i_num_fiber] = int(round(self.number_circle[i_num_fiber] * \
                                                self.probability_value[i_num_fiber]))
        total_number_fiber = sum(self.each_circle_zone_num[:])
        print 'total_number', total_number_fiber
        if total_number_fiber < self.number_fiber:
            print 'it is impossible to reach the required number of fiber'
            sys.exit(0)
        #distribute the fiber in the yarn:
        if regenerate:
            self.circle_file = open(filepath, "w")
            self.fib_centers_x = open(filepath1, "w")
            self.fib_centers_y = open(filepath2, "w")
            self.x_central = 0.
            self.y_central = 0.
            self.z = 0.
            index = 1
            #record the gmesh data of yarn
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
            index = index + 5
            #record the gmesh data of the fiber
            number_fiber_in_loop = self.number_fiber
            self.circle_loop = 0
            i_determine = 0
            while number_fiber_in_loop > 0:
                number_fiber_in_loop = number_fiber_in_loop - self.each_circle_zone_num[i_determine]
                i_determine += 1
                self.circle_loop += 1
            print 'the number of loop for distributing fibers', self.circle_loop
            #sys.exit(0)
            determine_generate = 0 #the number of generated fiber
            index_position = 0
            i_circle_number = 0
            number_fiber_in_loop = self.number_fiber
            while i_circle_number < self.circle_loop: #and determine_generate < self.number_fiber:
                #for i_circle_number in sp.arange(self.total_circles):
                if i_circle_number == 0:
                    self.x_position[determine_generate] = self.x_central
                    self.y_position[determine_generate] = self.y_central
                    self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                    self.x_central, self.y_central, self.z, 
                                    self.cellsize_centre))
                    self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                    self.x_central - self.radius_fiber[0], self.y_central, 
                                    self.z, self.cellsize_centre))
                    self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                    self.x_central, self.y_central + self.radius_fiber[0],
                                    self.z, self.cellsize_centre))
                    self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                    self.x_central + self.radius_fiber[0], self.y_central,
                                    self.z, self.cellsize_centre))                      
                    self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                    self.x_central, self.y_central - self.radius_fiber[0],
                                    self.z, self.cellsize_centre))                    
                    index += 5
                    index_position += 1
                    i_circle_number += 1
                    determine_generate += 1
                    number_fiber_in_loop = number_fiber_in_loop - 1
                elif i_circle_number == 1:
                    for i_index in sp.arange(self.number_circle[i_circle_number]):
                        self.x_central = self.x_position_vl[i_index + 1]
                        self.y_central = self.y_position_vl[i_index + 1]
                        self.x_position[determine_generate] = self.x_central
                        self.y_position[determine_generate] = self.y_central
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                        self.x_central , self.y_central, self.z, 
                                        self.cellsize_centre))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                        self.x_central - self.radius_fiber[0], self.y_central, 
                                        self.z, self.cellsize_centre))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                        self.x_central, self.y_central + self.radius_fiber[0],
                                        self.z, self.cellsize_centre))
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                        self.x_central + self.radius_fiber[0], self.y_central,
                                        self.z, self.cellsize_centre))                      
                        self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                        self.x_central, self.y_central - self.radius_fiber[0],
                                        self.z, self.cellsize_centre))
                        index_position += 1
                        index += 5
                        determine_generate += 1
                    i_circle_number += 1
                    number_fiber_in_loop = number_fiber_in_loop - 6
                    #determine_generate += self.number_circle[i_circle_number]
                    #print 'the number of circle value', determine_generate
                elif i_circle_number > 1 and i_circle_number < self.circle_loop - 1:
                    #location_number = []
                    print 'i_circle_number value', i_circle_number
                    location_number = sp.empty(self.each_circle_zone_num[i_circle_number], int)
                    #each_zone = 0
                    for i_index in sp.arange(self.each_circle_zone_num[i_circle_number]):
                        if i_index == 0:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[i_circle_number] )
                            random_position = int(a_position)
                            location_number[i_index] = random_position#.append(random_position)
                            a = random_position - 1 
                            b = random_position - 1
                            self.x_central = self.x_position_vl[a]#index_position + random_position]
                            self.y_central = self.y_position_vl[b]#index_position + random_position]
                            self.x_position[determine_generate] = self.x_central
                            self.y_position[determine_generate] = self.y_central
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                            self.x_central , self.y_central, self.z, 
                                            self.cellsize_centre))
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                            self.x_central - self.radius_fiber[0], self.y_central, 
                                            self.z, self.cellsize_centre))
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                            self.x_central, self.y_central + self.radius_fiber[0],
                                            self.z, self.cellsize_centre))
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                            self.x_central + self.radius_fiber[0], self.y_central,
                                            self.z, self.cellsize_centre))                      
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                            self.x_central, self.y_central - self.radius_fiber[0],
                                            self.z, self.cellsize_centre))
                            
                            determine_generate += 1
                            index += 5
                        else:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[i_circle_number])
                            random_position = int(a_position)
                            determine_value = (random_position == location_number)
                            print determine_value
                            while determine_value.any() == True:
                                a_position = np.random.uniform(index_position, 
                                            index_position + self.number_circle[i_circle_number] + 1)
                                random_position = int(a_position)
                                determine_value = (random_position == location_number)
                            else:
                                location_number[i_index] = random_position#.append(random_position)
                                self.x_central = self.x_position_vl[random_position - 1]
                                self.y_central = self.y_position_vl[random_position - 1]
                                self.x_position[determine_generate] = self.x_central
                                self.y_position[determine_generate] = self.y_central
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                                self.x_central , self.y_central, self.z, 
                                                self.cellsize_centre))
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                                self.x_central - self.radius_fiber[0], self.y_central, 
                                                self.z, self.cellsize_centre))
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                                self.x_central, self.y_central + self.radius_fiber[0],
                                                self.z, self.cellsize_centre))
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                                self.x_central + self.radius_fiber[0], self.y_central,
                                                self.z, self.cellsize_centre))                      
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                                self.x_central, self.y_central - self.radius_fiber[0],
                                                self.z, self.cellsize_centre))
                                index += 5
                                determine_generate += 1
                    index_position += self.number_circle[i_circle_number]
                    number_fiber_in_loop = number_fiber_in_loop - self.each_circle_zone_num[i_circle_number]
                    i_circle_number += 1
                elif i_circle_number == self.circle_loop - 1:
                    location_number = sp.empty(number_fiber_in_loop, int)
                    for i_index in sp.arange(number_fiber_in_loop):
                        if i_index == 0:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[i_circle_number] )
                            random_position = int(a_position)
                            location_number[i_index] = random_position#.append(random_position)
                            a = random_position - 1 
                            b = random_position - 1
                            self.x_central = self.x_position_vl[a]#index_position + random_position]
                            self.y_central = self.y_position_vl[b]#index_position + random_position]
                            self.x_position[determine_generate] = self.x_central
                            self.y_position[determine_generate] = self.y_central
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                            self.x_central , self.y_central, self.z, 
                                            self.cellsize_centre))
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                            self.x_central - self.radius_fiber[0], self.y_central, 
                                            self.z, self.cellsize_centre))
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                            self.x_central, self.y_central + self.radius_fiber[0],
                                            self.z, self.cellsize_centre))
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                            self.x_central + self.radius_fiber[0], self.y_central,
                                            self.z, self.cellsize_centre))                      
                            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                            self.x_central, self.y_central - self.radius_fiber[0],
                                            self.z, self.cellsize_centre))
                            determine_generate += 1
                            index += 5
                        else:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[i_circle_number])
                            random_position = int(a_position)
                            determine_value = (random_position == location_number)
                            print determine_value
                            while determine_value.any() == True:
                                a_position = np.random.uniform(index_position, 
                                            index_position + self.number_circle[i_circle_number] + 1)
                                random_position = int(a_position)
                                determine_value = (random_position == location_number)
                            else:
                                location_number[i_index] = random_position#.append(random_position)
                                self.x_central = self.x_position_vl[random_position - 1]
                                self.y_central = self.y_position_vl[random_position - 1]
                                self.x_position[determine_generate] = self.x_central
                                self.y_position[determine_generate] = self.y_central
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                                self.x_central , self.y_central, self.z, 
                                                self.cellsize_centre))
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                                                self.x_central - self.radius_fiber[0], self.y_central, 
                                                self.z, self.cellsize_centre))
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                                                self.x_central, self.y_central + self.radius_fiber[0],
                                                self.z, self.cellsize_centre))
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                                                self.x_central + self.radius_fiber[0], self.y_central,
                                                self.z, self.cellsize_centre))                      
                                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                                                self.x_central, self.y_central - self.radius_fiber[0],
                                                self.z, self.cellsize_centre))
                                index += 5
                                determine_generate += 1
                    #determine_generate += self.each_circle_zone_num[i_circle_number]
                    #print 'the number of circle of generated', determine_generate
                    index_position += self.number_circle[i_circle_number]
                    i_circle_number += 1
            #generate circle and joining part for gmsh
            if self.verbose:
                print "all the points have been generated"
            index_point_in_circle = 0 #the number of each point in Circle part
            for i1 in sp.arange(0, self.number_fiber + 1, 1):
                if i1 == 0:
                    #index = index + 1
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
            #elapsed = (time.clock() - start)
            if self.verbose:
                print "Time to generate the mesh file: %.3f" %(elapsed)
                #above part is for generating the surface loop in the yarn domain
                print "Mesh generation finished"
            
            index_loop_in_plane = index
            index = index + 1
            self.circle_file.write("Plane Surface(%d) = {%d};\n" %(index, index_loop_in_plane))
            self.circle_file.close()
            self.fib_centers_x.write(repr(self.x_position))
            self.fib_centers_y.write(repr(self.y_position))
            self.fib_centers_x.close()
            self.fib_centers_y.close()
            circledef = open(filepath, "r").readlines()
        else:
            circledef = open(filepath, "r").readlines()
        return ''.join(circledef)
        
            
    def mesh_new_generate(self, filename='yarn_new.geo', filename1 = 'fib_centers_x_new.geo',
        filename2 = 'fib_centers_y_new.geo' ,regenerate=True):
        self.mesh = Gmsh2D(self.create_circle_domain_gmesh(filename, filename1, filename2,
                            regenerate))
        return self.mesh
        
        
        