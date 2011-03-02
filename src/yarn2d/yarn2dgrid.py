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
from yarn2d.config import FIBERLAYOUTS

#-------------------------------------------------------------------------
#
# Yarn2dGrid class
#
#-------------------------------------------------------------------------

class Yarn2dGrid(object):
    def __init__(self, cfg):
        self.cfg = cfg
        self.fiberlayout = self.cfg.get('domain.fiberlayout_method')
        if not (self.fiberlayout in FIBERLAYOUTS):
            print 'ERROR: unkown fiber layout method %s' % self.fiberlayout
            sys.exit(0)
        
        self.Ry = self.cfg.get('domain.yarnradius')
        self.Rf = self.cfg.get('fiber.radius_fiber')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        #computational radius
        self.radius_yarn = self.scaleL * self.Ry
        self.radius_fiber =  [self.scaleL * rad for rad in self.Rf]
        self.radius_boundlayer = max(self.radius_fiber)/2.
        self.radius_domain = self.radius_yarn + self.radius_boundlayer
        self.cellsize_centre = self.cfg.get('domain.cellsize_centre')
        self.cellSize = self.cfg.get('domain.cellsize_fiber')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        self.number_fiber_blend = [int(round(val/100*self.number_fiber)) for val in self.blend]
        self.number_fiber_blend[-1] = self.number_fiber - sum(self.number_fiber_blend[:-1])
        print 'fibers per blend', self.number_fiber_blend,' total', self.number_fiber
        self.theta_value = self.cfg.get('domain.theta_value')
        self.beta_value = self.cfg.get('domain.beta_value')
        
        self.verbose = self.cfg.get('general.verbose')

    def create_circle_domain_gmsh(self, filename='yarn.geo', 
            layoutfile = 'layout.dat', regenerate=True):
        """
        Create gmsh file with the circle domain of yarn and the fibers in it
        returns string with defenition, path to file
        Layout is dumped to layout.dat
        """
        filepath = utils.OUTPUTDIR + os.sep + filename
        start = time.clock()
        self.x_central = 0.
        self.y_central = 0.
        self.z = 0.
        
        if not regenerate:
            npzfile = np.load(utils.OUTPUTDIR + os.sep + layoutfile)
            self.x_position = npzfile['fib_centers_x'], 
            self.y_position = npzfile['fib_centers_y'],
            self.all_radius_fibers = npzfile['fib_radius'],
            self.fiber_kind = npzfile['fiber_kind']
            circledef = open(filepath, "r").readlines()
        else:
            #first determine the different centerpoints of the fibers
            if self.fiberlayout == 'random':
                layoutfun = randomfiberlayout
            elif self.fiberlayout == 'virtloc':
                layoutfun = virtloclayout
            elif self.fiberlayout == 'virtlocoverlap':
                layoutfun = virtlocoverlaplayout
            else:
                print 'ERROR: not implemented fiberlayout %s' % self.fiberlayout
                sys.exit(0)
            
            self.x_position, self.y_position, self.all_radius_fibers, \
                    self.fiber_kind = layoutfun(self)
    
            #write data to files
            self.circle_file = open(filepath, "w")
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
                index = index + 1
                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                    self.x_position[i-1], self.y_position[i-1], 
                                    self.z,self.cellSize))
                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                        self.x_position[i-1] - self.all_radius_fibers[i-1],
                        self.y_position[i-1], self.z, self.cellSize))
                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                        self.x_position[i-1], 
                        self.y_position[i-1] + self.all_radius_fibers[i-1],
                        self.z, self.cellSize))
                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                        self.x_position[i-1] + self.all_radius_fibers[i-1],
                        self.y_position[i-1], self.z, self.cellSize))
                self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                        self.x_position[i-1], 
                        self.y_position[i-1] - self.all_radius_fibers[i-1],
                        self.z, self.cellSize))
                index = index + 4
            #above part is for generating the points of circle in the domain
            if self.verbose:
                print "all the points have been generated"
            index_point_in_circle = 0 #the number of each point in Circle part
            for i1 in sp.arange(0, self.number_fiber + 1, 1):
                index = index + 1
                if i1 == 0:
                    index_point_in_circle = index_point_in_circle + 1
                else:
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
            for i2 in sp.arange(0, 4*(self.number_fiber + 1)-1, 1):
                index_circle_for_loop +=  1
                self.circle_file.write("%d," %(index_circle_for_loop))
            index_circle_for_loop += 1
            self.circle_file.write("%d};\n" %(index_circle_for_loop))
            
            elapsed = (time.clock() - start)
            if self.verbose:
                print "Time to generate the mesh file: %.3f" %(elapsed)
                #above part is for generating the surface loop in the yarn domain
                print "Mesh generation finished"
            
            index_loop_in_plane = index
            index = index + 1
            self.circle_file.write("Plane Surface(%d) = {%d};\n" %(index, index_loop_in_plane))
            #flush to disk - slow!
            self.circle_file.flush()
            os.fsync(self.circle_file)
            self.circle_file.close()
            #dump the layout to disk
            self.x_position, self.y_position, self.all_radius_fibers, \
                    self.fiber_kind
            np.savez(utils.OUTPUTDIR + os.sep + layoutfile, 
                     fib_centers_x=self.x_position, 
                     fib_centers_y=self.y_position,
                     fib_radius=self.all_radius_fibers,
                     fiber_kind=self.fiber_kind)
            circledef = open(filepath, "r").readlines()
        return '\n'.join(circledef)
    
    def mesh_2d_generate(self, filename='yarn.geo', regenerate=True):
        """
        Return a Gmsh2D object from Fipy for the yarn 2D grid
        The gmsh file is written to filename.
        If regenerate is True, it is looked if the file already exists, and if
        so, the file is reused instead of generated.
        """
        self.mesh = Gmsh2D(self.create_circle_domain_gmsh(filename,
                            regenerate=regenerate))
        return self.mesh

def randomfiberlayout(options):
    """ Generate the fiber layout in the yarn in random fashion
    Options should contain attributes:
        number_fiber : amount of fibers to generate
        number_fiber_blend : number of fibers per blend type
        x_central: central point, 0. otherwise
        y_central: central point, 0. otherwise
        radius_yarn
        radius_fiber: array of length 1 with radius of the fiber
    
    Returns a tuple 
        (list of x coord center points,
         list of y coord center points,
         list of radius of fiber,
         list integers indicating kind of fiber)
    """
    #check options
    if hasattr(options, 'x_central'):
        x_central = options.x_central
    else:
        x_central = 0.
    if hasattr(options, 'y_central'):
        y_central = options.y_central
    else:
        y_central = 0.

    ##if len(options.radius_fiber) > 1:
    ##    print 'ERROR: randomfiberlayout can only handle one type of fiber'
    ##    sys.exit()

    #allocate space
    x_position = sp.empty(options.number_fiber, float)
    y_position = sp.empty(options.number_fiber, float)
    radius_fiber = sp.empty(options.number_fiber, float)
    fiber_kind = sp.empty(options.number_fiber, int)
    
    #preset the radius_fiber. This algorithm does first the largest fibers,
    #then the others
    ind = sp.arange(options.number_fiber)
    radius_fiber[:options.number_fiber_blend[0]] = options.radius_fiber[0]
    fiber_kind[:options.number_fiber_blend[0]] = 0
    for j in range(1, len(options.number_fiber_blend)):
        assert options.radius_fiber[j] <= options.radius_fiber[j-1], \
                'Give fibers in decreasing size in ini file'
        radius_fiber[sum(options.number_fiber_blend[:j]):sum(options.number_fiber_blend[:j+1])]\
            = options.radius_fiber[j]
        fiber_kind[:options.number_fiber_blend[0]] = j

    for i in sp.arange(0, options.number_fiber, 1):
        #generate the position of fiber
        a = np.random.uniform(-1, 1)
        b = np.random.uniform(-1, 1)
        #distance between the current point and center
        distance_center = sp.sqrt((a - x_central)**2 + (b - y_central)**2)
        #distance between the current point and existing points
        distance_each = sp.sqrt((a - x_position[:i])**2 + \
                                (b - y_position[:i])**2)
        while (distance_center + radius_fiber[i] >= options.radius_yarn or
               (i>0 and not (( distance_each
                >= (1.001*radius_fiber[i] + radius_fiber[:i])).all()))): 
            #   (i>0 and (np.min(distance_each)\
            #    <= (2*options.radius_fiber[0])+ 0.001 * options.radius_fiber[0]))):
            a = np.random.uniform(-1, 1)
            b = np.random.uniform(-1, 1)
            distance_center = sp.sqrt((a - x_central)**2 + \
                                      (b - y_central)**2)
            #distance between the current point and existing points
            distance_each = sp.sqrt((a - x_position[:i])**2 + \
                                (b - y_position[:i])**2)
        else:
            x_position[i] = a
            y_position[i] = b
    return (x_position, y_position, radius_fiber, fiber_kind)

def virtloclayout(options):
    """ Generate the fiber layout in the yarn in virtual locations
    Options should contain attributes:
        number_fiber : amount of fibers to generate
        radius_yarn
        radius_fiber: array of length 1 with radius of the fiber
    
    Returns a tuple 
        (list of x coord center points,
         list of y coord center points,
         list of radius of fiber,
         list integers indicating kind of fiber)
    """
    number_for_circles = int((options.radius_yarn - options.radius_fiber[0]) / (2 * 
                                    options.radius_fiber[0]))
    total_circles = number_for_circles + 1
    radius_circle = sp.empty(total_circles, float)
    radius_circle_1 = sp.empty(total_circles, float)
    x_position = sp.empty(options.number_fiber,float)
    y_position = sp.empty(options.number_fiber,float)
    number_circle = sp.empty(total_circles, int)
    radius_fiber = sp.empty(options.number_fiber, float)
    fiber_kind = sp.empty(options.number_fiber, int)
    
    ind = sp.arange(options.number_fiber)
    radius_fiber[:options.number_fiber_blend[0]] = options.radius_fiber[0]
    fiber_kind[:options.number_fiber_blend[0]] = 0
    for j in range(1, len(options.number_fiber_blend)):
        assert options.radius_fiber[j] <= options.radius_fiber[j-1], \
                'Give fibers in decreasing size in ini file'
        radius_fiber[sum(options.number_fiber_blend[:j]):sum(options.number_fiber_blend[:j+1])]\
            = options.radius_fiber[j]
        fiber_kind[:options.number_fiber_blend[0]] = j
        
    for i_circle in sp.arange(total_circles):
        if i_circle == 0:
            radius_circle[i_circle] = (i_circle + 1) * (options.radius_fiber[0] + 0.1 * 
                            options.radius_fiber[0])
            radius_circle_1[i_circle] = (i_circle + 1) * options.radius_fiber[0] * 1.1
            number_circle[i_circle] = 1
        elif i_circle > 0:
            radius_circle[i_circle] = i_circle * 2. * (options.radius_fiber[0] + 0.1 * 
                            options.radius_fiber[0])
            radius_circle_1[i_circle] = i_circle * 2. * options.radius_fiber[0] * 1.1
            number_circle[i_circle] = int(sp.pi / (1.01 * options.radius_fiber[0] 
                                            / radius_circle_1[i_circle]))
            #print 'the value', self.radius_fiber[0] /self.radius_circle[i_circle]
            #print 'the number in each circle', self.number_circle[i_circle]
    total_number_vl = sum(number_circle[:])
    if options.number_fiber > total_number_vl:
        print 'ERROR: the number of fiber is more than the virtual locations'
        sys.exit(0)
    #calculate the postion of each virtual location
    x_position_vl = []
    y_position_vl = []
    probability_value = sp.empty(total_circles, float)
    each_circle_zone_num = sp.empty(total_circles, int)
    for i_circle in sp.arange(total_circles):
        if i_circle == 0:
            x_position_vl.append(0.)
            y_position_vl.append(0.)
        else:
            each_circle = number_circle[i_circle]
            for i_position in sp.arange(each_circle):
                x_position_t = radius_circle[i_circle] * sp.cos(2 * i_position * 
                            sp.arcsin(1.01 * options.radius_fiber[0] / 
                            radius_circle_1[i_circle]))
                y_position_t = radius_circle[i_circle] * sp.sin(2 * i_position * 
                            sp.arcsin(1.01 * radius_fiber[0] /
                            radius_circle_1[i_circle]))
                x_position_vl.append(x_position_t)
                y_position_vl.append(y_position_t)
        #calculate the distribution value in each circle zone
        probability_value[i_circle] = (1 - 2 * options.theta_value) * sp.power(
                                    (sp.exp(1.) -sp.exp(radius_circle[i_circle] / 
                                    options.radius_yarn))/(sp.exp(1) - 1), 
                                    options.beta_value) + options.theta_value        
        each_circle_zone_num[i_circle] = int(round(number_circle[i_circle] * \
                                                probability_value[i_circle]))
        total_number_fiber = sum(each_circle_zone_num[:])
    #distribute the fiber in the yarn:
    if hasattr(options, 'x_central'):
        x_central = options.x_central
    else:
        x_central = 0.
    if hasattr(options, 'y_central'):
        y_central = options.y_central
    else:
        y_central = 0.
    number_fiber_in_loop = options.number_fiber
    circle_loop = 0
    i_determine = 0
    while number_fiber_in_loop > 0:
        number_fiber_in_loop = number_fiber_in_loop - each_circle_zone_num[i_determine]
        i_determine += 1
        circle_loop += 1
    print 'the number of loop for distributing fibers', circle_loop
    
    determine_generate = 0 #the number of generated fiber
    index_position = 0
    i_circle_number = 0
    number_fiber_in_loop = options.number_fiber
    while i_circle_number < circle_loop: 
        if i_circle_number == 0:
            x_position[determine_generate] = 0.0
            y_position[determine_generate] = 0.0
            index_position += 1
            i_circle_number += 1
            determine_generate += 1
            number_fiber_in_loop = number_fiber_in_loop - 1
        elif i_circle_number == 1:
            for i_index in sp.arange(number_circle[i_circle_number]):
                x_position[determine_generate] = x_position_vl[i_index + 1]
                y_position[determine_generate] = y_position_vl[i_index + 1]
                index_position += 1
                determine_generate += 1
            i_circle_number += 1
            number_fiber_in_loop = number_fiber_in_loop - 6
        elif i_circle_number > 1 and i_circle_number < circle_loop - 1:
            location_number = sp.empty(each_circle_zone_num[i_circle_number], int)
            for i_index in sp.arange(each_circle_zone_num[i_circle_number]):
                if i_index == 0:
                    a_position = np.random.uniform(index_position, 
                                index_position + number_circle[i_circle_number])
                    random_position = int(a_position)
                    location_number[i_index] = random_position
                    x_position[determine_generate] = x_position_vl[random_position]
                    y_position[determine_generate] = y_position_vl[random_position]
                    determine_generate += 1
                else:
                    a_position = np.random.uniform(index_position, 
                                index_position + number_circle[i_circle_number])
                    random_position = int(a_position)
                    determine_value = (random_position == location_number)
                    while determine_value.any() == True:
                        a_position = np.random.uniform(index_position, 
                                    index_position + number_circle[i_circle_number])
                        random_position = int(a_position)
                        determine_value = (random_position == location_number)
                    else:
                        location_number[i_index] = random_position
                        x_position[determine_generate] = x_position_vl[random_position]
                        y_position[determine_generate] = y_position_vl[random_position]
                        determine_generate += 1
            index_position += number_circle[i_circle_number]
            number_fiber_in_loop = number_fiber_in_loop - each_circle_zone_num[i_circle_number]
            i_circle_number += 1
            
        elif i_circle_number == circle_loop - 1:
            location_number = sp.empty(number_fiber_in_loop, int)
            for i_index in sp.arange(number_fiber_in_loop):
                if i_index == 0:
                    a_position = np.random.uniform(index_position, index_position +
                                number_circle[i_circle_number])
                    random_position = int(a_position)
                    location_number[i_index] = random_position
                    x_position[determine_generate] = x_position_vl[random_position]
                    y_position[determine_generate] = y_position_vl[random_position]
                    determine_generate += 1
                else:
                    a_position = np.random.uniform(index_position, index_position +
                                number_circle[i_circle_number])
                    random_position = int(a_position)
                    determine_value = (random_position == location_number)
                    while determine_value.any() == True:
                        a_position = np.random.uniform(index_position, 
                                    index_position + number_circle[i_circle_number])
                        random_position = int(a_position)
                        determine_value = (random_position == location_number)
                    else:
                        location_number[i_index] = random_position
                        x_position[determine_generate] = x_position_vl[random_position]
                        y_position[determine_generate] = y_position_vl[random_position]
                        determine_generate += 1
            index_position += number_circle[i_circle_number]
            i_circle_number += 1
    return (x_position, y_position, radius_fiber, fiber_kind)
    
def virtlocoverlaplayout(options):
    raise NotImplementedError
