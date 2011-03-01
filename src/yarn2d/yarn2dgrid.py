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
        self.number_fiber_blend = [val/100*self.number_fiber for val in self.blend]
        self.number_fiber_blend[-1] = self.number_fiber - sum(self.number_fiber_blend[:-1])
        print 'fibers per blend', self.number_fiber_blend,' total', self.number_fiber
        
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

    if len(options.radius_fiber) > 1:
        print 'ERROR: randomfiberlayout can only handle one type of fiber'
        sys.exit()

    #allocate space
    x_position = sp.empty(options.number_fiber, float)
    y_position = sp.empty(options.number_fiber, float)
    radius_fiber = sp.empty(options.number_fiber, float)
    fiber_kind = sp.zeros(options.number_fiber, int)

    for i in sp.arange(0, options.number_fiber, 1):
        #generate the position of fiber
        a = np.random.uniform(-1, 1)
        b = np.random.uniform(-1, 1)
        #distance between the current point and center
        distance_center = sp.sqrt((a - x_central)**2 + (b - y_central)**2)
        #distance between the current point and existing points
        distance_each = sp.sqrt((a - x_position[:i])**2 + \
                                (b - y_position[:i])**2)
        while (distance_center + options.radius_fiber[0] >= options.radius_yarn or 
               (i>0 and (np.min(distance_each)\
                <= (2*options.radius_fiber[0])+ 0.001 * options.radius_fiber[0]))): 
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
            radius_fiber[i] = options.radius_fiber[0]
    return (x_position, y_position, radius_fiber, fiber_kind)

def virtloclayout(options):
    raise NotImplementedError

def virtlocoverlaplayout(options):
    raise NotImplementedError
