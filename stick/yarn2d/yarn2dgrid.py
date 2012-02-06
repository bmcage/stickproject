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
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import pylab
import matplotlib
import sympy
from sympy.abc import x,y

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
from fipy import Gmsh2D
from fipy import *
from yarn.config import FIBERLAYOUTS
from fiber.config import FiberConfigManager, FIBER_FORM, CIRCLE, ELLIPSE

from virtlocgeom import *
from fiber_layout import *

#-------------------------------------------------------------------------
#
# Yarn2dGrid class
#
#-------------------------------------------------------------------------


NONTOUCH_FAC = 1.01
    
class Yarn2dGrid(object):
    def __init__(self, cfg):
        self.cfg = cfg
        self.verbose = self.cfg.get('general.verbose')
        self.fiberlayout = self.cfg.get('domain.fiberlayout_method')
        
        self.x_central = 0.
        self.y_central = 0.
        self.z = 0.
        
        if not (self.fiberlayout in FIBERLAYOUTS):
            print 'ERROR: unkown fiber layout method %s' % self.fiberlayout
            sys.exit(0)
        self.Ry = self.cfg.get('domain.yarnradius')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        #computational radius
        self.radius_yarn = self.scaleL * self.Ry
        self.cellsize_centre = self.cfg.get('domain.cellsize_centre')
        self.cellSize = self.cfg.get('domain.cellsize_fiber')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        if self.verbose:
            print 'Blend used:', len(self.blend)
        self.number_fiber_blend = [int(round(val/100*self.number_fiber)) for val in self.blend]
        self.number_fiber_blend[-1] = self.number_fiber - np.sum(self.number_fiber_blend[:-1])
        if self.verbose:
            print 'Fibers per blend', self.number_fiber_blend,' total', self.number_fiber
        self.theta_value = self.cfg.get('domain.theta_value')
        self.beta_value =  self.cfg.get('domain.beta_value')
        
        #obtain size of fibers
        self.Rf = []
        self.Rf_form = []
        self.ecc = []
        self.mean_deviation = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            cfg_fiber = FiberConfigManager.get_instance(filename)
            self.Rf.append(cfg_fiber.get('fiber.radius_pure_fiber'))
            for i in range(cfg_fiber.get('fiber.nrlayers')):
                section = 'fiberlayer_%i' % i
                self.Rf[-1] += cfg_fiber.get(section + '.thickness')
            #scaled mean deviation of the fiber
            self.mean_deviation.append(self.scaleL * cfg_fiber.get('fiber.mean_deviation'))
            self.Rf_form.append(FIBER_FORM[cfg_fiber.get('fiber.form')])
            self.ecc.append(cfg_fiber.get('fiber.eccentricity'))
        self.radius_fiber =  [self.scaleL * rad for rad in self.Rf]
        self.radius_boundlayer = max(self.radius_fiber)/2.
        self.radius_domain = self.radius_yarn + self.radius_boundlayer
        self.prob_area = eval (self.cfg.get('fiber.prob_area'))       
        
    def create_circle_domain_gmsh(self, filename='yarn.geo', 
            layoutfile = 'layout.dat', regenerate=True, plotyarn=False):
        """
        Create gmsh file with the circle domain of the yarn and the fibers in it
        Returns string with definition
        Gmsh geo is stored in OUTPUTDIR/filename
        Layout is stored in OUTPUTDIR/layout.dat as a numpy.dump
        If regenerate == True, then filename and layout are used to obtain the
        layout and the geo definition
        """
        filepath = utils.OUTPUTDIR + os.sep + filename
        #filepath_f = utils.OUTPUTDIR + os.sep +figure_file
        start = time.clock()
        
        if not regenerate:
            npzfile = np.load(utils.OUTPUTDIR + os.sep + layoutfile)
            self.x_position = npzfile['fib_centers_x'], 
            self.y_position = npzfile['fib_centers_y'],
            self.all_radius_fibers = npzfile['fib_radius'],
            self.fiber_kind = npzfile['fiber_kind']
            circledef = open(filepath, "r").readlines()
            return '\n'.join(circledef)

        #first set options that define the layout algorithm
        ouroptions = {
                'x_central' : self.x_central,
                'y_central' : self.y_central,
                'number_fiber' : self.number_fiber,
                'number_fiber_blend' : self.number_fiber_blend,
                'radius_fiber' : self.radius_fiber,
                'radius_yarn' : self.radius_yarn,
                }
        if self.fiberlayout == 'random':
            layoutfun = randomfiberlayout
        elif self.fiberlayout in ['virtloc', 'virtlocoverlap']:
            ouroptions.update({
                'theta_value' : self.theta_value,
                'beta_value' : self.beta_value,
                'mean_deviation': self.mean_deviation,
                'prob_area': self.prob_area
                })
            ouroptions['radius_first_center'] = self.cfg.get(
                                    'domain.radius_first_center_virtloc')
            if self.fiberlayout == 'virtloc':
                layoutfun = virtloclayout
            else:
                layoutfun = virtlocoverlaplayout
        else:
            print 'ERROR: not implemented fiberlayout %s' % self.fiberlayout
            sys.exit(0)
        self.x_position, self.y_position, self.all_radius_fibers, \
                    self.fiber_kind = layoutfun(ouroptions)
        #write data to files
        self.circle_file = open(filepath, "w")
        #The outer boundary
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
        nrcircle = 0
        nrellipse = 0
        for ind, number in enumerate(self.number_fiber_blend):
            if self.Rf_form[ind] == CIRCLE:
                nrcircle += number
            elif self.Rf_form[ind] == ELLIPSE:
                nrellipse += number
        index_circle = sp.empty(nrcircle)
        i_circle_p = 0
        index_ellipse = sp.empty(nrellipse)
        i_ellipse = 0
        #begin to account the number of each kind of fiber cross-section
        for i_shape in sp.arange(0, len(self.fiber_kind)):
            if self.Rf_form[self.fiber_kind[i_shape]] == CIRCLE:
                index_circle[i_circle_p] = i_shape
                i_circle_p += 1
            else:
                index_ellipse[i_ellipse] = i_shape
                i_ellipse += 1
        #generate all the points for circle cross-section
        for i_position_circle in index_circle:
            #center point and 4 points
            index += 1
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index,
                                self.x_position[i_position_circle], 
                                self.y_position[i_position_circle], 
                                self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+1,
                    self.x_position[i_position_circle] - \
                    self.all_radius_fibers[i_position_circle],
                    self.y_position[i_position_circle], self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+2,
                    self.x_position[i_position_circle], 
                    self.y_position[i_position_circle] + \
                    self.all_radius_fibers[i_position_circle],
                    self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+3,
                    self.x_position[i_position_circle] + \
                    self.all_radius_fibers[i_position_circle],
                    self.y_position[i_position_circle], self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g,%g,%g,%g};\n" %(index+4,
                    self.x_position[i_position_circle], 
                    self.y_position[i_position_circle] - \
                    self.all_radius_fibers[i_position_circle],
                    self.z, self.cellSize))
            index = index + 4
        #generate all the points for ellipse cross-section
        for i_position_ellipse in index_ellipse:
            index += 1
            #rotate each ellipse randomly
            rotate_theta = np.random.uniform(0,sp.pi * 2)
            #begin to write the points
            long_axis = self.all_radius_fibers[i_position_ellipse]
            short_axis = sqrt(self.all_radius_fibers[i_position_ellipse]**2 \
                    * (1-self.ecc[self.fiber_kind[i_position_ellipse]]**2))
            x_long_axis = long_axis * sp.cos(rotate_theta)
            y_long_axis = long_axis * sp.sin(rotate_theta)
            x_short_axis = short_axis * sp.sin(rotate_theta)
            y_short_axis = short_axis * sp.cos(rotate_theta)
            
            self.circle_file.write("Point(%d) = {%g, %g, %g, %g};\n" %(index, 
                                self.x_position[i_position_ellipse], 
                                self.y_position[i_position_ellipse], 
                                self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g, %g, %g, %g};\n" %(index + 1,
                                self.x_position[i_position_ellipse] + \
                                x_long_axis, self.y_position[i_position_ellipse] + 
                                y_long_axis, self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g, %g, %g, %g};\n" %(index + 2,
                                self.x_position[i_position_ellipse] - \
                                x_short_axis, self.y_position[i_position_ellipse] +
                                y_short_axis, self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g, %g, %g, %g};\n" %(index + 3, 
                                self.x_position[i_position_ellipse] - 
                                x_long_axis, self.y_position[i_position_ellipse] -
                                y_long_axis, self.z, self.cellSize))
            self.circle_file.write("Point(%d) = {%g, %g, %g, %g};\n" %(index + 4,
                                self.x_position[i_position_ellipse] +
                                x_short_axis, self.y_position[i_position_ellipse] -
                                y_short_axis, self.z, self.cellSize))
            index = index + 4
        index_point_in_circle = 0
        for i1 in sp.arange(0, len(index_circle) + 1, 1):
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
        index_point_in_ellipse = index_point_in_circle
        for i2 in sp.arange(0, len(index_ellipse), 1):
            index += 1
            index_point_in_ellipse += 5
            t1 = index_point_in_ellipse + 1
            t2 = index_point_in_ellipse + 2
            t3 = index_point_in_ellipse + 3
            t4 = index_point_in_ellipse + 4
            self.circle_file.write("Ellipse(%d) = {%d, %d, %d, %d};\n" %(index,
                                t1, index_point_in_ellipse, t2, t2))
            self.circle_file.write("Ellipse(%d) = {%d, %d, %d, %d};\n" %(index + 1, 
                                t2, index_point_in_ellipse, t3, t3))
            self.circle_file.write("Ellipse(%d) = {%d, %d, %d, %d};\n" %(index + 2, 
                                t3, index_point_in_ellipse, t4, t4))
            self.circle_file.write("Ellipse(%d) = {%d, %d, %d, %d};\n" %(index + 3,
                                t4, index_point_in_ellipse, t1, t1))
            index += 3
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

def test():
    #test if determine overlap works. 
    xpos = [0., 0.1, 0., 0.2]
    ypos = [0., 0., 0.1, 0.]
    rad =  [0.1, 0.05, 0.05, 0.05]
    kind = [0] * len(xpos)
    radyarn = 2.
    res = determine_overlap(xpos, ypos, rad)
    print 'res', res
    plot_yarn(xpos, ypos, rad)
    move_fibers_nonoverlap(xpos, ypos, rad, radyarn)
    plot_yarn(xpos, ypos, rad)
    pylab.show()
    

if __name__ == '__main__':
    test()
