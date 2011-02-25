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
# Yarn2dGrid new class to generate the domain including two kinds of fibers
# with different value of radius
#-------------------------------------------------------------------------

class Yarn2dNewGridRadius(object):
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
        
        self.type_fiber = self.cfg.get('fiber.number_type')
        self.blend = self.cfg.get('fiber.blend')

    def create_circle_domain_gmsh(self, filename='yarn_r_new.geo', filename1 = 'fib_centers_x_new_r.geo',
        filename2 = 'fib_centers_y_new_r.geo' ,regenerate=True):
        filepath = utils.OUTPUTDIR + os.sep + filename
        filepath1 = utils.OUTPUTDIR + os.sep + filename1
        filepath2 = utils.OUTPUTDIR + os.sep + filename2
        self.number_for_circles = sp.empty(self.type_fiber)
        self.total_circles = sp.empty(self.type_fiber)
        
        #the number of layers for each kind of fiber
        self.number_for_circles[:] = int((self.radius_yarn - self.radius_fiber[:]) / (2 * \
                                    self.radius_fiber[:]))
        self.total_circles[:] = self.number_for_circles[:] + 1
        #the radius of each layer for fiber_1
        self.radius_circle = sp.empty(self.total_circles[0], float)
        self.radius_circle_1 = sp.empty(self.total_circles[0], float)
        for i_circles in sp.arange(self.total_circles[0]):
            if i_circle == 0:
                self.radius_circle[i_circle] = (i_circle + 1) * (self.radius_fiber[0] + 0.1 * \
                                self.radius_fiber[0])
                self.radius_circle_1[i_circle] = (i_circle + 1) * self.radius_fiber[0]
            elif i_circle > 0:
                self.radius_circle[i_circle] = i_circle * 2. * (self.radius_fiber[0] + 0.1 * \
                                self.radius_fiber[0])
                self.radius_circle_1[i_circle] = i_circle * 2. * self.radius_fiber[0]
        #the radius of each layer for fiber_2
        self.radius_circle1 = sp.empty(self.total_circles[1], float)
        self.radius_circle_2 = sp.empty(self.total_circles[1], float)
        for i_circles in sp.arange(self.total_circles[1]):
            if i_circle == 0:
                self.radius_circle1[i_circle] = (i_circle + 1) * (self.radius_fiber[1] + 0.1 * \
                                self.radius_fiber[1])
                self.radius_circle_2[i_circle] = (i_circle + 1) * self.radius_fiber[1]
            elif i_circle > 0:
                self.radius_circle1[i_circle] = i_circle * 2. * (self.radius_fiber[1] + 0.1 * \
                                self.radius_fiber[1])
                self.radius_circle_2[i_circle] = i_circle * 2. * self.radius_fiber[1]
        """
        generate the virtual locations
        """
        #calculate the number of virtual location
        self.number_circle_1 = sp.empty(self.total_circles, int) #the number of virtual locations for fiber_1
        for i_circle in sp.arange(self.total_circles[0]):
            if i_circle == 0:
                self.number_circle_1[i_circle] = 1
            else:
                self.number_circle_1[i_circle] = int(sp.pi / (self.radius_fiber[0] \
                                                / self.radius_circle_1[i_circle]))
                print 'the value', self.radius_fiber[0] /self.radius_circle[i_circle]
                print 'the number in each circle', self.number_circle[i_circle]
        total_number_vl_1 = sum(self.number_circle_1[:])
        
        self.number_circle_2 = sp.empty(self.total_circles, int) #the number of virtual locations for fiber_2
        for i_circle in sp.arange(self.total_circles[1]):
            if i_circle == 0:
                self.number_circle_2[i_circle] = 1
            else:
                self.number_circle_2[i_circle] = int(sp.pi / (self.radius_fiber[1] \
                                                / self.radius_circle_2[i_circle]))
                print 'the value', self.radius_fiber[0] /self.radius_circle1[i_circle]
                print 'the number in each circle', self.number_circle_2[i_circle]
        total_number_vl_2 = sum(self.number_circle_2[:])
        
        #calculate the postion of each virtual location for fiber_1
        self.x_position_vl_1 = []
        self.y_position_vl_1 = []
        for i_circle in sp.arange(self.total_circles[0]):
            if i_circle == 0:
                self.x_position_vl_1.append(0.)
                self.y_position_vl_1.append(0.)
            else:
                each_circle = self.number_circle_1[i_circle]
                for i_position in sp.arange(each_circle):
                    x_position = self.radius_circle[i_circle] * sp.cos(2 * i_position * sp.arcsin(self.radius_fiber[0] \
                                / self.radius_circle_1[i_circle]))
                    y_position = self.radius_circle[i_circle] * sp.sin(2 * i_position * sp.arcsin(self.radius_fiber[0] \
                                /self.radius_circle_1[i_circle]))
                    self.x_position_vl_1.append(x_position)
                    self.y_position_vl_1.append(y_position)
        #calculate the position of each virtual location for fiber_2
        self.x_position_vl_2 = []
        self.y_position_vl_2 = []
        for i_circle in sp.arange(self.total_circles[1]):
            if i_circle ==0:
                self.x_position_vl_2.append(0.)
                self.y_position_vl_2.append(0.)
            else:
                each_circle = self.number_circle_2[i_circle]
                for i_position in sp.arange(each_circle):
                    x_position = self.radius_circle_1[i_circle] * sp.cos(2 * i_position * sp.arcsin(self.radius_fiber[1] \
                                / self.radius_circle_1[i_circle]))
                    y_position = self.radius_circle_1[i_circle] * sp.sin(2 * i_position * sp.arcsin(self.radius_fiber[1] \
                                /self.radius_circle_1[i_circle]))
                    self.x_position_vl_2.append(x_position)
                    self.y_position_vl_2.append(y_position)
        #calculate the probability value of fiber in each circle zone for fiber_1
        self.probability_value_1 = sp.empty(self.total_circles[0], float)
        for i_num_fiber in sp.arange(self.total_circles[0]):
            self.probability_value_1[i_num_fiber] = (1 - 2 * self.theta_value) * sp.power((sp.exp(1) - \
                                                sp.exp(self.radius_circle[i_num_fiber] / self.radius_yarn))\
                                               /(sp.exp(1) - 1), self.beta_value) + self.theta_value
        print 'the value of probability', self.probability_value_1
        #calculate the probability value of fiber in each circle zone for fiber_2
        self.probability_value_2 = sp.empty(self.total_circles[1], float)
        for i_num_fiber in sp.arange(self.total_circles[1]):
            self.probability_value_2[i_num_fiber] = (1 - 2 * self.theta_value) * sp.power((sp.exp(1) - \
                                                sp.exp(self.radius_circle1[i_num_fiber] / self.radius_yarn))\
                                               /(sp.exp(1) - 1), self.beta_value) + self.theta_value
        print 'the value of probability', self.probability_value_2
        
        #calculate how many fibers can be in each circle zone for fiber_1
        self.each_circle_zone_num1 = sp.empty(self.total_circles[0], int)
        for i_num_fiber in sp.arange(self.total_circles[0]):
            self.each_circle_zone_num1[i_num_fiber] = int(round(self.number_circle_1[i_num_fiber] * \
                                                self.probability_value_1[i_num_fiber]))
        total_number_fiber_1 = sum(self.each_circle_zone_num1[:])
        print 'total_number', total_number_fiber_1
        #calculate how many fibers can be in each circle zone for fiber_2
        self.each_circle_zone_num2 = sp.empty(self.total_circles[1], int)
        for i_num_fiber in sp.arange(self.total_circles[1]):
            self.each_circle_zone_num2[i_num_fiber] = int(round(self.number_circle_2[i_num_fiber] * \
                                                self.probability_value_2[i_num_fiber]))
        total_number_fiber_2 = sum(self.each_circle_zone_num2[:])
        print 'total_number', total_number_fiber_1
        total_number_fiber = total_number_fiber_1 + total_number_fiber_2
        if total_number_fiber < self.number_fiber:
            print 'it is impossible to reach the required number of fiber'
            sys.exit(0)
        
        #distribute two kinds of fibers in the yarn domain
        
        
        

        