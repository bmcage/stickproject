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
from pylab import *
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
# Yarn2dGrid new class to generate the domain by using overlapping method
#-------------------------------------------------------------------------

class Yarn2DOverlapping(object):
    NONTOUCH_FAC = 1.01
    
    def __init__(self, cfg):
        self.cfg = cfg
        
        self.Ry = self.cfg.get('domain.yarnradius')
        self.Rf = self.cfg.get('fiber.radius_fiber')
        self.scaleL = 1./self.Ry #get the scale factor for relative domain
        self.radius_yarn = self.scaleL * self.Ry
        self.radius_fiber = [self.scaleL * rad for rad in self.Rf]
        self.type_fiber = len(self.radius_fiber)
        self.blend = self.cfg.get('fiber.blend')
        assert len(self.blend) == len(self.radius_fiber)
        self.radius_boundlayer = max(self.radius_fiber)/2.
        self.radius_domain = self.radius_yarn + self.radius_boundlayer
        self.cellsize_centre = self.cfg.get('domain.cellsize_centre')
        self.cellSize = self.cfg.get('domain.cellsize_fiber')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.theta_value = self.cfg.get('domain.theta_value')
        self.beta_value = self.cfg.get('domain.beta_value')

        self.verbose = self.cfg.get('general.verbose')

    def determine_centers(self):
        #1. generate first the circle bands in which the virtual locations are 
        #   placed
        self.total_circles = sp.empty(self.type_fiber)
        self.total_circles_shift = sp.empty(self.type_fiber)
        self.number_fiber_type = sp.empty(len(self.blend))
        for ind in sp.arange(len(self.blend)):
            self.number_fiber_type[ind] = int(self.blend[ind] / 100. * self.number_fiber)
        
        #the number of circle bands for each kind of fiber
        print self.radius_fiber[:]
        for ind in sp.arange(len(self.radius_fiber)):
            self.total_circles[ind] = int((self.radius_yarn - 
                        self.radius_fiber[ind]*self.NONTOUCH_FAC) / \
                        (2 * self.radius_fiber[ind]*self.NONTOUCH_FAC)) + 1
            self.total_circles_shift[ind] = int((self.radius_yarn - 
                        self.radius_fiber[ind]*self.NONTOUCH_FAC/2) / \
                        (2 * self.radius_fiber[ind]*self.NONTOUCH_FAC)) + 1
        #the radius of midpoint of each circle bands with fibers
        print 'the number of circle layers', self.total_circles[:]
        self.radius_circle = [0] * self.type_fiber
        self.radius_circle_shift = [0] * self.type_fiber
        for ind in range(self.type_fiber):
            self.radius_circle[ind] = sp.empty(self.total_circles[ind], float)
            self.radius_circle_shift[ind] = sp.empty(self.total_circles_shift[ind], float)
            self.radius_circle[ind][0] = 0.
            self.radius_circle_shift[ind][0] = \
                    self.radius_fiber[ind]*self.NONTOUCH_FAC * (1/2+1/4)
            for i_circles in sp.arange(1, self.total_circles[ind]):
                self.radius_circle[ind][i_circles] =  \
                            self.radius_circle[ind][i_circles-1] +2. * \
                            (self.radius_fiber[ind]*self.NONTOUCH_FAC)
            for i_circles in sp.arange(self.total_circles_shift[ind]):
                self.radius_circle_shift[ind][i_circles] =  \
                            self.radius_circle_shift[ind][i_circles-1] + 2. * \
                            (self.radius_fiber[ind]*self.NONTOUCH_FAC)
        #2. generate the virtual locations
        #calculate the number of virtual location per band
        self.number_circle = [0] * self.type_fiber
        self.number_circle_shift = [0] * self.type_fiber
        self.total_number_virtloc = [0] * self.type_fiber
        self.total_number_virtloc_shift = [0] * self.type_fiber
        for ind in sp.arange(self.type_fiber):
            self.number_circle[ind] = sp.empty(self.total_circles[ind])
            self.number_circle_shift[ind] = sp.empty(self.total_circles_shift[ind])
            for i_circle in sp.arange(self.total_circles[ind]):
                if i_circle == 0:
                    self.number_circle[ind][i_circle] = 1
                else:
                    self.number_circle[ind][i_circle] = int(sp.pi / (self.radius_fiber[ind] \
                                                    / self.radius_circle[ind][i_circle]))
                    #print 'the value', self.radius_fiber[0] /self.radius_circle[i_circle]
                    print 'the number in each circle', self.number_circle[ind][i_circle]
            self.total_number_virtloc[ind] = sum(self.number_circle[ind][:])
        
        #the number of virtual locations for fiber_shift
            for i_circle in sp.arange(self.total_circles_shift[ind]):
                if i_circle == 0:
                    self.number_circle_shift[ind][i_circle] = 1
                else:
                    self.number_circle_shift[ind][i_circle] = int(sp.pi / (self.radius_fiber[ind] \
                                                    / self.radius_circle_shift[ind][i_circle]))
                    #print 'the value', self.radius_fiber[0] /self.radius_circle1[i_circle]
                    print 'the number in each circle', self.number_circle_shift[ind][i_circle]
            self.total_number_virtloc_shift[ind] = sum(self.number_circle_shift[ind][:])
        """
        for ind in range(self.type_fiber):
            for cnt, workarr in [(1, self.radius_circle[ind]), 
                            (2, self.radius_circle_shift[ind])]:
                if cnt == 1:
                    result = sp.empty(self.total_circles[ind], int)
                    radcirc = self.radius_circle[ind]
                    for index in sp.arange(self.total_circles[ind]):
                        result[index] = int(sp.pi * radcirc[ind] / self.radius_fiber[ind])
                else:
                    result = sp.empty(self.total_circles_shift[ind], int)
                    radcirc = self.radius_circle_shift[ind]
                    for index in sp.arange(self.total_circles[ind]):
                        result[index] = int(sp.pi * radcirc[ind] / self.radius_fiber[ind])
                #result[:] = int(sp.pi * radcirc[:] / self.radius_fiber[ind])
                print result
                
                if cnt == 1:
                    result[0] = 1
                    self.number_circle[ind] = result
                else:
                    self.number_circle_shift[ind] = result
            self.total_number_virtloc[ind] = sum(self.number_circle[ind])
            self.total_number_virtloc_shift[ind] = sum(self.number_circle_shift[ind])
            print 'the number_circle', self.number_circle[ind]
        """
        #3. calculate the postion of each virtual location for fiber
        self.x_position_vl = [0.] * self.type_fiber
        self.y_position_vl = [0.] * self.type_fiber 
        self.x_position_vl_shift = [0.] * self.type_fiber
        self.y_position_vl_shift = [0.] * self.type_fiber
        for ind in range(self.type_fiber):
            self.x_position_vl[ind] = sp.empty(self.total_number_virtloc[ind], float)
            self.y_position_vl[ind] = sp.empty(self.total_number_virtloc[ind], float)
            self.x_position_vl_shift[ind] = sp.empty(self.total_number_virtloc_shift[ind], float)
            self.y_position_vl_shift[ind] = sp.empty(self.total_number_virtloc_shift[ind], float)
            i_point = 0
            for i_layers in sp.arange(self.total_circles[ind]):
                if i_layers == 0:
                    self.x_position_vl[ind][i_point] = 0.
                    self.y_position_vl[ind][i_point] = 0.
                    i_point += 1
                else:
                    each_circle = self.number_circle[ind][i_layers]
                    for i_position in sp.arange(each_circle):
                        x_position = self.radius_circle[ind][i_layers] * sp.cos(2 * i_position * sp.arcsin(self.radius_fiber[ind] \
                                    / self.radius_circle[ind][i_layers]))
                        y_position = self.radius_circle[ind][i_layers] * sp.sin(2 * i_position * sp.arcsin(self.radius_fiber[ind] \
                                    /self.radius_circle[ind][i_layers]))
                        self.x_position_vl[ind][i_point] = x_position
                        self.y_position_vl[ind][i_point] = y_position
                        i_point += 1
            i_point_1 = 0
            for i_layers in sp.arange(self.total_circles_shift[ind]):
                if i_layers == 0:
                    self.x_position_vl_shift[ind][i_point_1] = 0.
                    self.y_position_vl_shift[ind][i_point_1] = 0.
                    i_point_1 += 1
                else:
                    each_circle = self.number_circle_shift[ind][i_layers]
                    print 'the layer number', i_layers
                    for i_position in sp.arange(each_circle):
                        x_position = self.radius_circle_shift[ind][i_layers] \
                                    * sp.cos(2 * i_position * sp.arcsin(self.radius_fiber[ind] \
                                    / self.radius_circle_shift[ind][i_layers])) \
                                    + self.radius_fiber[ind]*self.NONTOUCH_FAC * (1/2+1/4)
                        y_position = self.radius_circle_shift[ind][i_layers]* sp.sin(2 * i_position * sp.arcsin(self.radius_fiber[ind] \
                                    /self.radius_circle_shift[ind][i_layers]))
                        self.x_position_vl_shift[ind][i_point_1] = x_position
                        self.y_position_vl_shift[ind][i_point_1] = y_position
                        i_point_1 += 1
        print 'the central point of virtual location', self.x_position_vl[:]
        #calculate the probability value of fiber in each circle zone for fiber_1
        self.probability_value = [0] * self.type_fiber
        self.probability_value_shift = [0] * self.type_fiber
        for ind in range(self.type_fiber):
            self.probability_value[ind] = sp.empty(self.total_circles[ind], float)
            self.probability_value_shift[ind] = sp.empty(self.total_circles_shift[ind], float)
            for i_num_fiber in sp.arange(self.total_circles[ind]):
                self.probability_value[ind][i_num_fiber] = (1 - 2 * self.theta_value) \
                                                            * sp.power((sp.exp(1) - \
                                                            sp.exp(self.radius_circle[ind][i_num_fiber] \
                                                            / self.radius_yarn))\
                                                            /(sp.exp(1) - 1), self.beta_value) \
                                                            + self.theta_value
            print 'the value of probability', self.probability_value
            #calculate the probability value of fiber in each circle zone for fiber_shift
            for i_num_fiber in sp.arange(self.total_circles_shift[ind]):
                self.probability_value_shift[ind][i_num_fiber] = (1 - 2 * self.theta_value) * sp.power((sp.exp(1) - \
                                                    sp.exp(self.radius_circle_shift[ind][i_num_fiber] / self.radius_yarn))\
                                                   /(sp.exp(1) - 1), self.beta_value) + self.theta_value
            print 'the value of probability', self.probability_value_shift
        
        #calculate how many fibers can be in each circle zone for fiber_1
        self.each_circle_zone_num = [0] * self.type_fiber
        self.each_circle_zone_num_shift = [0] * self.type_fiber
        total_number_fiber = sp.empty(self.type_fiber)
        total_number_fiber_shift = sp.empty(self.type_fiber)
        for ind in range(self.type_fiber):
            self.each_circle_zone_num[ind] = sp.empty(self.total_circles[ind])
            self.each_circle_zone_num_shift[ind] = sp.empty(self.total_circles_shift[ind])
            for i_num_fiber in sp.arange(self.total_circles[ind]):
                self.each_circle_zone_num[ind][i_num_fiber] = int(round(self.number_circle[ind][i_num_fiber] * \
                                                    self.probability_value[ind][i_num_fiber]))
                print self.each_circle_zone_num[ind][i_num_fiber]
                print self.number_circle[ind][i_num_fiber]
            total_number_fiber[ind] = sum(self.each_circle_zone_num[ind][:])
            print 'total_number', total_number_fiber
            #calculate how many fibers can be in each circle zone for fiber_shift
            for i_num_fiber in sp.arange(self.total_circles_shift[ind]):
                self.each_circle_zone_num_shift[ind][i_num_fiber] = int(round(self.number_circle_shift[ind][i_num_fiber] * \
                                                    self.probability_value_shift[ind][i_num_fiber]))
            total_number_fiber_shift[ind] = sum(self.each_circle_zone_num_shift[ind][:])
            print 'total_number', total_number_fiber_shift
        """
        total_number_fiber = total_number_fiber_1 + total_number_fiber_2
        if total_number_fiber < self.number_fiber:
            print 'it is impossible to reach the required number of fiber'
            sys.exit(0)
        """
        #distribute the fibers in the virtual locations
        number_fiber_in_loop = self.number_fiber_type[:]
        self.circle_loop = sp.zeros(self.type_fiber)
        
        self.x_position = [0] * self.type_fiber
        self.y_position = [0] * self.type_fiber
        location_number = [0] * self.type_fiber
        #location_number_
        for ind in sp.arange(self.type_fiber):
            i_determine = 0
            print 'the number_fiber_in_loop', number_fiber_in_loop[ind]
            while number_fiber_in_loop[ind] > 0:
                number_fiber_in_loop[ind] = number_fiber_in_loop[ind] - \
                                            self.each_circle_zone_num[ind][i_determine]
                i_determine += 1
                self.circle_loop[ind] += 1
                print 'the length of the each zone', len(self.each_circle_zone_num[ind])
                print 'i_determine', i_determine
            print 'the number of loop for distributing fibers', self.circle_loop[ind]
        for ind in sp.arange(len(self.blend)):
            self.number_fiber_type[ind] = int(self.blend[ind] / 100. * self.number_fiber)
        print self.number_fiber_type
        number_fiber_in_loop = self.number_fiber_type[:]
        for ind in sp.arange(self.type_fiber):
            determine_generate = 0 #the number of generated fiber
            index_position = 0
            i_circle_number = 0
            
            print 'the value of number_fiber_type', self.number_fiber_type[ind]
            self.x_position[ind] = sp.empty(self.number_fiber_type[ind])
            self.y_position[ind] = sp.empty(self.number_fiber_type[ind])
            while i_circle_number < self.circle_loop[ind]:
                if i_circle_number == 0:
                    self.x_position[ind][determine_generate] = 0.
                    self.y_position[ind][determine_generate] = 0.
                    index_position += 1
                    i_circle_number += 1
                    determine_generate += 1
                    number_fiber_in_loop[ind] = number_fiber_in_loop[ind] - 1
                    
                elif i_circle_number == 1:
                    for i_index in sp.arange(self.number_circle[ind][i_circle_number]):
                        self.x_central = self.x_position_vl[ind][i_index + 1]
                        self.y_central = self.y_position_vl[ind][i_index + 1]
                        self.x_position[ind][determine_generate] = self.x_central
                        self.y_position[ind][determine_generate] = self.y_central
                        index_position += 1
                        determine_generate += 1
                    i_circle_number += 1
                    number_fiber_in_loop[ind] = number_fiber_in_loop[ind] - 6
                    
                elif i_circle_number > 1 and i_circle_number < self.circle_loop[ind] - 1:
                    location_number[ind] = sp.empty(self.each_circle_zone_num[ind][i_circle_number], int)
                    for i_index in sp.arange(self.each_circle_zone_num[ind][i_circle_number]):
                        if i_index == 0:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[ind][i_circle_number])
                            random_position = int(a_position)
                            location_number[ind][i_index] = random_position
                            a = random_position 
                            b = random_position
                            self.x_central = self.x_position_vl[ind][a]#index_position + random_position]
                            self.y_central = self.y_position_vl[ind][b]#index_position + random_position]
                            self.x_position[ind][determine_generate] = self.x_central
                            self.y_position[ind][determine_generate] = self.y_central
                            determine_generate += 1
                            #index += 5
                        else:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[ind][i_circle_number])
                            random_position = int(a_position)
                            determine_value = (random_position == location_number[ind])
                            while determine_value.any() == True:
                                a_position = np.random.uniform(index_position, 
                                            index_position + self.number_circle[ind][i_circle_number])
                                random_position = int(a_position)
                                determine_value = (random_position == location_number[ind])
                            else:
                                location_number[ind][i_index] = random_position#.append(random_position)
                                self.x_central = self.x_position_vl[ind][random_position]
                                self.y_central = self.y_position_vl[ind][random_position]
                                self.x_position[ind][determine_generate] = self.x_central
                                self.y_position[ind][determine_generate] = self.y_central
                                determine_generate += 1
                                
                    index_position += self.number_circle[ind][i_circle_number]
                    number_fiber_in_loop[ind] = number_fiber_in_loop[ind] - self.each_circle_zone_num[ind][i_circle_number]
                    i_circle_number += 1
                    
                elif i_circle_number == self.circle_loop[ind] - 1:
                    print 'the last layer', number_fiber_in_loop[ind]
                    location_number[ind] = sp.empty(number_fiber_in_loop[ind], int)
                    for i_index in sp.arange(number_fiber_in_loop[ind]):
                        if i_index == 0:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[ind][i_circle_number] )
                            random_position = int(a_position)
                            location_number[ind][i_index] = random_position#.append(random_position)
                            a = random_position
                            b = random_position
                            self.x_central = self.x_position_vl[ind][a]#index_position + random_position]
                            self.y_central = self.y_position_vl[ind][b]#index_position + random_position]
                            self.x_position[ind][determine_generate] = self.x_central
                            self.y_position[ind][determine_generate] = self.y_central
                            determine_generate += 1
                        else:
                            a_position = np.random.uniform(index_position, 
                                        index_position + self.number_circle[ind][i_circle_number])
                            random_position = int(a_position)
                            determine_value = (random_position == location_number[ind])
                            while determine_value.any() == True:
                                a_position = np.random.uniform(index_position, 
                                            index_position + self.number_circle[ind][i_circle_number])
                                random_position = int(a_position)
                                determine_value = (random_position == location_number[ind])
                            else:
                                location_number[ind][i_index] = random_position#.append(random_position)
                                self.x_central = self.x_position_vl[ind][random_position]
                                self.y_central = self.y_position_vl[ind][random_position]
                                self.x_position[ind][determine_generate] = self.x_central
                                self.y_position[ind][determine_generate] = self.y_central
                            determine_generate += 1
                    index_position += self.number_circle[ind][i_circle_number]
                    i_circle_number += 1
        print 'virtual location passed'
        print self.x_position
        
        #shift position part
        for ind in sp.arange(len(self.blend)):
            self.number_fiber_type[ind] = int(self.blend[ind] / 100. * self.number_fiber)
        #self.number_fiber_type[0] = 2
        number_fiber_in_loop_shift = self.number_fiber_type[:]
        print 'number of two fibers', self.number_fiber_type
        self.circle_loop_shift = sp.zeros(self.type_fiber)
        
        self.x_position_shift = [0] * self.type_fiber
        self.y_position_shift = [0] * self.type_fiber
        location_number_shift = [0] * self.type_fiber
        for ind in sp.arange(self.type_fiber):
            i_determine_shift = 0
            print 'the number_fiber_in_loop', number_fiber_in_loop_shift[ind]
            while number_fiber_in_loop_shift[ind] > 0:
                print 'the number of each circle', self.each_circle_zone_num_shift[ind]
                print 'the number of in loop', number_fiber_in_loop_shift[ind]
                number_fiber_in_loop_shift[ind] = number_fiber_in_loop_shift[ind] - \
                                            self.each_circle_zone_num_shift[ind][i_determine_shift]
                i_determine_shift += 1
                self.circle_loop_shift[ind] += 1
                print self.circle_loop_shift[ind]
                print 'the length of the each zone', len(self.each_circle_zone_num_shift[ind])
                print 'i_determine', i_determine
            print 'the number of loop for distributing fibers', self.circle_loop_shift[ind]
        for ind in sp.arange(len(self.blend)):
            self.number_fiber_type[ind] = int(self.blend[ind] / 100. * self.number_fiber)
        print self.number_fiber_type
        number_fiber_in_loop_shift = self.number_fiber_type[:]
        for ind in sp.arange(self.type_fiber):
            determine_generate_shift = 0 #the number of generated fiber
            index_position_shift = 0
            i_circle_number_shift = 0
            self.x_position_shift[ind] = sp.empty(self.number_fiber_type[ind], float)
            self.y_position_shift[ind] = sp.empty(self.number_fiber_type[ind], float)
            while i_circle_number_shift < self.circle_loop_shift[ind]:
                if i_circle_number_shift == 0:
                    self.x_position_shift[ind][determine_generate_shift] = self.radius_fiber[ind] *\
                                                                    self.NONTOUCH_FAC * (1/2+1/4)
                    self.y_position_shift[ind][determine_generate_shift] = 0.
                    print 'determine_generate value', determine_generate_shift
                    index_position_shift += 1
                    i_circle_number_shift += 1
                    determine_generate_shift += 1
                    number_fiber_in_loop_shift[ind] = number_fiber_in_loop_shift[ind] - 1
                    print 'finish the first point', number_fiber_in_loop_shift[ind]
                    print 'value of the self.circle_loop_shift', self.circle_loop_shift[ind]
                    
                elif i_circle_number_shift == 1:
                    print 'print the second layer', self.number_circle_shift[ind][i_circle_number_shift]
                    for i_index in sp.arange(self.number_circle_shift[ind][i_circle_number_shift]):
                        self.x_central = self.x_position_vl_shift[ind][i_index + 1]
                        self.y_central = self.y_position_vl_shift[ind][i_index + 1]
                        #print 'x_central for shift', x_central
                        self.x_position_shift[ind][determine_generate_shift] = self.x_central
                        self.y_position_shift[ind][determine_generate_shift] = self.y_central
                        index_position_shift += 1
                        determine_generate_shift += 1
                    i_circle_number_shift += 1
                    number_fiber_in_loop_shift[ind] = number_fiber_in_loop_shift[ind] - self.number_circle_shift[ind][i_circle_number_shift - 1]
                    
                elif i_circle_number_shift > 1 and i_circle_number_shift < self.circle_loop_shift[ind] - 1:
                    
                    location_number_shift[ind] = sp.empty(self.each_circle_zone_num_shift[ind][i_circle_number_shift], int)
                    for i_index in sp.arange(self.each_circle_zone_num_shift[ind][i_circle_number_shift]):
                        if i_index == 0:
                            a_position = np.random.uniform(index_position_shift, 
                                        index_position_shift + self.number_circle_shift[ind][i_circle_number_shift])
                            random_position = int(a_position)
                            location_number_shift[ind][i_index] = random_position
                            a = random_position 
                            b = random_position
                            self.x_central = self.x_position_vl_shift[ind][a]#index_position + random_position]
                            self.y_central = self.y_position_vl_shift[ind][b]#index_position + random_position]
                            self.x_position_shift[ind][determine_generate_shift] = self.x_central
                            self.y_position_shift[ind][determine_generate_shift] = self.y_central
                            determine_generate_shift += 1
                            print 'determine value when layer is larger than 2', determine_generate_shift
                            #index += 5
                        else:
                            a_position = np.random.uniform(index_position_shift, 
                                        index_position_shift + self.number_circle_shift[ind][i_circle_number_shift])
                            random_position = int(a_position)
                            determine_value = (random_position == location_number_shift[ind])
                            while determine_value.any() == True:
                                a_position = np.random.uniform(index_position_shift, 
                                            index_position_shift + self.number_circle_shift[ind][i_circle_number_shift])
                                random_position = int(a_position)
                                determine_value = (random_position == location_number_shift[ind])
                            else:
                                location_number_shift[ind][i_index] = random_position#.append(random_position)
                                self.x_central = self.x_position_vl_shift[ind][random_position]
                                self.y_central = self.y_position_vl_shift[ind][random_position]
                                self.x_position_shift[ind][determine_generate_shift] = self.x_central
                                self.y_position_shift[ind][determine_generate_shift] = self.y_central
                                determine_generate_shift += 1
                                print 'etermine value when the layer is larger 2_1', determine_generate_shift
                                
                    index_position_shift += self.number_circle_shift[ind][i_circle_number_shift]
                    number_fiber_in_loop_shift[ind] = number_fiber_in_loop_shift[ind] - self.each_circle_zone_num_shift[ind][i_circle_number_shift]
                    i_circle_number_shift += 1
                    print 'how many fibers are left', number_fiber_in_loop_shift[ind]
                    
                    
                elif i_circle_number_shift == self.circle_loop_shift[ind] - 1:
                    print 'the last layer', number_fiber_in_loop_shift[ind]
                    location_number_shift[ind] = sp.empty(number_fiber_in_loop_shift[ind], int)
                    for i_index in sp.arange(number_fiber_in_loop_shift[ind]):
                        if i_index == 0:
                            a_position = np.random.uniform(index_position_shift, 
                                        index_position_shift + self.number_circle_shift[ind][i_circle_number_shift] )
                            random_position = int(a_position)
                            location_number_shift[ind][i_index] = random_position#.append(random_position)
                            a = random_position
                            b = random_position
                            self.x_central = self.x_position_vl_shift[ind][a]#index_position + random_position]
                            self.y_central = self.y_position_vl_shift[ind][b]#index_position + random_position]
                            self.x_position_shift[ind][determine_generate_shift] = self.x_central
                            self.y_position_shift[ind][determine_generate_shift] = self.y_central
                            determine_generate_shift += 1
                        else:
                            a_position = np.random.uniform(index_position_shift, 
                                        index_position_shift + self.number_circle_shift[ind][i_circle_number_shift])
                            random_position = int(a_position)
                            determine_value = (random_position == location_number_shift[ind])
                            while determine_value.any() == True:
                                a_position = np.random.uniform(index_position_shift, 
                                            index_position_shift + self.number_circle_shift[ind][i_circle_number_shift])
                                random_position = int(a_position)
                                determine_value = (random_position == location_number_shift[ind])
                            else:
                                print 'the lengh of shift x position', len(self.x_position_shift[ind])
                                print 'determine generate shift', determine_generate_shift
                                location_number_shift[ind][i_index] = random_position#.append(random_position)
                                self.x_central = self.x_position_vl_shift[ind][random_position]
                                self.y_central = self.y_position_vl_shift[ind][random_position]
                                self.x_position_shift[ind][determine_generate_shift] = self.x_central
                                self.y_position_shift[ind][determine_generate_shift] = self.y_central
                            determine_generate_shift += 1
                    index_position_shift += self.number_circle_shift[ind][i_circle_number_shift]
                    i_circle_number_shift += 1
        print 'virtual location passed'
        print self.x_position_shift[0]
        
        #take 50% numbers from each other
        self.number_vl_overlap = sp.empty(len(self.number_fiber_type))
        for ind in sp.arange(len(self.blend)):
            self.number_fiber_type[ind] = int(self.blend[ind] / 100. * self.number_fiber)
        for ind in sp.arange(len(self.number_fiber_type)):
            self.number_vl_overlap[ind] = int(self.number_fiber_type[ind] / 2.)
        print 'overlap value', self.number_vl_overlap
        self.position_half = [0] * self.type_fiber
        self.position_half_shift = [0] * self.type_fiber
        for ind in range(self.type_fiber):
            self.position_half[ind] = sp.empty(self.number_vl_overlap[ind])
            self.position_half_shift[ind] = sp.empty(self.number_vl_overlap[ind])
            i_half = 0
            while i_half < self.number_vl_overlap[ind]:
                a_position = np.random.uniform(0., self.number_fiber_type[ind])
                position_random = int(a_position)
                if i_half == 0:
                    self.position_half[ind][i_half] = position_random
                else:
                    determine_value = (position_random == self.position_half[ind])
                    while determine_value.any() == True:
                        a_position = np.random.uniform(0., self.number_fiber_type[ind])
                        position_random = int(a_position)
                        determine_value = (position_random == self.position_half[ind])
                    else:
                        self.position_half[ind][i_half] = position_random
                        print 'self.position_random', position_random 
                i_half += 1
            
            i_half_1 = 0
            while i_half_1 < self.number_vl_overlap[ind]:
                b_position = np.random.uniform(0.,self.number_fiber_type[ind])
                position_random = int(b_position)
                if i_half_1 == 0:
                    self.position_half_shift[ind][i_half_1] = position_random
                else:
                    determine_value = (position_random == self.position_half_shift[ind])
                    while determine_value.any() == True:
                        b_position = np.random.uniform(0., self.number_fiber_type[ind])
                        position_random = int(b_position)
                        determine_value = (position_random == self.position_half_shift[ind])
                    else:
                        self.position_half_shift[ind][i_half_1] = position_random
                i_half_1 += 1
        print 'the random position value', self.position_half
        #draw the position of central points from combination
        x_position_random_1 = sp.empty(self.number_vl_overlap[0], float)
        y_position_random_1 = sp.empty(self.number_vl_overlap[0], float)
        x_position_random_1_shift = sp.empty(self.number_vl_overlap[0], float)
        y_position_random_1_shift = sp.empty(self.number_vl_overlap[0], float)
        for index in sp.arange(self.number_vl_overlap[0]):
            a_ind = self.position_half[0][index]
            b_ind = self.position_half_shift[0][index]
            print 'a_ind', a_ind
            print 'b_ind', b_ind
            x_position_random_1[index] = self.x_position[0][a_ind]
            y_position_random_1[index] = self.y_position[0][a_ind]
            x_position_random_1_shift[index] = self.x_position_shift[0][b_ind]
            y_position_random_1_shift[index] = self.y_position_shift[0][b_ind]
        print x_position_random_1_shift
        
        plt.plot(x_position_random_1, y_position_random_1, 'o',
            x_position_random_1_shift, y_position_random_1_shift, 's')
        plt.axis([-1,1, -1, 1])
        plt.show()
        
    def create_circle_domain_gmesh(self, filename='yarn_new_1.geo', 
                                   filename1 = 'fib_centers_x_new1.geo',
                                   filename2 = 'fib_centers_y_new1.geo',
                                   regenerate=True):
        filepath = utils.OUTPUTDIR + os.sep + filename
        filepath1 = utils.OUTPUTDIR + os.sep + filename1
        filepath2 = utils.OUTPUTDIR + os.sep + filename2
        
            