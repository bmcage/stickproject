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
Module with different methods for distributing the fibers in the yarn layout
1. integral from the area proabability distribution function;

2. distribute by using probability function for virtual locations
"""
from __future__ import division
import os.path
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
from matplotlib.patches import Circle, Wedge, Polygon, Ellipse
from matplotlib.collections import PatchCollection
import pylab
import matplotlib

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
from lib.utils.arraycompare import fullcompare_array, circledist
#from fipy import Gmsh2D
from fipy import *
from yarn.config import FIBERLAYOUTS
from virtlocgeom import *
from probability_area import *

def integration_layout(radius_fiber, radius_yarn, radius_each_circle, number_fiber_distribution, first_center,
                                        original_function_pro, number_circle_central, x_position_vl, y_position_vl):
#
#radius_fiber: float, the radius value of the fiber;
#radius_yarn: float, the radius of the yarn;
#radius_each_circle: array, the radius of each ring zone;
#number_fiber_distribution: the number of fibers disbributed in the domain;
#original_function_pro: function, the area probability distribution function;
#number_circle_central: array, the number of virtual locations in each ring zone;
#x_position_vl, y_position_vl: array, the virtual locations in the  coordinate
#
    delta_zone_r = 2. * radius_fiber / 1000.
    coefficient_integration = sp.around(radius_yarn / delta_zone_r) 
    delta_zone_r = radius_yarn / coefficient_integration

    rad_zones = sp.empty(len(radius_each_circle)+1,  float)
    rad_zones[1:]= radius_each_circle[:]
    rad_zones[0] = first_center
    rad_zones[-1] = 1.
    rad_intervals = rad_zones[1:] - rad_zones[:-1]
    each_num_pro = []
    each_num_integration = sp.zeros(len(rad_intervals), float)
    
    delta_r_i = 2. * radius_fiber / 1000.
    prev_int = 0.
    for i_circle in sp.arange(len(rad_intervals)):
        grid = sp.linspace(rad_zones[i_circle], rad_zones[i_circle+1], 
                    sp.around((rad_zones[i_circle+1]-rad_zones[i_circle]) 
                                            / delta_r_i))
        step = grid[1] - grid[0]
        disc_vals = original_function_pro(grid)

        each_num_pro.append(2 * step * sp.sum(disc_vals * grid))
        each_num_integration[i_circle] = prev_int + (each_num_pro[-1] / 
                                                   (radius_fiber ** 2.))
        prev_int = each_num_integration[i_circle]
    print 'the value from the integration', each_num_integration
    #raw_input("check whether the value is equal to the input")
    #each_num_integration = each_num_integration
    if abs(each_num_integration[-1] - 1) > 0.03:
        print 'the precision of the integration cannot reach'
        print each_num_integration[-1]
        assert False
    each_num_integration = sp.around(each_num_integration * 
                                        number_fiber_distribution)
    print 'The number of fiber up to zone', each_num_integration
    print 'The number per zone', each_num_integration[0], each_num_integration[1:] \
                                - each_num_integration[:-1]

    each_num_integration[-1] = number_fiber_distribution
    each_circle_zone_num = sp.ones(len(each_num_integration))
    each_circle_zone_num[0] = each_num_integration[0]
    each_circle_zone_num[1:] = each_num_integration[1:] - each_num_integration[:-1]
    print each_circle_zone_num
        
    circle_loop = len(each_circle_zone_num)
    diff_num = 0
    less_num = 0
    position_circle = []
    for i_circle in sp.arange(len(each_circle_zone_num)):
        if each_circle_zone_num[i_circle] > number_circle_central[i_circle]:
            diff_num += (each_circle_zone_num[i_circle] - number_circle_central[i_circle])
            each_circle_zone_num[i_circle] = number_circle_central[i_circle]
        elif each_circle_zone_num[i_circle] < number_circle_central[i_circle]:
            less_num += (number_circle_central[i_circle] - each_circle_zone_num[i_circle])
            position_circle.append(i_circle)
    for i_add in sp.arange(len(position_circle)):
        add_position = position_circle[i_add]
        difference = number_circle_central[add_position] - each_circle_zone_num[add_position]
        if difference <= diff_num:
            each_circle_zone_num[add_position] = number_circle_central[add_position]
            diff_num -= difference
        else:
            each_circle_zone_num[add_position] += diff_num 
            diff_num -= diff_num
    print 'afer rearranging the number of fibers', each_circle_zone_num
    filename = utils.OUTPUTDIR + os.sep + "virlayout.gz"    
    x_position = sp.empty(number_fiber_distribution)
    y_position = sp.empty(number_fiber_distribution)
    determine_generate = 0 #the number of generated fiber
    index_position = 0
    i_circle_number = 0
    number_fiber_in_loop = number_fiber_distribution
    #print 'begin to distribute the position', circle_loop
    while i_circle_number < circle_loop: 
        #print i_circle_number
        location_number = sp.zeros(each_circle_zone_num[i_circle_number]) - 1
        for i_index in sp.arange(each_circle_zone_num[i_circle_number]):
            a_position = np.random.uniform(index_position, 
                        index_position + number_circle_central[i_circle_number])
            #print a_position
            random_position = int(a_position)
            determine_value = (random_position == location_number)
            while determine_value.any() == True:
                a_position = np.random.uniform(index_position, 
                            index_position + number_circle_central[i_circle_number])
                random_position = int(a_position)
                determine_value = (random_position == location_number)
            else:
                x_position[determine_generate] = x_position_vl[random_position]
                y_position[determine_generate] = y_position_vl[random_position]
                location_number[i_index] = random_position
                determine_generate += 1
        index_position += number_circle_central[i_circle_number]
        number_fiber_in_loop = number_fiber_in_loop \
                                    - each_circle_zone_num[i_circle_number]
        i_circle_number += 1
    return (each_circle_zone_num, x_position, y_position)

def vl_distribution(num_fiber_VL, number_fiber_distribution, x_position_VL, y_position_VL):
#
#number_fiber_VL: array, the virtual locations number in each ring zone;
#number_fiber_distribution: integer, the number of the fiber distributed;
#x_position_VL, y_position_VL: array, the position value of virtual locations in the coordinates
#
    each_circle_zone_num = sp.zeros(len(num_fiber_VL))
    total_sum_fiber_VL = sp.sum(num_fiber_VL)
    if total_sum_fiber_VL < number_fiber_distribution:
        print "the number of virtual location for fibers should be larger than that \
            of .ini file"
        assert False
    x_position = sp.empty(number_fiber_distribution)
    y_position = sp.empty(number_fiber_distribution)
    determine_generate = 0 #the number of generated fiber
    index_position = 0
    i_circle_number = 0
    index_circle_num = 0
    circle_loop = 0
    number_fiber_in_loop = number_fiber_distribution
    while number_fiber_in_loop > 0:
        number_fiber_in_loop -= num_fiber_VL[index_circle_num]
        each_circle_zone_num[index_circle_num] = num_fiber_VL[index_circle_num]
        index_circle_num += 1
        circle_loop += 1
    each_circle_zone_num[index_circle_num - 1] += number_fiber_in_loop 
    number_fiber_in_loop = number_fiber_distribution
    while i_circle_number < circle_loop:
        if i_circle_number < circle_loop - 1:
            location_number = sp.zeros(num_fiber_VL[i_circle_number]) - 1
            for i_index in sp.arange(num_fiber_VL[i_circle_number]):
                a_position = np.random.uniform(index_position, index_position + 
                            num_fiber_VL[i_circle_number])
                random_position = int(a_position)
                determine_value = (random_position == location_number)
                while determine_value.any() == True:
                    a_position = np.random.uniform(index_position, 
                                index_position + num_fiber_VL[i_circle_number])
                    random_position = int(a_position)
                    determine_value = (random_position == location_number)
                else:
                    x_position[determine_generate] = x_position_VL[random_position]
                    y_position[determine_generate] = y_position_VL[random_position]
                    location_number[i_index] = random_position
                    determine_generate += 1
            index_position += num_fiber_VL[i_circle_number]
            number_fiber_in_loop = number_fiber_in_loop \
                                        - num_fiber_VL[i_circle_number]
        else:
            location_number = sp.zeros(number_fiber_in_loop)
            new_index = number_fiber_in_loop
            for i_index in sp.arange(new_index):
                a_position = np.random.uniform(index_position, index_position + 
                            num_fiber_VL[i_circle_number])
                random_position = int(a_position)
                determine_value = (random_position == location_number)
                while determine_value.any() == True:
                    a_position = np.random.uniform(index_position, 
                                index_position + num_fiber_VL[i_circle_number])
                    random_position = int(a_position)
                    determine_value = (random_position == location_number)
                else:
                    x_position[determine_generate] = x_position_VL[random_position]
                    y_position[determine_generate] = y_position_VL[random_position]
                    location_number[i_index] = random_position
                    determine_generate += 1
            index_position += num_fiber_VL[i_circle_number]
        i_circle_number += 1
    return (each_circle_zone_num, x_position, y_position)

def compare_relative_error(mean_value_each, mean_value_alpha, nrzones):
    ind = sp.arange(nrzones)
    width = 0.2
    
    plt.figure()
    propor_draw = plt.bar(ind, mean_value_each, width, color = 'b')
    alpha_draw = plt.bar(ind + width, mean_value_alpha, color = 'y')
    
    plt.ylabel(ur'$\Delta D = \frac{\sum\limits_{i=1}^{n}(p_{f}^{i} - p_{a}^{i})^2}{p_{f}^{i}}$')
    plt.xticks(ind + width, ('Zone1', 'Zone2', 'Zone3', 'Zone4', 'Zone4'))
    plt.legend((propor_draw[0], alpha_draw[0]), ('$\alpha = \frac{R_{f}^{m}}{R_{f}^{n}}$', '$\alpha = \frac{1}{2}$'))
    
    plt.show()
