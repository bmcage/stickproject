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
Module with precomputed functions needed for virt loc computation:
1.regular virt loc
2.shifted virt loc
3.combined non-overlapping virt loc
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
import lib.utils.utils as utils
from lib.utils.arraycompare import fullcompare_array, circledist
from fipy import Gmsh2D
from fipy import *
from yarn.config import FIBERLAYOUTS
from virtlocgeom import *
from arearatioprobability import *

from distribution_method import *

#-------------------------------------------------------------------------
#
# three types of layout
#
#-------------------------------------------------------------------------

NONTOUCH_FAC = 1.1
def randomfiberlayout(options):
    """ Generate the fiber layout in the yarn in random fashion
    Options should be a dictionary with keys:
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
    x_central = options.get('x_central', 0.)
    y_central = options.get('y_central', 0.)
    onumber_fiber = options.get('number_fiber', 1)
    onumber_fiber_blend = options.get('number_fiber_blend', [1])
    oradius_fiber = options.get('radius_fiber', [1.])
    oradius_yarn = options.get('radius_yarn', 2.)

    ##if len(options.radius_fiber) > 1:
    ##    print 'ERROR: randomfiberlayout can only handle one type of fiber'
    ##    sys.exit()

    #allocate space
    x_position = sp.empty(onumber_fiber, float)
    y_position = sp.empty(onumber_fiber, float)
    radius_fiber = sp.empty(onumber_fiber, float)
    fiber_kind = sp.empty(onumber_fiber, int)
    
    #preset the radius_fiber. This algorithm does first the largest fibers,
    #then the others
    ind = sp.arange(onumber_fiber)
    radius_fiber[:onumber_fiber_blend[0]] = oradius_fiber[0]
    fiber_kind[:onumber_fiber_blend[0]] = 0
    for j in range(1, len(onumber_fiber_blend)):
        assert oradius_fiber[j] <= oradius_fiber[j-1], \
                'Give fibers in decreasing size in ini file'
        radius_fiber[sum(onumber_fiber_blend[:j]):sum(onumber_fiber_blend[:j+1])]\
            = oradius_fiber[j]
        fiber_kind[onumber_fiber_blend[j-1]:onumber_fiber_blend[j] + onumber_fiber_blend[j-1]] = j
    print 'Putting', onumber_fiber, 'fibers'
    for i in sp.arange(0, onumber_fiber, 1):
        
        #generate the position of fiber
        a = np.random.uniform(-1, 1)
        b = np.random.uniform(-1, 1)
        #distance between the current point and center
        distance_center = sp.sqrt((a - x_central)**2 + (b - y_central)**2)
        #distance between the current point and existing points
        distance_each = sp.sqrt((a - x_position[:i])**2 + \
                                (b - y_position[:i])**2)
        while (distance_center + radius_fiber[i] >= oradius_yarn or
               (i>0 and not (( distance_each
                >= (1.001*radius_fiber[i] + radius_fiber[:i])).all()))): 
            #   (i>0 and (np.min(distance_each)\
            #    <= (2*oradius_fiber[0])+ 0.001 * oradius_fiber[0]))):
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
        print i
    plot_yarn(x_position, y_position, radius_fiber)
    return (x_position, y_position, radius_fiber, fiber_kind)

def virtloclayout(options):
    """ Generate the fiber layout in the yarn in virtual locations
    Options should contain attributes:
        number_fiber : amount of fibers to generate
        radius_yarn
        radius_fiber: array of length 1 with radius of the fiber
        radius_first_center: value where first virtloc center is, default=0.
    
    Returns a tuple 
        (list of x coord center points,
         list of y coord center points,
         list of radius of fiber,
         list integers indicating kind of fiber
         list of the number of fibers in each ring zones
         )
    """
    x_central = options.get('x_central', 0.)
    y_central = options.get('y_central', 0.)
    onumber_fiber = options.get('number_fiber', 1)
    onumber_fiber_blend = options.get('number_fiber_blend', [1])
    oradius_fiber = options.get('radius_fiber', [1.])
    oradius_yarn = options.get('radius_yarn', 2.)
    ofirst_center = options.get('radius_first_center', 0.0)
    
    otheta_value = options.get('theta_value', 0.1)
    obeta_value = options.get('beta_value', 0.05)
    
    omean_deviation = options.get('mean_deviation', 0.0)
    ##The method of the integration
    original_function_pro = options.get('prob_area',  lambda r: r**1 )
    #original_function_pro = original_function_pro[0]
    distribute_fiber = options.get('layout_distribution', 'virtual_location')
    
    filename_1 = utils.OUTPUTDIR + os.sep + "proportion_vl_value.gz"
    if len(oradius_fiber) > 1 or len(onumber_fiber_blend) > 1:
        print 'ERROR: for virtual location layout the number of fibers must be 1.'
        print 'Actual number of type of fiber %d' %len(oradius)
        assert False
    if onumber_fiber != onumber_fiber_blend[0]:
        print 'ERROR: number fiber and blend do not correspond'
        assert False
    print "two values:", oradius_fiber[0], omean_deviation
    oradius_fiber = oradius_fiber[0] + omean_deviation
    
    """
     the first"1" means the central part in the yarn; the second "1" means the 
     zone on the boundary yarn
    """
    if ofirst_center != 0. :
        total_circles = int((oradius_yarn - (ofirst_center + 2. * 
                                    oradius_fiber*NONTOUCH_FAC))
                                    / (2 * oradius_fiber*NONTOUCH_FAC)) + 1 
    else:
        total_circles = int((oradius_yarn - oradius_fiber)
                        / (2. * oradius_fiber)) + 1 
    
    number_circle_central = sp.empty(total_circles, int)
    radius_circle_central = sp.empty(total_circles, float)
    radius_each_circle = sp.empty(total_circles, float)
    radius_fiber = sp.ones(onumber_fiber_blend[0], float)*(oradius_fiber - omean_deviation)
    fiber_kind = sp.zeros(onumber_fiber_blend[0], int)
    ind = sp.arange(onumber_fiber)
    plus_angle = sp.zeros(total_circles, float)
    if ofirst_center == 0:
        for i_circle in sp.arange(total_circles):
            radius_circle_central[i_circle] =  (i_circle * oradius_fiber * 2
                                                 + ofirst_center)
            radius_each_circle[i_circle] = radius_circle_central[i_circle] + oradius_fiber * NONTOUCH_FAC
            number_circle_central[i_circle]= max([int(sp.pi
                                            * radius_circle_central[i_circle]
                                            / (oradius_fiber)), 1])
            plus_angle[i_circle] = 2. * sp.pi / number_circle_central[i_circle]
            print 'the number of virtual locations in each zone', number_circle_central[i_circle]
    else:
        for i_circle in sp.arange(total_circles):
            radius_circle_central[i_circle] = ((2 * i_circle + 1) * oradius_fiber 
                                                * NONTOUCH_FAC + ofirst_center)
            radius_each_circle[i_circle] = radius_circle_central[i_circle] + oradius_fiber * NONTOUCH_FAC
            number_circle_central[i_circle]= max([int(sp.pi
                                            * radius_circle_central[i_circle]
                                            / (oradius_fiber)), 1])
            plus_angle[i_circle] = 2. * sp.pi / number_circle_central[i_circle]
            print 'the number of virtual locations in each zone', number_circle_central[i_circle]
    total_number_vl = sum(number_circle_central[:])
    ##area of each ring zone
    radius_vl = sp.ones(total_number_vl)
    radius_vl[:] = oradius_fiber
    area_ring_zone = sp.empty(len(radius_circle_central), float)
    
    pro_fiber_VL = sp.empty(len(radius_circle_central), float)
    num_fiber_VL = sp.empty(len(radius_circle_central), float)
    
    if ofirst_center == 0:
        pre_area = 0.
        for i_ring_zone in sp.arange(len(radius_circle_central)):
            area_ring_zone[i_ring_zone] = sp.pi * sp.power(radius_each_circle[i_ring_zone],
                                        2.) - pre_area
            pre_area = sp. pi * sp.power(radius_each_circle[i_ring_zone], 2.)
            
            #just for the VL-fiber distribution
            if distribute_fiber == 'virtual_location':
                pro_fiber_VL[i_ring_zone] = prob_func_VL(radius_circle_central[i_ring_zone],
                                        oradius_yarn, otheta_value, obeta_value)
                num_fiber_VL[i_ring_zone] = max(sp.around(number_circle_central[i_ring_zone] * \
                                            pro_fiber_VL[i_ring_zone]), 1)
    else:
        area_shift = sp.pi * sp.power(ofirst_center, 2.)
        pre_area = area_shift
        for i_ring_zone in sp.arange(len(radius_circle_central)):
            area_ring_zone[i_ring_zone] = sp.pi * sp.power(radius_each_circle[i_ring_zone], 
                                        2.) - pre_area
            pre_area = sp.pi * sp.power(radius_each_circle[i_ring_zone],  2.)
            #just for the VL-fiber distribution
            if distribute_fiber == 'virtual_location':
                pro_fiber_VL[i_ring_zone] = prob_func_VL(radius_circle_central[i_ring_zone],
                                        oradius_yarn, otheta_value, obeta_value)
                num_fiber_VL[i_ring_zone] = sp.around(number_circle_central[i_ring_zone] * 
                                        pro_fiber_VL[i_ring_zone])
    #print 'sum of the locations from VL', sp.sum(num_fiber_VL)
    #raw_input("enter")
##    print 'the total number postions for the fiber', sp.sum(num_fiber_VL), onumber_fiber_blend[0], \
##            num_fiber_VL
##    print 'the number fiber in each ring zone', num_fiber_VL
    if onumber_fiber > total_number_vl:
        print 'ERROR: the number of fiber is more than the virtual locations'
        print 'the total number of virtual locations', total_number_vl
        sys.exit(0)    
        
    #each vitrual location position
    x_position_vl = []
    y_position_vl = []
    for i_circle in sp.arange(total_circles):
        each_circle = number_circle_central[i_circle]
        if each_circle == 1:
            x_position_vl.append(0.)
            y_position_vl.append(0.)
        else:    
            for i_position in sp.arange(each_circle):
                x_position_t = radius_circle_central[i_circle] * sp.cos(i_position * 
                            plus_angle[i_circle])
                y_position_t = radius_circle_central[i_circle] * sp.sin(i_position *
                            plus_angle[i_circle])
                x_position_vl.append(x_position_t)
                y_position_vl.append(y_position_t)

    number_fiber_each_zone = sp.empty(total_circles, integer)
    noapp_number_fiber = sp.empty(total_circles, float)
    delta_number = sp.zeros(total_circles - 1, float)

    if hasattr(options, 'x_central'):
        x_central = options.x_central
    else:
        x_central = 0.
    if hasattr(options, 'y_central'):
        y_central = options.y_central
    else:
        y_central = 0.
    circle_loop = 0
    i_determine = 0
    number_fiber_in_loop = onumber_fiber_blend[0]
    i_type = len(onumber_fiber_blend)
    print "choosing distributing the fibers", distribute_fiber
    #use the integral of area probability function to distribute the fibers
    if distribute_fiber == 'integral':
        each_circle_zone_num, x_position, y_position = integration_layout(oradius_fiber, 
                                        oradius_yarn, radius_each_circle, onumber_fiber_blend[0], 
                                        ofirst_center, original_function_pro, number_circle_central,
                                        x_position_vl, y_position_vl)
    elif distribute_fiber == 'virtual_location':
   #distribute the fiber into the VL from VL-fiber function
        print 'begin with distributing by VLM'
        each_circle_zone_num, x_position, y_position = vl_distribution(num_fiber_VL, onumber_fiber_blend[0], 
                                                                                x_position_vl, y_position_vl)

    #Calculate the ratio value in each ring zone
        zone_position, ratio_value = calculate_proportion(oradius_yarn, radius_fiber, 
                                    x_position, y_position)
        print 'the ratio value for the single kind of fiber', ratio_value
    #polynomial fitting
        zone_position[:] = zone_position[:] - zone_position[0]
        zone_position = np.append(zone_position, [1.0])
        ratio_value = np.append(ratio_value, [0.0])
        poly_single_kind = np.poly1d(np.polyfit(zone_position, ratio_value, 5))
        print 'the coefficient for the polynomial equations', np.polyfit(zone_position, ratio_value, 5)
##    zone_p = sp.linspace(0., 1., 50)
##    plt.figure()
##    plt.plot(zone_p, poly_single_kind(zone_p), '--')
##    plt.plot(zone_position, ratio_value, 's')
##    plt.xlim(0., 1.05)
##    plt.ylim(0., 0.8)
##    plt.xlabel('Relative position in the yarn domain')
##    plt.ylabel('Probability value')
##    plt.draw()
##    plt.show()
##    raw_input("enter")
#   output the data on fibers position
    oradius_fiber_array = sp.zeros(onumber_fiber_blend[0], float)
    oradius_fiber_array[:] = oradius_fiber - omean_deviation
#    dump.write({'x_position':x_position, 'y_position':y_position, 'radius_fiber':oradius_fiber,
#                },
#                filename = filename, extension = '.gz')
    #fig = plot_yarn(x_position, y_position, oradius_fiber_array)#, fiber_kind)
    return (x_position, y_position, radius_fiber, fiber_kind, x_position_vl, 
            y_position_vl,  each_circle_zone_num,  radius_circle_central,  area_ring_zone)

def virtlocoverlaplayout(options):
    """ Generate the fiber layout in the yarn in virtual locations with random distribution
    Options should contain attributes:
        number_fiber : amount of fibers to generate
        radius_yarn
        radius_fiber: array of length 1 with radius of the fiber
    
    Returns a tuple 
        (list of x coord center points,
         list of y coord center points,l
         list of radius of fiber,
         list integers indicating kind of fiber
        )
    """
    x_central = options.get('x_central', 0.)
    y_central = options.get('y_central', 0.)
    onumber_fiber = options.get('number_fiber', 1)
    onumber_fiber_blend = options.get('number_fiber_blend', [1])
    oradius_fiber = options.get('radius_fiber', [1.])

    oradius_yarn = options.get('radius_yarn', 2.)
    
    otheta_value = options.get('theta_value', 0.1)
    obeta_value = options.get('beta_value', 0.05)
    
    omean_deviation = options.get('mean_deviation', [0.0])
    oprob_area = options.get('prob_area',  [ lambda r: r ** 1 ] )
    
    distribute_fiber = options.get('layout_distribution', 'integral')

    onumber_fiber_each = sp.zeros(len(onumber_fiber_blend), float)
    onumber_fiber_each[:] = onumber_fiber_blend[:]
    type_fiber = len(onumber_fiber_blend)
    print 'generate fiber_kind value',type_fiber
    raw_input("check the fiber kind number")
    x_position = [0] * type_fiber
    y_position = [0] * type_fiber
    each_num_circle = [0] * type_fiber
    radius_fiber = [0] * type_fiber
    fiber_kind = [0] * type_fiber

    x_position_shift = [0] * type_fiber
    y_position_shift = [0] * type_fiber
    each_num_circle_shift = [0] * type_fiber
    radius_fiber_shift = [0] * type_fiber
    fiber_kind_shift = [0] * type_fiber
    
    position_half = [0] * type_fiber
    position_half_shift = [0] * type_fiber
    
    number_vl_overlap = sp.empty(len(onumber_fiber_blend))
    for ind in sp.arange(len(onumber_fiber_blend)):
        number_vl_overlap[ind] = int(onumber_fiber_blend[ind] / 2.)
        half_for_each = number_vl_overlap[ind]
    for i_type in sp.arange(type_fiber):
        position_half[i_type] = sp.empty(onumber_fiber_blend[i_type], int)
        x_position[i_type] = sp.zeros(onumber_fiber_blend[i_type], float)
        y_position[i_type] = sp.zeros(onumber_fiber_blend[i_type], float)
        
        radius_fiber[i_type] = sp.zeros(onumber_fiber_blend[i_type],float)
        fiber_kind[i_type] = sp.empty(onumber_fiber_blend[i_type], int)
        filename = utils.OUTPUTDIR + os.sep + "virlayout_%d.gz"%(i_type) 
        filename_1 = utils.OUTPUTDIR + os.sep + "virlayout_area%d.gz"%(i_type) 
        ouroptions = {
                'x_central' : x_central,
                'y_central' : y_central,
                'number_fiber' : onumber_fiber_blend[i_type],
                'number_fiber_blend' : [onumber_fiber_blend[i_type]],
                'radius_fiber' : [oradius_fiber[i_type]],
                'radius_yarn' : oradius_yarn,
                'theta_value' : otheta_value,
                'beta_value' : obeta_value,
                'mean_deviation': omean_deviation[i_type],
                'prob_area': oprob_area[i_type],
                'layout_distribution': distribute_fiber
                }
        ouroptions['radius_first_center'] = 0.             
        ax_position, ay_position, aradius_fiber, afiber_kind, x_vl, y_vl, \
        each_circle_zone_num, radius_circle_central, area_ring_zone = virtloclayout(ouroptions)
        proportion_no_sub = sp.pi * sp.power(oradius_fiber[0], 2.) * each_circle_zone_num / \
        area_ring_zone
        x_position[i_type][:] = ax_position[:]
        y_position[i_type][:] = ay_position[:]
        radius_fiber[i_type] = aradius_fiber[:]
        zone_position_one_kind, ratio_one_kind = calculate_proportion(oradius_yarn,
                                            aradius_fiber, ax_position, ay_position)
##        for i_rad in sp.arange(len(aradius_fiber)):
##            radius_fiber[i_type][i_rad] = np.random.normal(aradius_fiber[0], 
##                                omean_deviation[i_type])  
        fiber_kind[i_type][:] = afiber_kind[:]
        each_num_circle[i_type] = sp.empty(len(area_ring_zone), int)
        each_num_circle[i_type][:] = each_circle_zone_num[:]

        for i_point in sp.arange(len(x_position[i_type])):
            position_half[i_type][i_point] = int(i_point)
            
        dump.write({'x_position':ax_position, 'y_position':ay_position, 'radius_fiber':aradius_fiber,
                    },
                    filename = filename, extension = '.gz')
        
    #for the shift part
    number_vl_overlap_shift = onumber_fiber_blend - number_vl_overlap
    for i_type in sp.arange(type_fiber):
        x_position_shift[i_type] = sp.zeros(number_vl_overlap_shift[i_type], float)
        y_position_shift[i_type] = sp.zeros(number_vl_overlap_shift[i_type], float)
        radius_fiber_shift[i_type] = sp.zeros(number_vl_overlap_shift[i_type], float)
        fiber_kind_shift[i_type] = sp.empty(number_vl_overlap_shift[i_type])

        x_position_shift[i_type] = sp.zeros(onumber_fiber_blend[i_type], float)
        y_position_shift[i_type] = sp.zeros(onumber_fiber_blend[i_type], float)
        radius_fiber_shift[i_type] = sp.zeros(onumber_fiber_blend[i_type], float)
        fiber_kind_shift[i_type] = sp.empty(onumber_fiber_blend[i_type])
        filename = utils.OUTPUTDIR + os.sep + "virlayshift_%d.gz"%(i_type)
        filename_1 = utils.OUTPUTDIR + os.sep + "virlayshift_area%d.gz"%(i_type)
        ouroptions = {
                'x_central' : x_central,
                'y_central' : y_central,
                'number_fiber' : onumber_fiber_each[i_type],
                'number_fiber_blend' : [onumber_fiber_each[i_type]],
                'radius_fiber' : [oradius_fiber[i_type]],
                'radius_yarn' : oradius_yarn,
                'theta_value' : otheta_value,
                'beta_value' : obeta_value,
                'mean_deviation': omean_deviation[i_type],
                'prob_area': oprob_area[i_type],
                'layout_distribution': distribute_fiber
        }
        ouroptions['radius_first_center'] = 0.5 * oradius_fiber[i_type]
        ax_position_shift, ay_position_shift, aradius_fiber_shift, afiber_kind_shift,\
        x_vl_shift,y_vl_shift,  each_circle_zone_num_shift,  \
        radius_circle_central_shift, area_ring_zone_shift= virtloclayout(ouroptions)
        x_position_shift[i_type][:] = ax_position_shift[:]
        y_position_shift[i_type][:] = ay_position_shift[:]
        radius_fiber_shift[i_type][:] = aradius_fiber_shift[:]

##        for i_rad in sp.arange(len(aradius_fiber_shift)):
##            radius_fiber_shift[i_type][i_rad] = np.random.normal(aradius_fiber_shift[0],
##                                            omean_deviation[i_type])
        
        fiber_kind_shift[i_type][:] = afiber_kind_shift[:]
        each_num_circle_shift[i_type] = sp.empty(len(area_ring_zone_shift),int)
        each_num_circle_shift[i_type][:] = each_circle_zone_num_shift[:]
        zone_position_shifted, ratio_shifted = calculate_proportion(oradius_yarn,
                                            radius_fiber_shift[i_type], ax_position_shift,
                                            ay_position_shift)
        
        dump.write({'x_position':ax_position_shift, 'y_position':ay_position_shift, 
                    'radius_fiber':aradius_fiber_shift},
                    filename = filename, extension = '.gz')
    #take 50% points from each other
    position_half = [0] * type_fiber
    position_half_shift = [0] * type_fiber
    for ind in range(type_fiber):
        i_half = 0
        #no shift part
        position_each = []
        sum_no_shift = 0
        begin_point = 0
        end_point = 0
        for i_circle in sp.arange(len(each_num_circle[ind])):
            number_chosen = sp.around(each_num_circle[ind][i_circle] / 2. + 0.001)
            sum_no_shift += number_chosen
            i_half_each = 0
            end_point += each_num_circle[ind][i_circle]
            while i_half_each < number_chosen:
                a_position = np.random.uniform(begin_point, end_point)
                position_random = int(a_position)
                determine_value = (position_random == sp.array(position_each))
                while determine_value.any() == True:
                    a_position = np.random.uniform(begin_point, end_point)
                    position_random = int(a_position)
                    determine_value = (position_random == sp.array(position_each))
                else:
                    position_each.append(position_random)
                i_half_each += 1
            begin_point += each_num_circle[ind][i_circle]
        position_half[ind] = sp.empty(sum_no_shift, int)
        for i_position in sp.arange(len(position_half[ind])):
            position_half[ind][i_position] = position_each[i_position]
        left_num_shift = onumber_fiber_blend[ind] - sum_no_shift   
        #shift part
        sum_with_shift = 0
        begin_point_shift = 0
        end_point_shift = 0
        position_each_shift = []
        #left_num = left_num_shift
        if ind == 0:
            for i_circle in sp.arange(1, len(each_num_circle_shift[ind]), 1):
                number_chosen = sp.around(each_num_circle_shift[ind][i_circle] / 2. + 0.001)
                i_half_each = 0
                end_point_shift += each_num_circle_shift[ind][i_circle]
                left_num = left_num_shift
                left_num -= number_chosen
                if left_num >= 0:
                    while i_half_each < number_chosen:
                        a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                        position_random = int(a_position)
                        determine_value = (position_random == sp.array(position_each_shift))
                        while determine_value.any() == True:
                            a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                            position_random = int(a_position)
                            determine_value = (position_random == sp.array(position_each_shift))
                        else:
                            position_each_shift.append(position_random)
                        i_half_each += 1
                    begin_point_shift += each_num_circle_shift[ind][i_circle]
                
                else:
                    while i_half_each < (left_num + number_chosen):
                        a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                        position_random = int(a_position)
                        determine_value = (position_random == sp.array(position_each_shift))
                        while determine_value.any() == True:
                            a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                            position_random = int(a_position)
                            determine_value = (position_random == position_each_shift)
                        else:
                            position_each_shift.append(position_random)
                        i_half_each += 1
        else:
            for i_circle in sp.arange(len(each_num_circle_shift[ind])):
                number_chosen = sp.around(each_num_circle_shift[ind][i_circle] / 2. + 0.001)
                i_half_each = 0
                end_point_shift += each_num_circle_shift[ind][i_circle]
                left_num = left_num_shift
                left_num -= number_chosen
                if left_num >= 0:
                    while i_half_each < number_chosen:
                        a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                        position_random = int(a_position)
                        determine_value = (position_random == sp.array(position_each_shift))
                        while determine_value.any() == True:
                            a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                            position_random = int(a_position)
                            determine_value = (position_random == sp.array(position_each_shift))
                        else:
                            position_each_shift.append(position_random)
                        i_half_each += 1
                    begin_point_shift += each_num_circle_shift[ind][i_circle]
                else:
                    while i_half_each < (left_num + number_chosen):
                        a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                        position_random = int(a_position)
                        determine_value = (position_random == sp.array(position_each_shift))
                        while determine_value.any() == True:
                            a_position = np.random.uniform(begin_point_shift, 
                                    end_point_shift)
                            position_random = int(a_position)
                            determine_value = (position_random == position_each_shift)
                        else:
                            position_each_shift.append(position_random)
                        i_half_each += 1
        position_half_shift[ind] = sp.empty(left_num_shift,int)
        #position_half_shift[ind] = sp.empty(len(position_each_shift),int)
        #print len(position_half_shift[ind]), len(position_each_shift)
        for i_position in sp.arange(len(position_half_shift[ind])):
            position_half_shift[ind][i_position] = position_each_shift[i_position]   

    x_position_random = [0] * len(onumber_fiber_blend)
    y_position_random = [0] * len(onumber_fiber_blend)
    x_position_random_shift = [0] * len(onumber_fiber_blend)
    y_position_random_shift = [0] * len(onumber_fiber_blend)
   #print 'the blend number', onumber_fiber_blend
    
    for ind in sp.arange(len(onumber_fiber_blend)):
        x_position_random[ind] = x_position[ind][position_half[ind]]
        y_position_random[ind] = y_position[ind][position_half[ind]]
        x_position_random_shift[ind] = x_position_shift[ind][position_half_shift[ind]]
        y_position_random_shift[ind] = y_position_shift[ind][position_half_shift[ind]]
    x_p = sp.zeros(len(x_position_random[0]))
    y_p = sp.zeros(len(x_position_random[0]))
    x_p[:] = x_position_random[0][:]
    y_p[:] = y_position_random[0][:]
    x_p_s = sp.zeros(len(x_position_random_shift[0]))
    y_p_s = sp.zeros(len(x_position_random_shift[0]))
    x_p_s[:] = x_position_random_shift[0][:]
    y_p_s[:] = y_position_random_shift[0][:]

    x_position = sp.empty(onumber_fiber, float)
    y_position = sp.empty(onumber_fiber, float)
    all_radius_fiber = sp.empty(onumber_fiber, float)
    fiber_kind = sp.empty(onumber_fiber, int)
    start = 0
    fiber_kind[:onumber_fiber_blend[0]] = 0
    for ind in sp.arange(len(onumber_fiber_blend)):
        for (t1array, t2array, t3array)  in \
                [(x_position, x_position_random, x_position_random_shift), 
                 (y_position, y_position_random, y_position_random_shift)]:
            t1array[start:start+len(x_position_random[ind])]\
                    = t2array[ind][:]
            t1array[start+len(x_position_random[ind]):start+onumber_fiber_blend[ind]]\
                    = t3array[ind][:]
        start = start+onumber_fiber_blend[ind]
    for j in range(len(onumber_fiber_blend)):
        all_radius_fiber[np.sum(onumber_fiber_blend[:j]):np.sum(onumber_fiber_blend[:j+1])]\
            = radius_fiber[j]
        fiber_kind[onumber_fiber_blend[j-1]:onumber_fiber_blend[j] + onumber_fiber_blend[j-1]] = j
    #plot starting overlap
    #plot_yarn(x_position, y_position, radius_fiber)
    filename = utils.OUTPUTDIR + os.sep + "combine.gz"
    # we now make sure the points are no longer overlapping, we detect overlap,
    # and use a dynamic scheme to correct the positions
    x_position_alpha = sp.empty(len(x_position), float)
    y_position_alpha = sp.empty(len(y_position), float)
    x_position_alpha[:] = x_position[:]
    y_position_alpha[:] = y_position[:]

    move_fibers_alpha(x_position_alpha, y_position_alpha, all_radius_fiber, 
                    oradius_yarn, omean_deviation)
    move_fibers_nonoverlap(x_position, y_position, all_radius_fiber, oradius_yarn, 
                        fiber_kind, omean_deviation)
   
##    x_polyester = []
##    y_polyester = []
##    x_cotton = []
##    y_cotton = []
##    radius_poly = []
##    radius_cotton = []
##    for i_fib in sp.arange(len(x_position)):
##        if fiber_kind[i_fib] == 0:
##            x_polyester.append(x_position[i_fib])
##            y_polyester.append(y_position[i_fib])
##            radius_poly.append(all_radius_fiber[i_fib])
##        else:
##            x_cotton.append(x_position[i_fib])
##            y_cotton.append(y_position[i_fib])
##            radius_cotton.append(all_radius_fiber[i_fib])
##    x_polyester = np.array(x_polyester)
##    y_polyester = np.array(y_polyester)
##    radius_poly = np.array(radius_poly)
##    x_cotton = np.array(x_cotton)
##    y_cotton = np.array(y_cotton)
##    radius_cotton = np.array(radius_cotton)
##    fig = pylab.figure()
##    ax = fig.add_subplot(111, xlim = (-1.1, 1.1), ylim = (-1.1, 1.1))
    
##    patches_1 = []
##    patches_2 = []
##    for x_center, y_center, radii in zip(x_polyester[:], y_polyester[:], radius_poly[:]):
##        circle = Circle((x_center, y_center), radii, facecolor = 'g', alpha = 0.4)
##        patches_1.append(circle)
##    angle = sp.zeros(len(x_cotton), float)
##    for i_cotton in sp.arange(len(x_cotton)):
##        angle[i_cotton] = np.random.uniform(0.0, 180.0)
##    for x_center, y_center, radii, angle in zip(x_cotton[:], y_cotton[:], radius_cotton[:], angle[:]):
##        ellipse = Ellipse((x_center, y_center), radii, radii*1.5, angle = angle, alpha = 0.4)
##        patches_2.append(ellipse)
##    p_1 = PatchCollection(patches_1, facecolor = 'red', cmap = matplotlib.cm.jet, alpha = 0.4)
##    p_2 = PatchCollection(patches_2, facecolor = 'black', cmap = matplotlib.cm.jet, alpha = 0.4)
##    ax.add_collection(p_1)
##    ax.add_collection(p_2)    

    filename_2 = utils.OUTPUTDIR + os.sep + "proportion_vl_overlap.gz"
    filename_3 = utils.OUTPUTDIR + os.sep + "proportion_vl_overlap_alpha.gz"

##    zone_position_ov, ratio_vl_ov = calculate_proportion(oradius_yarn, all_radius_fiber, 
##                                    x_position, y_position)
##    zone_position_ov_1, ratio_vl_ov_1 = calculate_proportion(oradius_yarn, all_radius_fiber,
##                                    x_position_alpha, y_position_alpha)

##    dump.write({'zone_position': zone_position_ov, 'ratio_value':ratio_vl_ov},
##                filename = filename_2, extension = '.gz')
##    dump.write({'zone_position_alpha': zone_position_ov_1, 'ratio_value_alpha':
##                ratio_vl_ov_1}, filename = filename_3, extension = '.gz')
    
#record the ratio value from the new layout with removing the overlap
    ratio_each_kind = [0] * len(onumber_fiber_blend)
    ratio_each_alpha= [0] * len(onumber_fiber_blend)
    zone_position_each = [0] 
    for i_kind in range(len(onumber_fiber_blend)):
        x_each_kind = []
        y_each_kind = []
        x_each_kind_alpha = []
        y_each_kind_alpha = []
        radius_each_kind = []
        radius_each_kind_alpha = []
        for i_fiber in range(len(fiber_kind)):
            if fiber_kind[i_fiber] == i_kind:
                x_each_kind.append(x_position[i_fiber])
                y_each_kind.append(y_position[i_fiber])
                radius_each_kind.append(all_radius_fiber[i_fiber])
                ##for the alpha value is equal to -1
                x_each_kind_alpha.append(x_position_alpha[i_fiber])
                y_each_kind_alpha.append(y_position_alpha[i_fiber])
                radius_each_kind_alpha.append(all_radius_fiber[i_fiber])
        x_each_kind = sp.array(x_each_kind)
        y_each_kind = sp.array(y_each_kind)
        x_each_kind_alpha = sp.array(x_each_kind_alpha)
        y_each_kind_alpha = sp.array(y_each_kind_alpha)
        radius_each_kind = np.array(radius_each_kind)
        sum_each_kind = sp.sum(sp.pi * sp.power(radius_each_kind, 2.))
##        print 'the radius value in the calculation', radius_each_kind[0]
##        print 'the sum of each kind fiber is', sum_each_kind
        #raw_input("check the value from the simulation, especially the first one")
        zone_position_ov, ratio_vl_ov = calculate_proportion(oradius_yarn, 
                                        radius_each_kind, x_each_kind, y_each_kind)
        zone_position_ov_alpha, ratio_ov_alpha = calculate_proportion(oradius_yarn, 
                                        radius_each_kind_alpha, x_each_kind_alpha,
                                        y_each_kind_alpha)
        ratio_each_kind[i_kind] = sp.zeros(len(ratio_vl_ov), float)
        ratio_each_alpha[i_kind] = sp.zeros(len(ratio_ov_alpha), float)
        ratio_each_kind[i_kind][:] = ratio_vl_ov[:]
        ratio_each_alpha[i_kind][:] = ratio_ov_alpha[:]
        zone_position_each = sp.zeros(len(zone_position_ov), float)
        zone_position_each[:] = zone_position_ov[:]
##        print 'the ratio value for each kind', ratio_each_kind
##        raw_input('stop or continue')
##        print 'length of the position', len(zone_position_ov)
##        zone_position_fit = sp.ones(len(zone_position_ov) + 1, float)
##        zone_position_fit[0] = 0.
##        zone_position_fit[1:] = zone_position_ov[:]
##        zone_position_fit[-1] = 1.0
##        
##        ratio_to_fit = sp.zeros(len(ratio_vl_ov) + 1, float)
##        ratio_to_fit[:-1] = ratio_vl_ov[:]
##        ratio_to_fit[-1] = 0. 
##        poly_fit_random = np.poly1d(np.polyfit(zone_position_fit, ratio_to_fit, 4))
##        print 'the polynomial fitting for the layout generated randomly', np.polyfit(zone_position_fit, ratio_to_fit, 4)
##        draw_fit = sp.linspace(0., 1.0, 50)
##        pylab.figure()
##        pylab.plot(draw_fit, poly_fit_random(draw_fit), '--')
##        pylab.plot(zone_position_fit, ratio_to_fit, 'o')
##        pylab.xlim(0., 1.01)
##        pylab.ylim(0., 0.8)
##        pylab.axis()
##        pylab.show()
##        raw_input( "check the fitting")

        dump.write({'zone_position':zone_position_ov, 'ratio_value': ratio_vl_ov},
                    filename = utils.OUTPUTDIR + os.sep + "each_kind_ratio_%g"%(i_kind), 
                    extension = '.gz')
        dump.write({'zone_position': zone_position_ov_alpha, 'ratio_value':ratio_ov_alpha},
                    filename = utils.OUTPUTDIR + os.sep + "each_kind_ratio_alpha_%g"%(i_kind),
                    extension = '.gz')
    ##raw_input("wait")
##    print 'before return the value', fiber_kind
##    raw_input("check the fiber kind value")
    return (x_position, y_position, all_radius_fiber, fiber_kind, type_fiber, ratio_each_kind, 
                zone_position_each, ratio_each_alpha)

def determine_overlap(xpos, ypos, radin, average_mean_deviation):
    """
    Determine if there is overlap between circles
    """
    #create array with the points to compare
    coord = sp.empty((len(xpos), 2))
    coord[:, 0] = xpos[:]
    coord[:, 1] = ypos[:]
    ind, res = fullcompare_array(coord, func=circledist, funcdata_a=radin*(1+ average_mean_deviation / 2.)) ##(NONTOUCH_FAC-1)/2))
    return ind, res[0], res[1]

def move_fibers_nonoverlap(xpos, ypos, radin, rad_yarn, fiber_kind, mean_deviation):
    ok = False
    nrmoves = 0
    nrmovesmax = 5000
    nrmovesmaxyarn = 10000
    average_mean_deviation = sp.average(mean_deviation)
    while not ok:
        ok = True
        nrmoves += 1
        (tooclose, dist, distreq) = determine_overlap(xpos, ypos, 
                                    radin, average_mean_deviation)
        #move the fibers
        vx = sp.zeros(len(tooclose), float)
        vy = sp.zeros(len(tooclose), float)
##        fig = pylab.figure()
##        ax = fig.add_subplot(111, xlim = (-1.1, 1.1), ylim = (-1.1, 1.1))
##        patches_1 = []
##        patches_2 = []
##        for x_center, y_center, radii, fib_type in zip(xpos, ypos, radin, fiber_kind):
##            if fib_type == 0:
##                circle = Circle((x_center, y_center), radii)
##                patches_1.append(circle)
##            else:
##                circle = Circle((x_center, y_center), radii)
##                patches_2.append(circle)
##        ax.xaxis.set_ticks(np.arange(-1., 1.1, 0.5))
##        ax.xaxis.set_ticklabels(["-$R$","-0.5$R$","0","0.5$R$", "R"])
##        ax.yaxis.set_ticks(np.arange(-1., 1.1, 0.5))
##        ax.yaxis.set_ticklabels(["-$R$","-0.5$R$","0","0.5$R$", "R"])
##        circle = Circle((0., 0.), 1.0)
##        patches_1.append(circle)
##        p1 = PatchCollection(patches_1, facecolor = 'red', cmap = matplotlib.cm.jet, alpha = 0.4)
##        p2 = PatchCollection(patches_2, facecolor = 'black', cmap = matplotlib.cm.jet, alpha = 0.4)
##        ax.add_collection(p1)
##        ax.add_collection(p2)
##        pylab.savefig('/home/lipei/PhDwork/porous_media/my_article/save_fig/geom%04d.png'%(nrmoves))
        for ind in sp.arange(len(tooclose)):
            for ov_ind, ov_distreal, ov_distreq in zip(tooclose[ind],
                                dist[ind], distreq[ind]):
                if ov_ind == ind:
                    continue
                #see if overlap is sufficiently small
                if (ov_distreq - ov_distreal)/ov_distreq > (NONTOUCH_FAC-1)/2:
                    #sufficient separation to fit in a fiber
                    ok = False
                if ov_distreal == 0.:
                    if ind < ov_ind:
                        delta_dir_x = xpos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                        delta_dir_y = ypos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                    else:
                        delta_dir_x = -xpos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                        delta_dir_y = -ypos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                    dirx = (ov_distreal - ov_distreq) * \
                            (radin[ov_ind] / (radin[ov_ind] + radin[ind]))* \
                            delta_dir_x
                    diry = (ov_distreal - ov_distreq) * \
                            (radin[ov_ind] / (radin[ov_ind] + radin[ind]))* \
                            delta_dir_y
                else:
                    dirx = (ov_distreal - ov_distreq) * \
                            (radin[ov_ind] / (radin[ov_ind] + radin[ind]))* \
                            (xpos[ov_ind] - xpos[ind])/ (ov_distreal)
                    diry = (ov_distreal - ov_distreq) * \
                            (radin[ov_ind] / (radin[ov_ind] + radin[ind]))* \
                            (ypos[ov_ind] - ypos[ind])/ (ov_distreal)
                
                vx[ind] += dirx
                vy[ind] += diry
            distance_central = sp.sqrt(xpos[ind]**2 + 
                                ypos[ind]**2)
            if distance_central > rad_yarn-radin[ind]:
                ok = False
                delta_central =  rad_yarn - radin[ind] - distance_central
                dirx = delta_central * xpos[ind] / sp.sqrt(xpos[ind]**2
                        + ypos[ind]**2.)
                diry = delta_central * ypos[ind] / sp.sqrt(xpos[ind]**2 + 
                        ypos[ind]**2)
                vx[ind] += dirx
                vy[ind] += diry
        
        if ok:
            print 'Good fiber layout found after', nrmoves, 'moves'
            break
        FRACVAL = 1.05 # a good value ??
        for ind in sp.arange(len(tooclose)):
            lambda_value = np.random.uniform(0,0.001)
            xpos[ind] = xpos[ind] + vx[ind] * (FRACVAL + lambda_value)
            ypos[ind] = ypos[ind] + vy[ind] * (FRACVAL + lambda_value)
        
        #break if too many iterations
        if nrmoves == nrmovesmax:
            print 'ERROR: no good solution found, breaking loop'
        

def move_fibers_alpha(xpos, ypos, radin, rad_yarn, mean_deviation):
    ##change the alpha value
    ok = False
    alpha = -0.
    alphafactor = 1/(1-alpha)
    nrmoves = 0
    nrmovesmax = 5000
    nrmovesmaxyarn = 10000
    average_mean_deviation = sp.average(mean_deviation)
    while not ok:
        ok = True
        nrmoves += 1
        (tooclose, dist, distreq) = determine_overlap(xpos, ypos, 
                                    radin, average_mean_deviation)
        vx = sp.zeros(len(tooclose), float)
        vy = sp.zeros(len(tooclose), float)
        for ind in sp.arange(len(tooclose)):
            for ov_ind, ov_distreal, ov_distreq in zip(tooclose[ind],
                                dist[ind], distreq[ind]):
                if ov_ind == ind:
                    continue
                #see if overlap is sufficiently small
                if (ov_distreq - ov_distreal)/ov_distreq > (NONTOUCH_FAC-1)/2:
                    #sufficient separation to fit in a fiber
                    ok = False
                if ov_distreal == 0.:
                    if ind < ov_ind:
                        delta_dir_x = xpos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                        delta_dir_y = ypos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                    else:
                        delta_dir_x = -xpos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                        delta_dir_y = -ypos[ov_ind] / sp.sqrt(sp.power(xpos[ov_ind], 2.) + 
                                                    sp.power(ypos[ov_ind], 2.))
                    dirx = (ov_distreal - ov_distreq) * alphafactor *\
                            delta_dir_x
                    diry = (ov_distreal - ov_distreq) * alphafactor *\
                            delta_dir_y
                else:
                    dirx = (ov_distreal - ov_distreq) * alphafactor *\
                            (xpos[ov_ind] - xpos[ind])/ (ov_distreal)
                    diry = (ov_distreal - ov_distreq) * alphafactor *\
                            (ypos[ov_ind] - ypos[ind])/ (ov_distreal)
                
                vx[ind] += dirx
                vy[ind] += diry
            distance_central = sp.sqrt(xpos[ind]**2 + 
                                ypos[ind]**2)
            if distance_central > rad_yarn-radin[ind]:
                ok = False
                delta_central =  rad_yarn - radin[ind] - distance_central
                dirx = delta_central * xpos[ind] / sp.sqrt(xpos[ind]**2
                        + ypos[ind]**2.)
                diry = delta_central * ypos[ind] / sp.sqrt(xpos[ind]**2 + 
                        ypos[ind]**2)
                vx[ind] += dirx
                vy[ind] += diry
        
        if ok:
            print 'After changing the Alpha value, Good fiber layout found after', nrmoves, 'moves'
            break
        FRACVAL = 1.05 # a good value ??
        for ind in sp.arange(len(tooclose)):
            lambda_value = np.random.uniform(0,0.001)
            xpos[ind] = xpos[ind] + vx[ind] * (FRACVAL + lambda_value)
            ypos[ind] = ypos[ind] + vy[ind] * (FRACVAL + lambda_value)
        
        #break if too many iterations
        if nrmoves == nrmovesmax:
            print 'ERROR: no good solution found, breaking loop'
            break
        
def prob_func_VL(r, radius_yarn, theta, beta):
    return (1 - 2 * theta) * sp.power((sp.exp(1.) - sp.exp(r / radius_yarn)) /
                            (sp.exp(1.) - 1.), beta) + theta

def plot_yarn(x_position, y_position, radius_fiber, fiber_kind=None, title=None):
    """
    Function to make a nice plot of the yarn with all the fibers
    """
    colors = ['blue', 'red', 'black','green']
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlim = (-1.1, 1.1), ylim = (-1.1, 1.1))
    max_kind = 0
    if fiber_kind is not None:
        max_kind = np.max(fiber_kind)
    patches = [None] * (max_kind+1)
    #each fiber is drawn
    if fiber_kind is not None:
        for x_center, y_center, radii, kind in zip(x_position, y_position, 
                                                   radius_fiber, fiber_kind):
            circle = Circle((x_center, y_center), radii)
            patch = patches[kind]
            if patch is None:
                patches[kind] = []
                patch = patches[kind]
            patch.append(circle)
    else: 
        patches[0] = []
        patch = patches[0]
        for x_center, y_center, radii in zip(x_position, y_position, radius_fiber):
            circle = Circle((x_center, y_center), radii)
            patch.append(circle)
    #add the yarn
    ax.xaxis.set_ticks(np.arange(-1., 1.1, 0.5))
    ax.xaxis.set_ticklabels(["-$R$","-0.5$R$","0","0.5$R$", "R"])
    ax.yaxis.set_ticks(np.arange(-1., 1.1, 0.5))
    ax.yaxis.set_ticklabels(["-$R$","-0.5$R$","0","0.5$R$", "R"])
##    fixup_subplot(ax.color)
##    fig.get_gca()
##    fig.clean()
##    fig.set_xlim(-1.1, 1.1)
##    fig.xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
##    fig.xticks_label('R')
##    fig.axes.set_xlabel('$X$')
    circle_1 = Circle((0., 0.), 1.0)
    patch.append(circle_1)
    for ind, patch in enumerate(patches):
        p = PatchCollection(patch, cmap = matplotlib.cm.jet, alpha = 0.4,
                            color=colors[ind])
        ax.add_collection(p)
    if title:
        pylab.title(title)
    pylab.ion()
    pylab.show()
