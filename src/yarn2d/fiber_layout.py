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
from matplotlib.patches import Circle, Wedge, Polygon
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
from yarn2d.config import FIBERLAYOUTS, FIBERSHAPE
from virtlocgeom import *
from arearatioprobability import *

#-------------------------------------------------------------------------
#
# three types of layout
#
#-------------------------------------------------------------------------

NONTOUCH_FAC = 1.01
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
    print 'the fiber kind value', fiber_kind

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
         list integers indicating kind of fiber)
    """
    #type_fiber = len(options.number_fiber_blend)
    #print 'kind of fiber', type_fiber
    x_central = options.get('x_central', 0.)
    y_central = options.get('y_central', 0.)
    onumber_fiber = options.get('number_fiber', 1)
    onumber_fiber_blend = options.get('number_fiber_blend', [1])
    oradius_fiber = options.get('radius_fiber', [1.])
    oradius_yarn = options.get('radius_yarn', 2.)
    ofirst_center = options.get('radius_first_center', 0.0)
    otheta_value = options.get('theta_value', 0.1)
    obeta_value = options.get('beta_value', 0.05)
    filename_1 = utils.OUTPUTDIR + os.sep + "proportion_vl_value.gz"
    if len(oradius_fiber) > 1 or len(onumber_fiber_blend) > 1:
        print 'ERROR: for virtual location layout the number of fibers must be 1.'
        print 'Actual number of type of fiber %d' %len(oradius)
        assert False
    if onumber_fiber != onumber_fiber_blend[0]:
        print 'ERROR: number fiber and blend do not correspond'
        assert False
    oradius_fiber = oradius_fiber[0] 
    
    #number_for_circles = [0] 
    #total_circles = [0]  
    #print 'radius_circle_central', len(radius_circle_central)
    #number_circle = [0] 
    #radius_fiber = [0] * type_fiber
    #total_number_vl = [0] * type_fiber
    #for ind in sp.arange(type_fiber):
    print 'ofirst_center', ofirst_center
    print 'oradius_fiber', oradius_fiber
    if ofirst_center != 0. :
        total_circles = int((oradius_yarn - (ofirst_center + 2. * 
                                    oradius_fiber*NONTOUCH_FAC))
                                    / (2 * oradius_fiber*NONTOUCH_FAC)) + 1
    else:
        total_circles = int((oradius_yarn - oradius_fiber * NONTOUCH_FAC)
                        / (2. * oradius_fiber * NONTOUCH_FAC)) + 1
    
    number_circle_central = sp.empty(total_circles, int)
    radius_circle_central = sp.empty(total_circles, float)
    
    #one fiber, fixed kind and radius
    radius_fiber = sp.ones(onumber_fiber_blend[0], float)*oradius_fiber
    fiber_kind = sp.zeros(onumber_fiber_blend[0], int)
    ind = sp.arange(onumber_fiber)
    if ofirst_center == 0:
        for i_circle in sp.arange(total_circles):
            radius_circle_central[i_circle] =  (i_circle * oradius_fiber * 2
                                                * NONTOUCH_FAC + ofirst_center)
            number_circle_central[i_circle]= max([int(sp.pi
                                            * radius_circle_central[i_circle]
                                            / (NONTOUCH_FAC * oradius_fiber)), 1])
            print 'the number of virtual locations in each zone', number_circle_central[i_circle]
    else:
        for i_circle in sp.arange(total_circles):
            radius_circle_central[i_circle] = ((2 * i_circle + 1) * oradius_fiber 
                                                * NONTOUCH_FAC + ofirst_center)
            number_circle_central[i_circle]= max([int(sp.pi
                                            * radius_circle_central[i_circle]
                                            / (NONTOUCH_FAC * oradius_fiber)), 1])
            print 'the number of virtual locations in each zone', number_circle_central[i_circle]
    total_number_vl = sum(number_circle_central[:])
    print 'the radius of each circle zone', radius_circle_central
    
    if onumber_fiber > total_number_vl:
        print 'ERROR: the number of fiber is more than the virtual locations'
        print 'the total number of virtual locations', total_number_vl
        sys.exit(0)
    #calculate the postion of each virtual location
    #x_position_vl = [0] * type_fiber
    #y_position_vl = [0] * type_fiber
    #probability_value = [0] * type_fiber
    #each_circle_zone_num = [0] * type_fiber
    #total_number_fiber = [0] * type_fiber

    probability_value = sp.ones(total_circles, float)
    each_circle_zone_num = sp.zeros(total_circles, int)
    x_position_vl = []
    y_position_vl = []
    for i_circle in sp.arange(total_circles):
        #if i_circle == 0:
        #    x_position_vl.append(0.) 
        #    y_position_vl.append(0.)
        #else:
        each_circle = number_circle_central[i_circle]
        if each_circle == 1:
            x_position_vl.append(0.)
            y_position_vl.append(0.)
        else:    
            for i_position in sp.arange(each_circle):
                x_position_t = radius_circle_central[i_circle] * sp.cos(2 * i_position * 
                            sp.arcsin(1.01 * oradius_fiber / radius_circle_central[i_circle]))
                y_position_t = radius_circle_central[i_circle] * sp.sin(2 * i_position * 
                            sp.arcsin(1.01 * oradius_fiber / radius_circle_central[i_circle]))
                x_position_vl.append(x_position_t)
                y_position_vl.append(y_position_t)
        #calculate the distribution value in each circle zone
        probability_value[i_circle] =  (1 - 2 * otheta_value) * sp.power(
                                    (sp.exp(1.) -sp.exp(radius_circle_central[i_circle] / 
                                    oradius_yarn))/(sp.exp(1) - 1), 
                                    obeta_value) + otheta_value        
        each_circle_zone_num[i_circle] = int(round(number_circle_central[i_circle] * 
                                                probability_value[i_circle]))
    print 'the theta_value', otheta_value
    print 'the beta value', obeta_value
    total_number_fiber = sum(each_circle_zone_num[:])
    #distribute the fiber in the yarn:
    if onumber_fiber > total_number_fiber:
        print 'total number fiber', total_number_fiber
        print 'ERROR: the probability function can not distribute all the fibers in the circle zones'
        sys.exit(0)
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
    
    while number_fiber_in_loop > 0:
        number_fiber_in_loop = number_fiber_in_loop - each_circle_zone_num[i_determine]
        i_determine += 1
        circle_loop += 1
    #make the fiber reach the radius of yarn
    i_type = len(onumber_fiber_blend)
    ##if i_type == 1:
    while circle_loop < total_circles:
        obeta_value += obeta_value 
        print 'inc beta', obeta_value
        circle_loop = 0
        i_determine = 0
        number_fiber_in_loop = onumber_fiber_blend[0]
        for i_circle in sp.arange(total_circles):
            probability_value[i_circle] =  (1 - 2 * otheta_value) * sp.power(
                            (sp.exp(1.) -sp.exp(radius_circle_central[i_circle] / 
                            oradius_yarn))/(sp.exp(1) - 1), 
                            obeta_value) + otheta_value
            each_circle_zone_num[i_circle] = int(round(number_circle_central[i_circle] * 
                                                    probability_value[i_circle]))
        total_sum_fiber = sum(each_circle_zone_num[:])
        while total_sum_fiber < onumber_fiber:
            obeta_value = obeta_value * 3. / 4. 
            print 'dec beta', obeta_value
            for i_circle in sp.arange(total_circles):
                probability_value[i_circle] =  (1 - 2 * otheta_value) * sp.power(
                                (sp.exp(1.) -sp.exp(radius_circle_central[i_circle] / 
                                oradius_yarn))/(sp.exp(1) - 1), 
                                obeta_value) + otheta_value
                each_circle_zone_num[i_circle] = int(round(number_circle_central[i_circle] * 
                                                        probability_value[i_circle]))
            total_sum_fiber = sum(each_circle_zone_num[:])
        else:
            while number_fiber_in_loop > 0:
                number_fiber_in_loop = number_fiber_in_loop - each_circle_zone_num[i_determine]
                i_determine += 1
                circle_loop += 1
    else:
        print 'the fiber distribution can reach the radius of yarn'
    #write the locations
    print 'the probability_value is', probability_value
    print 'each zone has the number of fibers', each_circle_zone_num
    print 'Finally the beta value', obeta_value 
    print 'the theta value', otheta_value
    filename = utils.OUTPUTDIR + os.sep + "virlayout.gz"    
    x_position = sp.empty(onumber_fiber_blend[0])
    y_position = sp.empty(onumber_fiber_blend[0])
    determine_generate = 0 #the number of generated fiber
    index_position = 0
    i_circle_number = 0
    number_fiber_in_loop = onumber_fiber_blend[0]
    while i_circle_number < circle_loop: 
        if i_circle_number < circle_loop - 1:
            location_number = sp.zeros(each_circle_zone_num[i_circle_number]) - 1
            for i_index in sp.arange(each_circle_zone_num[i_circle_number]):
                if i_index == 0:
                    a_position = np.random.uniform(index_position, 
                                index_position + number_circle_central[i_circle_number])
                    random_position = int(a_position)
                    x_position[determine_generate] = x_position_vl[random_position]
                    y_position[determine_generate] = y_position_vl[random_position]
                else:
                    a_position = np.random.uniform(index_position, 
                                index_position + number_circle_central[i_circle_number])
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
            
        elif i_circle_number == circle_loop - 1:
            location_number = sp.zeros(number_fiber_in_loop) - 1
            for i_index in sp.arange(number_fiber_in_loop):
                #print 'number_fiber_in_loop', number_fiber_in_loop[i_type]
                if i_index == 0:
                    a_position = np.random.uniform(index_position, index_position +
                                number_circle_central[i_circle_number])
                    random_position = int(a_position)
                    location_number[i_index] = random_position
                    x_position[determine_generate] = x_position_vl[random_position]
                    y_position[determine_generate] = y_position_vl[random_position]
                    #determine_generate += 1
                else:
                    a_position = np.random.uniform(index_position, index_position +
                                number_circle_central[i_circle_number])
                    random_position = int(a_position)
                    determine_value = (random_position == location_number)
                    while determine_value.any() == True:
                        a_position = np.random.uniform(index_position, 
                                    index_position + number_circle_central[i_circle_number])
                        random_position = int(a_position)
                        determine_value = (random_position == location_number)
                    else:
                        location_number[i_index] = random_position
                        x_position[determine_generate] = x_position_vl[random_position]
                        y_position[determine_generate] = y_position_vl[random_position]
                determine_generate += 1
            index_position += number_circle_central[i_circle_number]
            i_circle_number += 1
    oradius_fiber_array = sp.zeros(onumber_fiber_blend[0], float)
    oradius_fiber_array[:] = oradius_fiber
    dump.write({'x_position':x_position, 'y_position':y_position, 'radius_fiber':oradius_fiber,
                },
                filename = filename, extension = '.gz')
    fig = plot_yarn(x_position, y_position, oradius_fiber_array)#, fiber_kind)
    #print x_position
    return (x_position, y_position, radius_fiber, fiber_kind, x_position_vl, 
            y_position_vl)

def virtlocoverlaplayout(options):
    """ Generate the fiber layout in the yarn in virtual locations with random distribution
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
    x_central = options.get('x_central', 0.)
    y_central = options.get('y_central', 0.)
    onumber_fiber = options.get('number_fiber', 1)
    onumber_fiber_blend = options.get('number_fiber_blend', [1])
    oradius_fiber = options.get('radius_fiber', [1.])
    oradius_yarn = options.get('radius_yarn', 2.)
    otheta_value = options.get('theta_value', 0.1)
    obeta_value = options.get('beta_value', [0.05])
    onumber_fiber_each = sp.zeros(len(onumber_fiber_blend), float)
    onumber_fiber_each[:] = onumber_fiber_blend[:]
    type_fiber = len(onumber_fiber_blend)
    
    x_position = [0] * type_fiber
    y_position = [0] * type_fiber
    radius_fiber = [0] * type_fiber
    fiber_kind = [0] * type_fiber

    x_position_shift = [0] * type_fiber
    y_position_shift = [0] * type_fiber
    radius_fiber_shift = [0] * type_fiber
    fiber_kind_shift = [0] * type_fiber
    
    for i_type in sp.arange(type_fiber):
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
                'beta_value' : obeta_value[i_type],
                }
        ouroptions['radius_first_center'] = 0.             
        ax_position, ay_position, aradius_fiber, afiber_kind, x_vl, y_vl = virtloclayout(ouroptions)
        x_position[i_type][:] = ax_position[:]
        y_position[i_type][:] = ay_position[:]
        radius_fiber[i_type] = aradius_fiber[:]
        fiber_kind[i_type][:] = afiber_kind[:]
        dump.write({'x_position':ax_position, 'y_position':ay_position, 'radius_fiber':aradius_fiber,
                    },
                    filename = filename, extension = '.gz')
        zone_position, ratio_fiber = calculate_proportion(oradius_yarn, aradius_fiber,
                                ax_position, ay_position)
        dump.write({'zone_position': zone_position, 'ratio_value': ratio_fiber},
                    filename = filename_1, extension = '.gz')
        
    #for the shift part
    for i_type in sp.arange(type_fiber):
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
                'beta_value' : obeta_value[i_type],
        }
        ouroptions['radius_first_center'] = 0.5*oradius_fiber[i_type]
        ax_position_shift, ay_position_shift, aradius_fiber_shift, afiber_kind_shift, x_vl_shift,y_vl_shift = virtloclayout(ouroptions)
        x_position_shift[i_type][:] = ax_position_shift[:]
        y_position_shift[i_type][:] = ay_position_shift[:]
        radius_fiber_shift[i_type][:] = aradius_fiber_shift[:]
        fiber_kind_shift[i_type][:] = afiber_kind_shift[:]
        dump.write({'x_position':ax_position_shift, 'y_position':ay_position_shift, 
                    'radius_fiber':aradius_fiber_shift},
                    filename = filename, extension = '.gz')
        zone_position_shift, ratio_fiber_shift = calculate_proportion(oradius_yarn, aradius_fiber_shift,
                                ax_position_shift, ay_position_shift)
        dump.write({'zone_position': zone_position_shift, 'ratio_value': ratio_fiber_shift},
                    filename = filename_1, extension = '.gz')

    #take 50% points from each other
    number_vl_overlap = sp.empty(len(onumber_fiber_blend))
    for ind in sp.arange(len(onumber_fiber_blend)):
        number_vl_overlap[ind] = int(onumber_fiber_blend[ind] / 2.)
    position_half = [0] * type_fiber
    position_half_shift = [0] * type_fiber
    for ind in range(type_fiber):
        position_half[ind] = sp.empty(number_vl_overlap[ind], int)
        position_half_shift[ind] = sp.empty(onumber_fiber_blend[ind]-number_vl_overlap[ind], int)
        i_half = 0
        while i_half < number_vl_overlap[ind]:
            a_position = np.random.uniform(0., onumber_fiber_blend[ind])
            position_random = int(a_position)
            if i_half == 0:
                position_half[ind][i_half] = position_random
                #print 'self.position_random', position_random
            else:
                determine_value = (position_random == position_half[ind])
                while determine_value.any() == True:
                    a_position = np.random.uniform(0., onumber_fiber_blend[ind])
                    position_random = int(a_position)
                    determine_value = (position_random == position_half[ind])
                else:
                    position_half[ind][i_half] = position_random
                    #print 'self.position_random', position_random 
            i_half += 1
        
        i_half_1 = 0
        while i_half_1 < onumber_fiber_blend[ind]-number_vl_overlap[ind]:
            b_position = np.random.uniform(0., onumber_fiber_blend[ind])
            position_random = int(b_position)
            if i_half_1 == 0:
                position_half_shift[ind][i_half_1] = position_random
            else:
                determine_value = (position_random == position_half_shift[ind])
                while determine_value.any() == True:
                    b_position = np.random.uniform(0., onumber_fiber_blend[ind])
                    position_random = int(b_position)
                    determine_value = (position_random == position_half_shift[ind])
                else:
                    position_half_shift[ind][i_half_1] = position_random
            i_half_1 += 1
    x_position_random = [0] * len(onumber_fiber_blend)
    y_position_random = [0] * len(onumber_fiber_blend)
    x_position_random_shift = [0] * len(onumber_fiber_blend)
    y_position_random_shift = [0] * len(onumber_fiber_blend)
    for ind in sp.arange(len(onumber_fiber_blend)):
        x_position_random[ind] = x_position[ind][position_half[ind]]
        y_position_random[ind] = y_position[ind][position_half[ind]]
        x_position_random_shift[ind] = x_position_shift[ind][position_half_shift[ind]]
        y_position_random_shift[ind] = y_position_shift[ind][position_half_shift[ind]]
    
    plt.figure()
    plt.plot(x_position_random[0], y_position_random[0], 'o',
        x_position_random_shift[0], y_position_random_shift[0], 's')
    plt.axis([-1,1, -1, 1])
    plt.draw()
    
    # merge back into 1 array
    x_position = sp.empty(onumber_fiber, float)
    y_position = sp.empty(onumber_fiber, float)
    radius_fiber = sp.empty(onumber_fiber, float)
    fiber_kind = sp.empty(onumber_fiber, int)
    start = 0
    fiber_kind[:onumber_fiber_blend[0]] = 0
    for ind in sp.arange(len(onumber_fiber_blend)):
        for (t1array, t2array, t3array)  in \
                [(x_position, x_position_random, x_position_random_shift), 
                 (y_position, y_position_random, y_position_random_shift)]:
            t1array[start:start+number_vl_overlap[ind]]\
                    = t2array[ind][:]
            t1array[start+number_vl_overlap[ind]:start+onumber_fiber_blend[ind]]\
                    = t3array[ind][:]
        start = start+onumber_fiber_blend[ind]
    for j in range(len(onumber_fiber_blend)):
        radius_fiber[sum(onumber_fiber_blend[:j]):sum(onumber_fiber_blend[:j+1])]\
            = oradius_fiber[j]
        fiber_kind[onumber_fiber_blend[j-1]:onumber_fiber_blend[j] + onumber_fiber_blend[j-1]] = j
    #plot starting overlap
    plot_yarn(x_position, y_position, radius_fiber)
    filename = utils.OUTPUTDIR + os.sep + "combine.gz"
    # we now make sure the points are no longer overlapping, we detect overlap,
    # and use a dynamic scheme to correct the positions
    x_position_alpha = sp.empty(len(x_position), float)
    y_position_alpha = sp.empty(len(y_position), float)
    x_position_alpha[:] = x_position[:]
    y_position_alpha[:] = y_position[:]
    move_fibers_alpha(x_position_alpha, y_position_alpha, radius_fiber, oradius_yarn)
    move_fibers_nonoverlap(x_position, y_position, radius_fiber, oradius_yarn)
    
    dump.write({'x_position':x_position, 'y_position':y_position, 'radius_fiber':radius_fiber},
                filename = filename, extension = '.gz')
##    plot_yarn(x_position, y_position, radius_fiber)
##    pylab.show()
##    plot_yarn(x_position_alpha, y_position_alpha, radius_fiber)
##    pylab.show()
    filename_2 = utils.OUTPUTDIR + os.sep + "proportion_vl_overlap.gz"
    filename_3 = utils.OUTPUTDIR + os.sep + "proportion_vl_overlap_alpha.gz"
    zone_position_ov, ratio_vl_ov = calculate_proportion(oradius_yarn, radius_fiber, 
                                    x_position, y_position)
    zone_position_ov_1, ratio_vl_ov_1 = calculate_proportion(oradius_yarn, radius_fiber,
                                    x_position_alpha, y_position_alpha)
    dump.write({'zone_position': zone_position_ov, 'ratio_value':ratio_vl_ov},
                filename = filename_2, extension = '.gz')
    dump.write({'zone_position_alpha': zone_position_ov_1, 'ratio_value_alpha':
                ratio_vl_ov_1}, filename = filename_3, extension = '.gz')
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
                radius_each_kind.append(radius_fiber[i_fiber])
                ##for the alpha value is equal to -1
                x_each_kind_alpha.append(x_position_alpha[i_fiber])
                y_each_kind_alpha.append(y_position_alpha[i_fiber])
                radius_each_kind_alpha.append(radius_fiber[i_fiber])
        x_each_kind = sp.array(x_each_kind)
        y_each_kind = sp.array(y_each_kind)
        x_each_kind_alpha = sp.array(x_each_kind_alpha)
        y_each_kind_alpha = sp.array(y_each_kind_alpha)
        zone_position_ov, ratio_vl_ov = calculate_proportion(oradius_yarn, 
                                        radius_each_kind, x_each_kind, y_each_kind)
        zone_position_ov_alpha, ratio_ov_alpha = calculate_proportion(oradius_yarn, 
                                        radius_each_kind_alpha, x_each_kind_alpha,
                                        y_each_kind_alpha)
        dump.write({'zone_position':zone_position_ov, 'ratio_value': ratio_vl_ov},
                    filename = utils.OUTPUTDIR + os.sep + "each_kind_ratio_%g"%(i_kind), 
                    extension = '.gz')
        dump.write({'zone_position': zone_position_ov_alpha, 'ratio_value':ratio_ov_alpha},
                    filename = utils.OUTPUTDIR + os.sep + "each_kind_ratio_alpha_%g"%(i_kind),
                    extension = '.gz')
    raw_input("wait")
    return (x_position, y_position, radius_fiber, fiber_kind)

def determine_overlap(xpos, ypos, radin):
    """
    Determine if there is overlap between circles
    """
    #create array with the points to compare
    coord = sp.empty((len(xpos), 2))
    coord[:, 0] = xpos[:]
    coord[:, 1] = ypos[:]
    ind, res = fullcompare_array(coord, func=circledist, funcdata_a=radin*(1+(NONTOUCH_FAC-1)/2))
    return ind, res[0], res[1]

def move_fibers_nonoverlap(xpos, ypos, radin, rad_yarn):
    ok = False
    nrmoves = 0
    nrmovesmax = 500
    nrmovesmaxyarn = 1000
    while not ok:
        ok = True
        nrmoves += 1
        (tooclose, dist, distreq) = determine_overlap(xpos, ypos, 
                                    radin)
        #move the fibers
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
            break

def move_fibers_alpha(xpos, ypos, radin, rad_yarn):
    ##change the alpha value
    ok = False
    nrmoves = 0
    nrmovesmax = 500
    nrmovesmaxyarn = 1000
    while not ok:
        ok = True
        nrmoves += 1
        (tooclose, dist, distreq) = determine_overlap(xpos, ypos, 
                                    radin)
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
                    dirx = (ov_distreal - ov_distreq) * \
                            (1. / 2.)* \
                            delta_dir_x
                    diry = (ov_distreal - ov_distreq) * \
                            (1. / 2.)* \
                            delta_dir_y
                else:
                    dirx = (ov_distreal - ov_distreq) * \
                            (1. / 2.)* \
                            (xpos[ov_ind] - xpos[ind])/ (ov_distreal)
                    diry = (ov_distreal - ov_distreq) * \
                            (1. / 2.)* \
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

def plot_yarn(x_position, y_position, radius_fiber):#, fiber_kind):
    """
    Function to make a nice plot of the yarn with all the fibers
    """
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlim = (-1.1, 1.1), ylim = (-1.1, 1.1))
    patches = []
    #each fiber is drawn
    for x_center, y_center, radii in zip(x_position, y_position, radius_fiber):
        circle = Circle((x_center, y_center), radii)
        patches.append(circle)
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
    circle = Circle((0., 0.), 1.0)
    patches.append(circle)
    p = PatchCollection(patches, cmap = matplotlib.cm.jet, alpha = 0.4)
    ax.add_collection(p)
    pylab.draw()

def test():
    pass

if __name__ == '__main__': 
    test()
