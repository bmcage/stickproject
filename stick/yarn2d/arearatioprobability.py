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
import const
import numpy as np
import scipy as sp
import time


#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
from fipy import Gmsh2D
from fipy import *
from virtlocgeom import *

NONTOUCH_FAC = 1.01

def calculate_proportion(rad_yarn, rad_fib, x_fib, y_fib, nrzones=5):
    #divide the yarn zone to five concentric zones
    zone_radius = sp.zeros(nrzones, float)
    zone_width = sp.zeros(nrzones, float)
    
    width_zone = 2. * rad_yarn / ((nrzones-1)*2 + 1)
    print rad_yarn
    for i_circle in sp.arange(len(zone_radius)):
        zone_radius[i_circle] += width_zone * i_circle + width_zone / 2.
        zone_width[i_circle] = width_zone
    total_zone = sum(zone_width) - zone_radius[0]
    print 'the radius of each zone is (zone_radius)', zone_radius
    #zone_radius[-1] = rad_yarn * NONTOUCH_FAC
    #zone_width[-1] = zone_width[-1] + rad_yarn * (NONTOUCH_FAC -1.)
    print 'the sum of five zones', total_zone
    
    if total_zone > rad_yarn * NONTOUCH_FAC:
        assert False
    while abs(rad_yarn - total_zone) > 1.0e-4:
        diff = abs(rad_yarn - total_zone)
        if rad_yarn > total_zone:
            zone_radius[:] = zone_radius[:] + diff / nrzones
        else:
            zone_radius[:] = zone_radius[:] - diff / nrzones
    #count how many virtual locations in each zone
    #calculate the area of fiber cross section in each domain
    distan_fib_central = sp.sqrt(x_fib**2. + y_fib**2.)
    #print 'the distance from fibre to central', distan_fib_central
    area_fib_zone = sp.zeros(len(zone_radius))
    count_number = 0
    filename = utils.OUTPUTDIR + os.sep + "proportion_value.gz"
    i_fiber_calculation = 0
    for i_circle in sp.arange(len(zone_radius)):
        if i_circle == 0:
            for i_fib in sp.arange(len(distan_fib_central)):
                if distan_fib_central[i_fib] + rad_fib[i_fib] <= zone_radius[i_circle]:
                    area_fib_zone[i_circle] += sp.pi * (rad_fib[i_fib] ** 2.0)
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
                elif distan_fib_central[i_fib] < zone_radius[i_circle] and\
                 (distan_fib_central[i_fib] + rad_fib[i_fib]) > zone_radius[i_circle]:
                    #print 'most in the ring zone'
                    solution_points = intersect_circles(x_fib[i_fib],y_fib[i_fib],
                                zone_radius[i_circle],rad_fib[i_fib])
                    x_in_one = float(solution_points[0][0])
                    x_in_secon = float(solution_points[1][0])
                    y_in_one = float(solution_points[0][1])
                    y_in_secon = float(solution_points[1][1])
                    #print 'the solution for the formula', solution_points
                    square_distan = (x_in_one - x_in_secon)**2. + \
                                    (y_in_one - y_in_secon)**2.
                    distan_two_points = sp.sqrt(square_distan)
                    alpha_r = sp.arccos((2. * rad_fib[i_fib] ** 2. - 
                            distan_two_points ** 2.)
                     / (2.* rad_fib[i_fib] ** 2.))
                    beta_r_zone = sp.arccos((2. * zone_radius[i_circle] ** 2. - square_distan) 
                                / (2.* zone_radius[i_circle] ** 2.))
                    #print x_fib[i_fib], y_fib[i_fib]
                    #print distan_two_points, square_distan, rad_fib[i_fib], zone_radius[i_circle]
                    #print alpha_r, beta_r_zone
##                    if alpha_r <= sp.pi:
##                        area_circ = alpha_r / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
##                        area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
##                        area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.) - \
##                                    1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
##                        area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
##                                                area_curve
##                        part_fib = area_circ - area_triangle + area_curve
##                    else:
                    area_circ = (2. * sp.pi - alpha_r) / (2. * sp.pi) * sp.pi * (rad_fib[i_fib] ** 2.)
                    area_triangle = 1./2. * sp.sin(alpha_r) * (rad_fib[i_fib] ** 2.)
                    area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * (zone_radius[i_circle] ** 2.) - \
                        1. / 2. * (zone_radius[i_circle] ** 2.) * sp.sin(beta_r_zone)
                    area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ + area_triangle + \
                        area_curve
                    part_fib = area_circ + area_triangle + area_curve
                    #print area_circ, area_triangle, area_curve
                    area_fib_zone[i_circle +1] += sp.pi * (rad_fib[i_fib] ** 2.) - part_fib
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
                elif distan_fib_central[i_fib] > zone_radius[i_circle] and distan_fib_central[i_fib] \
                        - rad_fib[i_fib] < zone_radius[i_circle]:
                    #print 'less part in the ring zone'
                    solution_points = intersect_circles(x_fib[i_fib],y_fib[i_fib],
                                zone_radius[i_circle],rad_fib[i_fib])
                    x_in_one = solution_points[0][0]
                    x_in_secon = solution_points[1][0]
                    y_in_one = solution_points[0][1]
                    y_in_secon = solution_points[1][1]
                    square_distan = (x_in_one - x_in_secon)**2+ \
                                    (y_in_one - y_in_secon)**2
                    square_distan = float(square_distan)
                    distan_two_points = sp.sqrt(square_distan)
                    alpha_r = sp.arccos((2. * sp.power(rad_fib[i_fib], 2.) - sp.power(distan_two_points, 2.))
                     / (2.* sp.power(rad_fib[i_fib], 2.)))
                    beta_r_zone = sp.arccos((2. * (zone_radius[i_circle] ** 2.) - (distan_two_points ** 2.))/
                        (2. * (zone_radius[i_circle] ** 2.)))
                    #print alpha_r, beta_r_zone
                    #sp.power(rad_fib[i_fib], 2.)
                    area_circ = alpha_r / (2. * sp.pi) * sp.pi * (rad_fib[i_fib] ** 2.)
                    area_triangle = 1./2. * sp.sin(alpha_r) * (rad_fib[i_fib] ** 2.)
                    area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * (zone_radius[i_circle] ** 2.0) - \
                        1. / 2. * (zone_radius[i_circle] ** 2.) * sp.sin(beta_r_zone)
                    area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
                        area_curve
                    #print area_circ, area_triangle, area_curve
                    part_fib = area_circ - area_triangle + area_curve
                    area_fib_zone[i_circle + 1] += sp.pi * (rad_fib[i_fib] ** 2.) - part_fib
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
                    #print 'less than half of fiber in the central zone', area_fib_zone[i_circle]
                elif distan_fib_central[i_fib] == zone_radius[i_circle]:
                    #print 'center on the ring'
                    solution_points = intersect_circles(x_fib[i_fib],y_fib[i_fib],
                                zone_radius[i_circle],rad_fib[i_fib])
                    x_in_one = solution_points[0][0]
                    x_in_secon = solution_points[1][0]
                    y_in_one = solution_points[0][1]
                    y_in_secon = solution_points[1][1]
                    square_distan = (x_in_one - x_in_secon)**2. + \
                                    (y_in_one - y_in_secon)**2.
                    #square_distan = float(square_distan)
                    distan_two_points = sp.sqrt(square_distan)
                    alpha_r = sp.arccos((2. * (rad_fib[i_fib] ** 2.) - (distan_two_points ** 2.))
                     / (2.* (rad_fib[i_fib] ** 2.)))
                    beta_r_zone = sp.arccos((2. * (zone_radius[i_circle] ** 2.) - (distan_two_points ** 2.))
                    / (2. * (zone_radius[i_circle] ** 2)))
                    #sp.power(rad_fib[i_fib], 2.)
                    #print alpha_r, beta_r_zone
                    area_circ = (2. * sp.pi - alpha_r) / (2. * sp.pi) * sp.pi * (rad_fib[i_fib] ** 2.)
                    area_triangle = 1./2. * sp.sin(alpha_r) * (rad_fib[i_fib] ** 2.)
                    area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * (zone_radius[i_circle] ** 2.0) - \
                        1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
                    area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
                        area_curve
                    #print area_circ, area_triangle, area_curve
                    part_fib = area_circ - area_triangle + area_curve 
                    area_fib_zone[i_circle +1] += sp.pi * sp.power(rad_fib[i_fib], 2.) - part_fib
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
        else:
            for i_fib in sp.arange(len(distan_fib_central)):
                if distan_fib_central[i_fib] - rad_fib[i_fib]>= zone_radius[i_circle -1] and \
                    distan_fib_central[i_fib] + rad_fib[i_fib] <= zone_radius[i_circle]:
                    area_fib_zone[i_circle] += sp.pi * sp.power(rad_fib[i_fib], 2.)
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
                    #print 'the value in the zone area', area_fib_zone[i_circle]
                elif distan_fib_central[i_fib] < zone_radius[i_circle] and \
                    (distan_fib_central[i_fib] + rad_fib[i_fib]) > zone_radius[i_circle]:
                    solution_points = intersect_circles(x_fib[i_fib],y_fib[i_fib],
                                zone_radius[i_circle],rad_fib[i_fib])
                    x_in_one = float(solution_points[0][0])
                    x_in_secon = float(solution_points[1][0])
                    y_in_one = float(solution_points[0][1])
                    y_in_secon = float(solution_points[1][1])
                    square_distan = (x_in_one - x_in_secon)**2. + \
                                    (y_in_one - y_in_secon)**2.
                    distan_two_points = sp.sqrt(square_distan)
                    alpha_r = sp.arccos((2. * sp.power(rad_fib[i_fib], 2.) - sp.power(distan_two_points, 2.))
                     / (2.* sp.power(rad_fib[i_fib], 2.)))
                    beta_r_zone = sp.arccos((2. * sp.power(zone_radius[i_circle], 2.) 
                                - sp.power(distan_two_points, 2.)) / (2.* 
                                sp.power(zone_radius[i_circle], 2.)))
##                    if alpha_r <= sp.pi:
##                        area_circ = alpha_r / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
##                        area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
##                        area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.) - \
##                                    1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
##                        area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
##                                                area_curve
##                        part_fib = area_circ - area_triangle + area_curve
##                    else:
                    area_circ = (2. * sp.pi - alpha_r) / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                    area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                    area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.) - \
                        1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
                    area_fib_zone[i_circle] += area_circ + area_triangle + \
                        area_curve
                    part_fib = area_circ + area_triangle + area_curve
                    #print i_circle, len(area_fib_zone)
                    if i_circle + 1 >= len(zone_radius):
                        print 'the position of fiber', x_fib[i_fib], y_fib[i_fib]
                        print 'the index of the fiber', i_fib
                        print 'the distance to then center', distan_fib_central[i_fib] 
                        print 'plus the radius of the fiber', distan_fib_central[i_fib] + rad_fib[0]
                    if (distan_fib_central[i_fib] + rad_fib[i_fib]) < rad_yarn:
                        area_fib_zone[i_circle +1] += sp.pi * sp.power(rad_fib[i_fib], 2.) - part_fib
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
                elif distan_fib_central[i_fib] > zone_radius[i_circle] and distan_fib_central[i_fib] \
                        - rad_fib[i_fib] < zone_radius[i_circle]:
                    solution_points = intersect_circles(x_fib[i_fib],y_fib[i_fib],
                                zone_radius[i_circle],rad_fib[i_fib])
                    x_in_one = solution_points[0][0]
                    x_in_secon = solution_points[1][0]
                    x_in_secon = solution_points[1][0]
                    y_in_one = solution_points[0][1]
                    y_in_secon = solution_points[1][1]
                    square_distan = (x_in_one - x_in_secon)**2.+ \
                                    (y_in_one - y_in_secon)**2.
                    square_distan = float(square_distan)
                    distan_two_points = sp.sqrt(square_distan)
                    alpha_r = sp.arccos((2. * sp.power(rad_fib[i_fib], 2.) - sp.power(distan_two_points, 2.))
                     / (2.* sp.power(rad_fib[i_fib], 2.)))
                    beta_r_zone = sp.arccos((2. * sp.power(zone_radius[i_circle], 2.) - sp.power(distan_two_points, 2.))/
                        (2. * sp.power(zone_radius[i_circle], 2.)))
                    sp.power(rad_fib[i_fib], 2.)
                    area_circ = alpha_r / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                    area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                    area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.0) - \
                        1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)

                    area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
                        area_curve
                    part_fib = area_circ - area_triangle + area_curve
                    area_fib_zone[i_circle + 1] += sp.pi * sp.power(rad_fib[i_fib], 2.) - part_fib
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
                elif distan_fib_central[i_fib] == zone_radius[i_circle]:
                    solution_points = intersect_circles(x_fib[i_fib],y_fib[i_fib],
                                zone_radius[i_circle],rad_fib[i_fib])
                    x_in_one = solution_points[0][0]
                    x_in_secon = solution_points[1][0]
                    y_in_one = solution_points[0][1]
                    y_in_secon = solution_points[1][1]
                    square_distan = (x_in_one - x_in_secon)**2. + \
                                    (y_in_one - y_in_secon)**2.
                    distan_two_points = sp.sqrt(square_distan)
                    alpha_r = sp.arccos((2. * sp.power(rad_fib[i_fib], 2.) - sp.power(distan_two_points, 2.))
                     / (2.* sp.power(rad_fib[i_fib], 2.)))
                    beta_r_zone = sp.arccos((2. * sp.power(zone_radius[i_circle], 2.) - sp.power(distan_two_points, 2.))
                    / (2.*sp.power(zone_radius[i_circle], 2.)))
                    sp.power(rad_fib[i_fib], 2.)
                    area_circ = (2. * sp.pi - alpha_r) / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                    area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                    area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.0) - \
                        1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
                    area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
                        area_curve
                    part_fib = area_circ - area_triangle + area_curve
                    area_fib_zone[i_circle +1] += sp.pi * sp.power(rad_fib[i_fib], 2.) - part_fib
                    i_fiber_calculation += 1
                    #print 'area_fib_zone', area_fib_zone
                    #raw_input("confirm their is a value")
    total_fiber_area = sp.sum(area_fib_zone)
    print 'in each zone the area of fibers cross section', area_fib_zone
    print 'the total area of fiber in the yarn domain and calculation', total_fiber_area
    print 'the area sum from the fiber input', sp.sum(sp.pi * sp.power(rad_fib, 2.))
    #raw_input ("next time step <return>....")
    print 'the value of area of fiber in each zone', area_fib_zone
    #calculate the area value of each concentric zone
    print 'the number of fiber in the area calculation', i_fiber_calculation
    zone_area = sp.zeros(len(zone_radius), float)
    for i_zone in sp.arange(len(zone_radius)):
        if i_zone == 0:
            zone_area[i_zone] = sp.pi * sp.power(zone_radius[i_zone], 2.)
        else:
            zone_area[i_zone] = sp.pi * sp.power(zone_radius[i_zone], 2.) - \
                                sp.sum(zone_area[:i_zone])#sp.pi * sp.power(zone_radius[i_zone - 1], 2.)
    ratio_each = area_fib_zone / zone_area
    total_zone_area = sp.sum(zone_area)
    index_zone = sp.arange(1, len(zone_radius)+1, 1)
    print 'the area of the yarn cross-section', sp.pi * sp.power(rad_yarn, 2.)
    print 'the area of the sum of ring zones', total_zone_area
    print 'the value of zone area is', zone_area
    print 'the ratio value in each zone', ratio_each
    zone_point = sp.zeros(len(zone_radius))
    for i_circle in sp.arange(len(zone_radius)):
        zone_point[i_circle] = width_zone / 2. + i_circle * width_zone
    #print 'the zone central point value', zone_point
    #i_zone[:] = i_zone[:] + 1
##    dump.write({'zone_number': zone_point, 'ratio_value': ratio_each}, filename = 
##                filename, extension = '.gz')
    return (zone_point, ratio_each)

def plot_ratio_function(zone_position, each_time_ratio, type_fiber, 
                        probability_function):
    """ A function to nicely plot area proportions found after the iteration """
    iteration_times = int(len(each_time_ratio) / 2)
    each_kind_ratio = [0] * type_fiber
    each_coefficients = [0.0233, 0.0273]
    mean_ratio_subplot = [0] * type_fiber
    
    for i_type in sp.arange(type_fiber):
        each_kind_ratio[i_type] = []
        
        for i_iteration in sp.arange(iteration_times):
            index_same_kind = type_fiber * i_iteration + i_type
            each_kind_ratio[i_type].append(each_time_ratio[index_same_kind])
        each_kind_ratio[i_type] = sp.array(each_kind_ratio[i_type])
        mean_ratio = sp.zeros(len(zone_position), float)
        mean_ratio_subplot[i_type] = sp.zeros(len(zone_position), float)
        value_in_iteration = sp.zeros(iteration_times, float)
        for i_position in sp.arange(len(zone_position) - 1):
            for i_iteration in sp.arange(iteration_times):
                value_in_iteration[i_iteration] = each_kind_ratio[i_type][i_iteration][i_position]
            mean_ratio[i_position] = np.mean(value_in_iteration)
            mean_ratio_subplot[i_type][i_position] = mean_ratio[i_position]
        position_for_function = sp.linspace(0, zone_position[-1], 50)
        each_probability_function = probability_function[i_type]
        coefficient = each_coefficients[i_type]
        import pylab
        pylab.figure()
        pylab.plot(position_for_function, each_probability_function(position_for_function) / coefficient, '--')
        pylab.plot(zone_position, mean_ratio, 'o')
        pylab.plot()
        pylab.xlabel('Relative position in the yarn domain')
        pylab.ylabel('Probability value')
        pylab.xlim(0., 1.05)
        pylab.ylim(0., 0.90)
        pylab.title('Area probability for %d th kind of fiber'% (i_type+1))
    #part for drawing in subplot way
    pylab.ion()
    pylab.figure()
    for i_type in sp.arange(type_fiber):
        mean_for_plot = sp.zeros(len(mean_ratio_subplot[i_type]))
        mean_for_plot[:] = mean_ratio_subplot[i_type][:]
        each_probability_function = probability_function[i_type]
        pylab.subplot(type_fiber, 1, i_type + 1)
        pylab.title('Area probability for %d th kind of fiber'% (i_type+1))
        pylab.plot(zone_position, mean_for_plot, 'o')
        pylab.xlim(0., 1.05)
        pylab.ylim(0., 0.8)
        pylab.xlabel('Relative position in the yarn domain')
        pylab.ylabel('Probalility value')
        pylab.plot(position_for_function, each_probability_function(position_for_function) / coefficient, '--')
        pylab.axis()
    pylab.show()
        
        
        
    
        