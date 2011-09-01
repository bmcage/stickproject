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
import matplotlib.pyplot as plt
import time
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import pylab
import matplotlib
##import sympy
##from sympy.abc import x,y

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
from fipy import Gmsh2D
from fipy import *
from yarn2d.config import FIBERLAYOUTS
from yarn2d.config import FIBERSHAPE
from virtlocgeom import *

NONTOUCH_FAC = 1.01
def calculate_proportion(rad_yarn, rad_fib, x_fib, y_fib,):
    #divide the yarn zone to five concentric zones
    zone_radius = sp.zeros(8, float)
    zone_width = sp.zeros(8, float)
    width_zone = rad_yarn / 7.5
    print rad_yarn
    for i_circle in sp.arange(len(zone_radius)):
        zone_radius[i_circle] += width_zone * i_circle + width_zone / 2.
        zone_width[i_circle] = width_zone
    total_zone = sum(zone_width) - zone_radius[0]
    #zone_radius[-1] = rad_yarn * NONTOUCH_FAC
    #zone_width[-1] = zone_width[-1] + rad_yarn * (NONTOUCH_FAC -1.)
    print 'the sum of five zones', total_zone
    if total_zone > rad_yarn * NONTOUCH_FAC:
        assert False
    while abs(rad_yarn - total_zone) > 1.0e-4:
        diff = abs(rad_yarn - total_zone)
        if rad_yarn > total_zone:
            zone_radius[:] = zone_radius[:] + diff / 5.0
        else:
            zone_radius[:] = zone_radius[:] - diff / 5.0
    #count how many virtual locations in each zone
    #calculate the area of fiber cross section in each domain
    print 'zone_radius value', zone_radius
    distan_fib_central = sp.sqrt(x_fib**2. + y_fib**2.)
    print 'the distance from fibre to central', distan_fib_central
    area_fib_zone = sp.zeros(len(zone_radius))
    count_number = 0
    filename = utils.OUTPUTDIR + os.sep + "proportion_value.gz"
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlim = (-0.15, 0.15), ylim = (-0.15, 0.15))
    patches = []
    for x_center, y_center, radii in zip(x_fib[:4], y_fib[:4], rad_fib[:4]):
        circle = Circle((x_center, y_center), radii)
        patches.append(circle)
    circle = Circle((0., 0.), zone_radius[0])
    patches.append(circle)
    p = PatchCollection(patches, cmap = matplotlib.cm.jet, alpha = 0.4)
    ax.add_collection(p)
    ##pylab.draw()
    i_fiber_calculation = 0
    for i_circle in sp.arange(len(zone_radius)):
        print 'GOING %g th zone' %(i_circle)
        if i_circle == 0:
            for i_fib in sp.arange(len(distan_fib_central)):
                if distan_fib_central[i_fib] + rad_fib[i_fib] <= zone_radius[i_circle]:
                    area_fib_zone[i_circle] += sp.pi * sp.power(rad_fib[i_fib], 2.0)
                    i_fiber_calculation += 1
                    print 'fiber totally in the central zone', area_fib_zone[i_circle]
                elif distan_fib_central[i_fib] < zone_radius[i_circle] and (distan_fib_central[i_fib] + rad_fib[i_fib]) > zone_radius[i_circle]:
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
                    beta_r_zone = sp.arccos((2. * sp.power(zone_radius[i_circle], 2.) - sp.power(distan_two_points, 2.)) 
                                / (2.* sp.power(zone_radius[i_circle], 2.)))
                    print 'the beta_r_zone', beta_r_zone
                    if alpha_r <= sp.pi:
                        area_circ = alpha_r / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                        area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                        area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.) - \
                                    1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
                        area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
                                                area_curve
                        part_fib = area_circ - area_triangle + area_curve
                    else:
                        area_circ = (2. * sp.pi - alpha_r) / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                        area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                        area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.) - \
                            1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
                        area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ + area_triangle + \
                            area_curve
                        part_fib = area_circ + area_triangle + area_curve 
                    area_fib_zone[i_circle +1] += sp.pi * sp.power(rad_fib[i_fib], 2.) - part_fib
                    i_fiber_calculation += 1
                elif distan_fib_central[i_fib] > zone_radius[i_circle] and distan_fib_central[i_fib] \
                        - rad_fib[i_fib] < zone_radius[i_circle]:
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
                    #print 'less than half of fiber in the central zone', area_fib_zone[i_circle]
                elif distan_fib_central[i_fib] == zone_radius[i_circle]:
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
                    alpha_r = sp.arccos((2. * sp.power(rad_fib[i_fib], 2.) - sp.power(distan_two_points, 2.))
                     / (2.* sp.power(rad_fib[i_fib], 2.)))
                    beta_r_zone = sp.arccos((2. * sp.power(zone_radius[i_circle], 2.) - sp.power(distan_two_points, 2.))
                    / (2.*sp.power(zone_radius[i_circle, 2.])))
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
        else:
            for i_fib in sp.arange(len(distan_fib_central)):
                if distan_fib_central[i_fib] - rad_fib[i_fib]>= zone_radius[i_circle -1] and \
                    distan_fib_central[i_fib] + rad_fib[i_fib] <= zone_radius[i_circle]:
                    area_fib_zone[i_circle] += sp.pi * sp.power(rad_fib[i_fib], 2.)
                    i_fiber_calculation += 1
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
                    print 'the beta_r_zone', beta_r_zone
                    if alpha_r <= sp.pi:
                        area_circ = alpha_r / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                        area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                        area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.) - \
                                    1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
                        area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ - area_triangle + \
                                                area_curve
                        part_fib = area_circ - area_triangle + area_curve
                    else:
                        area_circ = (2. * sp.pi - alpha_r) / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                        area_triangle = 1./2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                        area_curve = beta_r_zone / (2. * sp.pi) * sp.pi * sp.power(zone_radius[i_circle], 2.) - \
                            1. / 2. * sp.power(zone_radius[i_circle], 2.) * sp.sin(beta_r_zone)
                        area_fib_zone[i_circle] = area_fib_zone[i_circle] + area_circ + area_triangle + \
                            area_curve
                        part_fib = area_circ + area_triangle + area_curve 
                    area_fib_zone[i_circle +1] += sp.pi * sp.power(rad_fib[i_fib], 2.) - part_fib
                    i_fiber_calculation += 1
                    #print 'the value in the zone area', area_fib_zone[i_circle]  
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
                    #print 'the value in the zone area', area_fib_zone[i_circle]
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
                    / (2.*sp.power(zone_radius[i_circle, 2.])))
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
                    #print 'the value in the zone area', area_fib_zone[i_circle]
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
    print 'the value of zone area is', zone_area
    print 'the ratio value in each zone', ratio_each
    zone_point = sp.zeros(len(zone_radius))
    for i_circle in sp.arange(len(zone_radius)):
        zone_point[i_circle] = width_zone / 2. + i_circle * width_zone
    print 'the zone central point value', zone_point
    #i_zone[:] = i_zone[:] + 1
##    dump.write({'zone_number': zone_point, 'ratio_value': ratio_each}, filename = 
##                filename, extension = '.gz')
    return (zone_point, ratio_each)