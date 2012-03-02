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
import stick.const
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
import stick.lib.utils.utils as utils
from fipy import Gmsh2D
from fipy import *
from virtlocgeom import *

def calculation_area_proportion(rad_yarn, rad_fib, x_fib, y_fib, mean_value):
    """
     divide the yarn zone into the ring zones which cut the fiber cross-section
     into three pieces
    """
    #(1)divide the domain into the ring zones which contain the virtual locations
    yarn_center = 0.
    print 'the length of rad_fib', len(rad_fib)
    rad_virtual = rad_fib[0] + mean_value
    total_circles = int((rad_yarn - rad_virtual) / (2 * rad_virtual)) + 1
    print 'the total circles number', total_circles
    #(2)divide each ring zones into three piece and calculate the center position
    #and the radius value for each zone
    delta_ring_zone = rad_virtual / 3.0 * 2.0
    total_area_zone_num = (total_circles - 1) * 3 + 2
    print 'the total number of ring zones', total_area_zone_num
    each_zone_radius = sp.zeros(total_area_zone_num, float)
    each_zone_center = sp.empty(total_area_zone_num, float)
    width_each_zone = sp.zeros(total_area_zone_num, float)
    each_zone_area = sp.zeros(total_area_zone_num, float)
    for i_circle in sp.arange(total_area_zone_num):
        if i_circle < 2:
            each_zone_radius[i_circle] = rad_virtual / 2. *(i_circle + 1)
            each_zone_center[i_circle] = yarn_center + (3. /4. * rad_virtual) *\
                                        i_circle
            width_each_zone[i_circle] = rad_virtual / 2.
        else:
            each_zone_radius[i_circle] = each_zone_radius[i_circle - 1] + \
                                        delta_ring_zone
            each_zone_center[i_circle] = each_zone_radius[i_circle - 1] + \
                                        delta_ring_zone / 2.
            width_each_zone[i_circle] = delta_ring_zone
    for i_circle in sp.arange(total_area_zone_num):
        if i_circle == 0:
            each_zone_area[i_circle] = sp.pi * sp.power(each_zone_radius[i_circle],
                                        2.)
        else:
            each_zone_area[i_circle] = sp.pi * sp.power(each_zone_radius[i_circle],
                                        2.) - sp.pi * sp.power(each_zone_radius[i_circle - 1],
                                        2.)
    print 'the area of the each zone', each_zone_area
    #(3)Calculate the distance of each fiber to the center
    distance_centra_each = sp.zeros(len(x_fib))
    for i_fib in sp.arange(len(x_fib)):
        distance_centra_each[i_fib] = sp.sqrt(x_fib[i_fib] ** 2 + y_fib[i_fib] ** 2)
    #distance_centra_each = sp.sqrt(sp.power(x_fib, 2.) + sp.power(y_fib, 2.))
    print 'x position', x_fib
    print 'distance value', distance_centra_each
    #(4)begin to calculate the cross-section area of each fiber in each ring zone
    fib_area_zone = sp.zeros(total_area_zone_num, float)
    proportion_value_zone = sp.zeros(total_area_zone_num, float)
    print 'the zone radius are in use', len(each_zone_radius)
    print 'length of radius of fiber', len(rad_fib)
    print 'each ring zone radius', each_zone_radius
    
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlim = (-1.1, 1.1), ylim = (-1.1, 1.1))
    patches_1 = []
    patches_2 = []
    for x_center, y_center, radii in zip(x_fib[:], y_fib[:], rad_fib[:]):
        circle = Circle((x_center, y_center), radii)
        patches_1.append(circle)
    for i_radii in sp.arange(len(each_zone_radius)):
        circle = Circle((0., 0.), each_zone_radius[i_radii])
        patches_2.append(circle)
    circle = Circle((0., 0.), rad_yarn)
    patches_1.append(circle)
    p_1 = PatchCollection(patches_1, facecolor = 'red', cmap = matplotlib.cm.jet, alpha = 0.4)
    p_2 = PatchCollection(patches_2, edgecolor = 'black', cmap = matplotlib.cm.jet, alpha = 0.1)    
    ax.add_collection(p_1)
    ax.add_collection(p_2)
    
    pylab.draw()
    for i_circle in sp.arange(total_area_zone_num - 1):
        if i_circle == 0:
            fib_area_zone[i_circle] = sp.pi * sp.power(each_zone_radius[i_circle],
                                    2.0)
        elif i_circle == 1:
            fib_area_zone[i_circle] = sp.pi * sp.power(rad_fib[0], 2.) - \
                                    fib_area_zone[i_circle - 1]
        else:
            for i_fib in sp.arange(len(x_fib)):                
                if distance_centra_each[i_fib] > each_zone_radius[i_circle - 1] and \
                    distance_centra_each[i_fib] < each_zone_radius[i_circle]:
                    #calculate the central part of fiber
                    #intersection point with (i-1)th ring
                    #print 'the index value for the loop in the new ring zones', i_circle
                    #print 'distance to the center',distance_centra_each[i_fib]
                    #print 'the smaller radius of ring zone', each_zone_radius[i_circle - 1]
                    solution_points_1 = intersect_circles(x_fib[i_fib],y_fib[i_fib],
                                    each_zone_radius[i_circle - 1],rad_fib[i_fib])
                    #print 'solution for the smalller circle', solution_points_1
                    x_first_1 = solution_points_1[0][0]
                    x_secon_1 = solution_points_1[1][0]
                    y_first_1 = solution_points_1[0][1]
                    y_secon_1 = solution_points_1[1][1]
                    #intersection point with ith ring
                    solution_points_2 = intersect_circles(x_fib[i_fib], y_fib[i_fib],
                                    each_zone_radius[i_circle], rad_fib[i_fib])
                    #print 'solution for the current ring', solution_points_2
                    x_first_2 = solution_points_2[0][0]
                    x_secon_2 = solution_points_2[1][0]
                    y_first_2 = solution_points_2[0][1]
                    y_secon_2 = solution_points_2[1][1]
                    distance_sq_1 = (x_first_1 - x_secon_1)**2 + \
                                    (y_first_1 - y_secon_1)**2 
                    distance_sq_2 = (x_first_2 - x_secon_2)**2 + \
                                    (y_first_2 - y_secon_2)**2 
                                    
                    alpha_r_1 = sp.arccos((2 * sp.power(rad_fib[i_fib], 2.) - distance_sq_1) /\
                                (2 * sp.power(rad_fib[i_fib], 2.)))
                    beta_ring_1 = sp.arccos((2 * sp.power(each_zone_radius[i_circle - 1], 2.) 
                                    - distance_sq_1) / (2 * sp.power(each_zone_radius
                                    [i_circle-1], 2.)))
                    piece_fib_1 = alpha_r_1 / (2 * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)\
                                - 1. / 2. * sp.sin(alpha_r_1) * sp.power(rad_fib[i_fib], 2.)
                    piece_ring_1 = beta_ring_1 / (2 * sp.pi) * sp.pi * sp.power(each_zone_radius
                                    [i_circle - 1], 2.) - 1. / 2. * sp.sin(beta_ring_1) * \
                                    sp.power(each_zone_radius[i_circle - 1], 2.)
                    total_piece_1 = piece_fib_1 + piece_ring_1

                    alpha_r_2 = sp.arccos((2 * sp.power(rad_fib[i_fib], 2.) - distance_sq_2) / 
                                (2 * sp.power(rad_fib[i_fib], 2.)))
                    beta_ring_2 = sp.arccos((2 * sp.power(each_zone_radius[i_circle], 2.) - 
                                distance_sq_2) / (2 * sp.power(each_zone_radius[i_circle], 2.)))
                    print 'the sin value of the angle in the center', sp.sin(beta_ring_2)
                    piece_fib_2 = alpha_r_2 / (2 * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)\
                                - 1. / 2. * sp.sin(alpha_r_2) * sp.power(rad_fib[i_fib], 2.)
                    piece_ring_2 = beta_ring_2 / (2 * sp.pi) * sp.pi * sp.power(each_zone_radius
                                [i_circle], 2.) - 1./ 2. * sp.sin(beta_ring_2) * sp. power(
                                each_zone_radius[i_circle], 2.)
                    #while piece_fib_2 < 0:
                    #    sys.exit()
                    total_piece_2 = -piece_ring_2 + piece_fib_2
                    area_in_zone = sp.pi * sp.power(rad_fib[i_fib], 2.) - (total_piece_1 + 
                                    total_piece_2)
                    #print 'the area in the zone', area_in_zone
                    fib_area_zone[i_circle] = fib_area_zone[i_circle] + area_in_zone
                    fib_area_zone[i_circle - 1] = fib_area_zone[i_circle - 1] + total_piece_1
                    fib_area_zone[i_circle + 1] = fib_area_zone[i_circle + 1] + total_piece_2	
        
        if fib_area_zone[i_circle] < 0:
            print 'the index of the ring zone', i_circle
##            print 'the wrong area', fib_area_zone[i_circle]
##            print i_circle
##            sys.exit()
##        print 'fiber cross-section area', fib_area_zone[i_circle]
    #check whether the calculation includes all the fibers
    total_area_in_ring = sp.sum(each_zone_area)
    #total_arae_fib = 0.
    #for i_fib in sp.arange(len(rad_fib)):
    total_area_fib = sp.sum(sp.pi * sp.power(rad_fib, 2.))
##    if np.abs(total_area_in_ring - total_area_fib) > 0.001:
##        assert "The error exceeds the limitation for area calculation in the ring zones"
    for i_circle in sp.arange(len(each_zone_area)):
        proportion_value_zone[i_circle] = fib_area_zone[i_circle] / \
                                        each_zone_area[i_circle]
    print 'the value of proportion in each ring zone', proportion_value_zone
    polyfitting = np.poly1d(np.polyfit(each_zone_center, proportion_value_zone, 
                4))
    zone_centra = sp.linspace(0.0, 0.95, 50)
    pylab.figure()
    pylab.plot(each_zone_center, proportion_value_zone, '-')
    pylab.plot(zone_centra, polyfitting(zone_centra), '.-')
    pylab.draw()
    raw_input("Enter for continue")
    return (proportion_value_zone, each_zone_area, fib_area_zone, each_zone_center)

def calculation_area_pro_shift(rad_yarn, rad_fib, x_fib, y_fib, 
                            mean_value):
    #currently the delta_r = R_{virl} / 2.
    #(1)After shifting the central, how many ring zones are in the domain
    yarn_center = 0.0
    print 'the length of rad_fib', len(rad_fib)
    rad_virtual = rad_fib[0] + mean_value
    delta_r = rad_virtual / 2.
    total_circles_shift = int((rad_yarn - delta_r) / (2 * rad_virtual))
    total_circles = int((rad_yarn - rad_virtual) / (2 * rad_virtual)) + 1
    #delta_r = rad_virtual / 2.
    #(2)divide each ring zones into three piece and calculate the center position
    #and the radius value for each zone
    delta_ring_zone = rad_virtual / 3.0 * 2.
    total_area_zone_num = (total_circles - 1) * 3 + 2
    each_zone_radius = sp.zeros(total_area_zone_num, float)
    each_zone_center = sp.empty(total_area_zone_num, float)
    width_each_zone = sp.zeros(total_area_zone_num, float)
    each_zone_area = sp.zeros(total_area_zone_num, float)
    for i_circle in sp.arange(total_area_zone_num):
        if i_circle < 2:
            each_zone_radius[i_circle] = rad_virtual / 2. *(i_circle + 1)
            each_zone_center[i_circle] = yarn_center + (3. /4. * rad_virtual) *\
                                        i_circle
            width_each_zone[i_circle] = rad_virtual / 2.
        else:
            each_zone_radius[i_circle] = each_zone_radius[i_circle - 1] + \
                                        delta_ring_zone
            each_zone_center[i_circle] = each_zone_radius[i_circle - 1] + \
                                        delta_ring_zone / 2.
            width_each_zone[i_circle] = delta_ring_zone
    print 'the radius of each zone', each_zone_radius
    for i_circle in sp.arange(total_area_zone_num):
        if i_circle == 0:
            each_zone_area[i_circle] = sp.pi * sp.power(each_zone_radius[i_circle],
                                        2.)
            print 'each zone area', each_zone_area[i_circle]
        else:
            each_zone_area[i_circle] = sp.pi * sp.power(each_zone_radius[i_circle],
                                        2.) - sp.pi * sp.power(each_zone_radius[i_circle - 1],
                                        2.)
            print 'each zone area', each_zone_area[i_circle]
    #(3)Calculate the distance of each fiber to the center
    distance_centra_each = sp.sqrt(sp.power(x_fib, 2.0) + sp.power(y_fib, 2.0))
    
    #(4)begin to calculate the cross-section area of each fiber in each ring zone
    fib_area_zone = sp.zeros(total_area_zone_num, float)
    #proportion_value_zone = sp.zeros(total_area_zone_num, float)
    #(5)calculate the determining radius for angle 
    determining_rad = sp.zeros(total_circles - 1, float)
    print 'the length of the determining_rad', len(determining_rad)
##    for i_shift in sp.arange(total_circles - 1):
##        if i_shift == 0:
##            determining_rad[i_shift] = sp.sqrt(sp.power(each_zone_radius[i_shift + 2], 
##                                    2.) - sp.power(3. / 2. * rad_virtual, 2.))
##                                    
##        else:
##            print 'i_shift value', i_shift
##            determining_rad[i_shift] = sp.sqrt(sp.power(each_zone_radius[i_shift + 2 * 
##                                    (i_shift + 1)], 2.) - sp.power((3. / 2. + 2. * 
##                                    i_shift) * rad_virtual, 2.))
    i_determining = 0
    for i_circle in sp.arange(total_area_zone_num): 
                      
        if i_circle > 0:
            print 'the circle value', i_circle 
            for i_fib in sp.arange(len(x_fib)):
##                if distance_centra_each[i_fib] - rad_fib[i_fib] > each_zone_radius\
##                    [i_circle -1] and distance_centra_each[i_fib] > each_zone_radius\
##                    [i_circle] and distance_centra_each[i_fib] < each_zone_radius\
##                    [i_circle + 1]:
##                    alpha_r = a
##                elif distance_centra_each[i_fib] > each_zone_radius[i_circle] and \
##                    distance_centra_each[i_fib + 1]:
##                    aph
##                    
##                elif distance_centra_each[i_fib] + rad_fib[i_fib] > each_zone_radius\
##                    [i_circle + 1] and distance_centra_each[i_fib] < each_zone_radius\
##                    [i_circle + 2]:
                #calculate the area proportional value for 'shift'
                if distance_centra_each[i_fib] > each_zone_radius[i_circle] and \
                    distance_centra_each[i_fib] - rad_fib[i_fib] < each_zone_radius[i_circle]:
                    solution_points_1 = intersect_circles(x_fib[i_fib], y_fib[i_fib],
                                each_zone_radius[i_circle], rad_fib[i_fib])
                    x_first = solution_points_1[0][0]
                    x_secon = solution_points_1[1][0]
                    y_first = solution_points_1[0][1]
                    y_secon = solution_points_1[1][1]
                    distance_sq = (x_first - x_secon) ** 2 + (y_first - y_secon) ** 2
                    print 'the distance_sq', distance_sq
                    alpha_r = sp.arccos((2. * sp.power(rad_fib[i_fib], 2.) - distance_sq) / 
                            (2. * sp.power(rad_fib[i_fib], 2.)))
                    print 'alpha_r', alpha_r
                    beta_centra = sp.arccos ((2. * sp.power(each_zone_radius[i_circle], 2.) - 
                            distance_sq) / (2. * sp.power(each_zone_radius[i_circle], 2.)))
                    print 'beta_centra', beta_centra
                    piece_fib = alpha_r / (2. * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)
                    piece_ring = beta_centra / (2. * sp.pi) * sp.pi * sp.power(each_zone_radius
                                [i_circle], 2.)
                    triangle_fib = 1. / 2. * sp.sin(alpha_r) * sp.power(rad_fib[i_fib], 2.)
                    triangle_ring = 1. / 2. * sp.sin(beta_centra) * sp.power(each_zone_radius
                                    [i_circle], 2.)
                    area_fib = piece_fib - triangle_fib
                    area_ring = piece_ring - triangle_ring
                    total_piece = area_fib + area_ring
                    fib_area_zone[i_circle] = fib_area_zone[i_circle] + total_piece
                elif distance_centra_each[i_fib] < each_zone_radius[i_circle] and \
                    distance_centra_each[i_fib] > each_zone_radius[i_circle - 1]:
                    #print 'previous determining value is', i_determining
                    i_determining += 1
                    #print 'determining value is', i_determining
                    solution_points_1 = intersect_circles(x_fib[i_fib], y_fib[i_fib],
                                each_zone_radius[i_circle -1], rad_fib[i_fib])
                    
                    solution_points_2 = intersect_circles(x_fib[i_fib], y_fib[i_fib], 
                                each_zone_radius[i_circle], rad_fib[i_fib])
                    #intersection point with (i-1)th ring
                    x_first_1 = solution_points_1[0][0]
                    x_secon_1 = solution_points_1[1][0]
                    y_first_1 = solution_points_1[0][1]
                    y_secon_1 = solution_points_1[1][1]
                    #intersection point with ith ring
                    x_first_2 = solution_points_2[0][0]
                    x_secon_2 = solution_points_2[1][0]
                    y_first_2 = solution_points_2[0][1]
                    y_secon_2 = solution_points_2[1][1]
                    
                    distance_sq_1 = (x_first_1 - x_secon_1) ** 2 + \
                                (y_first_1 - y_secon_1) ** 2 
                    distance_sq_2 = (x_first_2 - x_secon_2) ** 2 + \
                                (y_first_2 - y_secon_2) ** 2
                    alpha_r_1 = sp.arccos((2 * sp.power(rad_fib[i_fib], 2.) - distance_sq_1) /\
                            (2 * sp.power(rad_fib[i_fib], 2.)))
                    beta_ring_1 = sp.arccos((2 * sp.power(each_zone_radius[i_circle - 1], 2.) 
                                - distance_sq_1) / (2 * sp.power(each_zone_radius
                                [i_circle - 1], 2.)))
                    piece_fib_1 = alpha_r_1 / (2 * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)\
                            - 1. / 2. * sp.sin(alpha_r_1) * sp.power(rad_fib[i_fib], 2.)
                    piece_ring_1 = beta_ring_1 / (2 * sp.pi) * sp.pi * sp.power(each_zone_radius
                                [i_circle - 1], 2.) - 1. / 2. * sp.sin(beta_ring_1) * \
                                sp.power(each_zone_radius[i_circle - 1], 2.)
                    total_piece_1 = piece_fib_1 + piece_ring_1
                    alpha_r_2 = sp.arccos((2 * sp.power(rad_fib[i_fib], 2.) 
                                - distance_sq_2) / (2 * sp.power(rad_fib[i_fib]
                                , 2.)))
                    beta_ring_2 = (2 * sp.power(each_zone_radius[i_circle], 
                                2.) - distance_sq_2) / (2. * sp.power(each_zone_radius[i_circle],
                                2.))
                    #if rad_fib[i_fib] <= determining_rad[i_circle - 2 * i_determining]:
                    piece_fib_2 = alpha_r_2 / (2 * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)\
                        - 1. / 2. * sp.sin(alpha_r_2) * sp.power(rad_fib[i_fib], 2.)
                    piece_ring_2 = beta_ring_1 / (2 * sp.pi) * sp.pi * sp.power(each_zone_radius
                            [i_circle], 2.) - 1. / 2. * sp.sin(beta_ring_2) * \
                            sp.power(each_zone_radius[i_circle], 2.)
                    total_piece_2 = piece_fib_2 + piece_ring_2
                    #elif rad_fib[i_fib] > determining_rad[i_circle - 2 * i_determining]:
##                        piece_fib_2 = alpha_r_2 / (2. * sp.pi) * sp.pi * sp.power(
##                                    rad_fib[i_fib], 2.)
##                        piece_ring_2 = beta_ring_2 / (2. * sp.pi) * sp.pi * \
##                                    sp.power(each_zone_radius[i_circle], 2.) - 1./ 2. * \
##                                    sp.sin(beta_ring_2) * sp.power(each_zone_radius[i_circle],
##                                    2.) - 1. / 2. * sp.sin(alpha_r_2) * sp.power(rad_fib[i_fib], 2.)
##                        total_piece_2 = piece_fib_2 - piece_ring_2
                    left_part = sp.pi * sp.power(rad_fib[i_fib], 2.) - (total_piece_2 + 
                                total_piece_1)
                    fib_area_zone[i_circle] = fib_area_zone[i_circle] + left_part
                elif distance_centra_each[i_fib] + rad_fib[i_fib] > each_zone_radius[i_circle] and \
                    distance_centra_each[i_fib] < each_zone_radius[i_circle - 1]:
                    solution_points_1 = intersect_circles(x_fib[i_fib], y_fib[i_fib],
                                each_zone_radius[i_circle -1], rad_fib[i_fib])
                    solution_points_2 = intersect_circles(x_fib[i_fib], y_fib[i_fib], 
                                each_zone_radius[i_circle], rad_fib[i_fib])
                    #intersection point with (i-1)th ring
                    x_first_1 = solution_points_1[0][0]
                    x_secon_1 = solution_points_1[1][0]
                    y_first_1 = solution_points_1[0][1]
                    y_secon_1 = solution_points_1[1][1]
                    #intersection point with ith ring
                    x_first_2 = solution_points_2[0][0]
                    x_secon_2 = solution_points_2[1][0]
                    y_first_2 = solution_points_2[0][1]
                    y_secon_2 = solution_points_2[1][1]
                    #if x_first_2.imag == 0:
                    distance_sq_1 = (x_first_1 - x_secon_1) ** 2 + \
                                (y_first_1 - y_secon_1) ** 2
                    distance_sq_2 = (x_first_2 - x_secon_2) ** 2 + \
                                (y_first_2 - y_secon_2) ** 2 
                    alpha_r_2 = sp.arccos((2 * sp.power(rad_fib[i_fib], 2.) - distance_sq_2) /\
                            (2 * sp.power(rad_fib[i_fib], 2.)))
                    beta_ring_2 = sp.arccos((2 * sp.power(each_zone_radius[i_circle], 2.) 
                                - distance_sq_2) / (2 * sp.power(each_zone_radius
                                [i_circle], 2.)))
                    piece_fib_2 = alpha_r_2 / (2 * sp.pi) * sp.pi * sp.power(rad_fib[i_fib], 2.)\
                            - 1. / 2. * sp.sin(alpha_r_2) * sp.power(rad_fib[i_fib], 2.)
                    piece_ring_2 = beta_ring_2 / (2 * sp.pi) * sp.pi * sp.power(each_zone_radius
                                [i_circle], 2.) - 1. / 2. * sp.sin(beta_ring_2) * \
                                sp.power(each_zone_radius[i_circle], 2.)
                    total_piece_2 = piece_fib_2 - piece_ring_2
                    alpha_r_1 = sp.arccos((2 * sp.power(rad_fib[i_fib], 
                                2.) - distance_sq_1) /(2 * sp.power(rad_fib[i_fib], 
                                2.)))
                    beta_ring_1 = sp.arccos((2 * sp.power(each_zone_radius[i_circle - 1], 
                                2.) - distance_sq_1) / (2 * sp.power(each_zone_radius
                                [i_circle - 1], 2.)))
##                        if rad_fib[i_fib] <= determining_rad[i_circle - 2 * i_determining - 1]:
##                            piece_fib_1 = (2. * sp.pi - alpha_r_1) / (2. * sp.pi) * \
##                                        sp.pi * sp.power(rad_fib[i_fib], 2.) + \
##                                        1. / 2. * sp.sin(alpha_r_1) * sp.power(
##                                        rad_fib[i_fib], 2.)
##                            piece_ring_1 = beta_ring_1 / (2. * sp.pi) * sp.pi * \
##                                        sp.power(each_zone_radius[i_circle - 1], 2.) - \
##                                        1. / 2. * sp.sin(beta_ring_1) * sp.power(
##                                        each_zone_radius[i_circle - 1])
##                            total_piece_1 = piece_fib_1 + piece_ring_1
##                        elif rad_fib[i_fib] > determining_rad[i_circle - 2 * i_determining - 1]:
                    piece_fib_1 = alpha_r_1 / (2. * sp.pi) * sp.pi * \
                                sp.power(rad_fib[i_fib], 2.) 
                    piece_ring_1 = beta_ring_1 / (2. * sp.pi) * sp.pi * \
                                sp.power(each_zone_radius[i_circle - 1], 2.) - \
                                1. / 2. * sp.sin(beta_ring_1) * sp.power(
                                each_zone_radius[i_circle - 1], 2.) - 1. / 2. * \
                                sp.sin(alpha_r_1) * sp.power(rad_fib[i_fib], 2.)
                    total_piece_1 = piece_fib_1 - piece_ring_1
                    left_part = sp.pi * sp.power(rad_fib[i_fib], 2.) - (total_piece_1 + 
                                total_piece_2)
                    fib_area_zone[i_circle + 1] = fib_area_zone[i_circle + 1] + \
                                                total_piece_2
##                    else:
##                        alpha_r_1 = sp.arccos((2. * sp.power(rad_fib[i_fib], 2.) 
##                                    - distance_sq_1) / (2. * sp.power(rad_fib[i_fib], 2.)))
##                        beta_ring_1 = sp.arccos((2. * sp.power(each_zone_radius[i_circle - 1], 
##                                    2.) - distance_sq_1) / (2. * 
##                                    sp.power(each_zone_radius[i_circle - 1], 2.)) )
##                        if rad_fib[i_fib] <= determining_rad[i_circle - 2 * i_determining - 1]:
##                            piece_fib_1 = alpha_r_1 / (2. * sp.pi) * sp.pi * \
##                                        sp.power(rad_fib[i_fib], 2.) - \
##                                        1. / 2. * sp.sin(alpha_r_1) * \
##                                        sp.power(rad_fib[i_fib], 2.)
##                            piece_ring_1 = beta_ring_1 / (2. * sp.pi) * sp.pi * \
##                                        sp.power(each_zone_radius[i_circle - 1], 2.) \
##                                        - 1./ 2. * sp.sin(beta_ring_1) * sp.power(
##                                        each_zone_radius[i_circle - 1], 2.)
##                            left_part = piece_fib_1 - piece_ring_1
##                        else:
##                            piece_fib_1 = (2. * sp.pi - alpha_r_1) / (2. * sp.pi) * \
##                                        sp.pi * sp.power(rad_fib[i_fib], 2.)
##                            piece_ring_1 = beta_ring_1 / (2. * sp.pi) * sp.pi * \
##                                        sp.power(each_zone_radius[i_circle - 1], 2.) - \
##                                        1. / 2. * sp.sin(beta_ring_1) * sp. power(
##                                        each_zone_radius[i_circle - 1], 2.) - \
##                                        1./ 2. * sp.sin(alpha_r_1) * sp.power(
##                                        rad_fib[i_fib], 2.)
##                            left_part = piece_fib_1 - piece_ring_1
                    fib_area_zone[i_circle] = fib_area_zone[i_circle] + left_part
    sum_each_zone = sp.sum(each_zone_area)
    area_each_fib = sp.pi * sp.power(rad_fib[:], 2.)
    sum_fib = sp.sum(area_each_fib)
    print 'the total value of fib area in each sub-ring zone', sum_each_zone
    print 'the total value of the fib cross-section', sum_fib
    proportion_each_zone = fib_area_zone / each_zone_area
    print 'the proportional value for each sub-ring zone', proportion_each_zone
    return (proportion_each_zone, fib_area_zone, each_zone_center, each_zone_area)
##                elif distance_centra_each[i_fib] + rad_fib[i_fib] > each_zone_radius[i_circle - 1] and \
##                    distance_centra_each[i_fib] + rad_fib[i_fib] <= each_zone_radius[i_circle] + \
##                    1. / 6. * rad_fib[i_fib]:
                    