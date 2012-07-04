"""
Example showing how a fiber-yarn layout can be constructed from a given
area probability function
"""

import os, sys
import numpy as np
import scipy as sp
import pylab
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon, Ellipse
from matplotlib.collections import PatchCollection
from yarn2d.distribution_method import *


# We start with inifile settings in 
ini_yarn = u"""
[general]
method = 'FVM'
submethod = 'fipy'
read = False
verbose = False
[domain]
fiberlayout_method = 'virtlocoverlap'
distribute_fiber = 'integral'
cellsize_centre = 1.0e-1
cellsize_fiber = 2.50e-2
yarnradius = 1.0
[fiber]
number_type = 2
number_fiber = 720
blend = [16.7, 83.3]
eps_value = 0.001
fiber_config = ['tmpfiber1.ini', 'tmpfiber2.ini']
prob_area = '[lambda r: 1.*(-7.278*r**6 + 6.510* r ** 5 + 5.727 * r ** 4 - 6.763 * r ** 3 + 1.137*r**2 + 0.405*r + 0.262), lambda r: (12.423*r**6 - 12.339*r**5 - 10.*r**4 + 115.424*r**3 - 52.966*r**2 + 7.109*r + 0.765)]'
[plot]
maxval = 0.0005
plotevery = 10
[writeout]
writeevery = 100
"""

ini_fiber1 = u"""
[general]
method = 'FVM'
submethod = 'odew'
read = False
verbose = True
[fiber]
radius_pure_fiber = 0.0543
form = 'circle'
nrlayers = 2
mean_deviation = 0.0
[fiberlayer_0]
thickness = 0.00085
[fiberlayer_1]
thickness = 0.00085
[plot]
plotevery = 1
"""
ini_fiber2 = u"""
[general]
method = 'FVM'
submethod = 'odew'
read = False
verbose = False
[fiber]
radius_pure_fiber = 0.0181
form = 'ellipse'
eccentricity = 0.7
nrlayers = 2
mean_deviation = 0.0
[fiberlayer_0]
thickness = 0.001
[fiberlayer_1]
thickness = 0.001
[plot]
plotevery = 10
"""

#read the statistic data from a real blended yarn from which the probability 
#distribution in the ini file has been derived
path = os.path.dirname(os.path.abspath(__file__))
data_polyester = np.loadtxt(os.path.join(path,'../fiber_polyester.csv'))
data_cotton = np.loadtxt(os.path.join(path,'../fiber_cotton.csv'))

x_position_real_fiber = []
y_position_real_fiber = []
x_position_polyester = []
y_position_polyester = []
x_position_cotton = []
y_position_cotton = []

radius_real_fiber = []
radius_polyester = []
radius_cotton = []
for i_polyester in sp.arange(len(data_polyester)):
    x_position_real_fiber.append(data_polyester[i_polyester][0])
    y_position_real_fiber.append(data_polyester[i_polyester][1])
    x_position_polyester.append(data_polyester[i_polyester][0])
    y_position_polyester.append(data_polyester[i_polyester][1])
    radius_real_fiber.append(data_polyester[i_polyester][2])
    radius_polyester.append(data_polyester[i_polyester][2])
for i_cotton in sp.arange(len(data_cotton)):
    x_position_real_fiber.append(data_cotton[i_cotton][0])
    y_position_real_fiber.append(data_cotton[i_cotton][1])
    x_position_cotton.append(data_cotton[i_cotton][0])
    y_position_cotton.append(data_cotton[i_cotton][1])
    radius_real_fiber.append(data_cotton[i_cotton][2])
    radius_cotton.append(data_cotton[i_cotton][2])
x_position_real_fiber = np.array(x_position_real_fiber)
y_position_real_fiber = np.array(y_position_real_fiber)
x_position_polyester = np.array(x_position_polyester)
y_position_polyester = np.array(y_position_polyester)
x_position_cotton = np.array(x_position_cotton)
y_position_cotton = np.array(y_position_cotton)
radius_real_fiber = np.array(radius_real_fiber)

radius_polyester = np.array(radius_polyester)
sum_surface_poly = sp.sum(sp.pi * 2. * radius_polyester)
sum_area_poly = sp.sum(sp.pi * sp.power(radius_polyester, 2.))
print 'the total area of polyester cross-section in the real yarn', sum_area_poly
print 'the surface value in 2D dimension for polyester', sum_surface_poly
radius_cotton = np.array(radius_cotton)
sum_surface_cotton = sp.sum(sp.pi * 2. * radius_cotton)
sum_area_cotton = sp.sum(sp.pi * sp.power(radius_cotton, 2.))
print 'the total area of cotton cross-section in the real yarn', sum_area_cotton
print 'the surface value in 2D dimension for cotton', sum_surface_cotton
raw_input("finish calculating the sum of cross-section area value  for each kind of fiber")

fiber_kind = np.zeros(len(radius_real_fiber), int)
fiber_kind[len(data_polyester):] = 1

from yarn2d.fiber_layout import plot_yarn
from yarn2d.arearatioprobability import (calculate_proportion, 
            plot_ratio_function, compare_relative_error)

plot_yarn(x_position_real_fiber, y_position_real_fiber, radius_real_fiber, 
          fiber_kind, title='Real fiber-yarn layout')

prob_real = calculate_proportion(1.0, radius_real_fiber, x_position_real_fiber, y_position_real_fiber,  
           nrzones=5.)

prob_poly = calculate_proportion(1.0, radius_polyester, x_position_polyester, y_position_polyester,
            nrzones = 5.)

prob_cotton = calculate_proportion(1.0, radius_cotton, x_position_cotton, y_position_cotton,
            nrzones = 5.)

zone_position_real = sp.ones(len(prob_real[0]) + 1, float)
zone_position_real[0] = 0.
zone_position_real[1:] = prob_real[0][:]
zone_position_real[-1] = 1.

zone_position_polyester = sp.zeros(len(prob_real[0]) + 2, float)
zone_position_polyester[1:-1] = prob_real[0][:]
print 'the coordinate value is', prob_real[0][:]
zone_position_polyester[-1] = 1.
zone_position_polyester[-2] = 1. - zone_position_polyester[1] / 2.

ratio_poly = sp.zeros(len(zone_position_polyester), float)
ratio_cotton = sp.zeros(len(prob_cotton[-1])+1 , float)
ratio_poly[1:-1] = prob_poly[-1][:]
ratio_cotton[:-1] = prob_cotton[-1][:]
#print len(zone_position_real), len(ratio_poly), len(x_position_polyester), len(x_position_cotton)
poly_polyester = np.poly1d(np.polyfit(zone_position_polyester, ratio_poly, 4))
poly_cotton = np.poly1d(np.polyfit(zone_position_real, ratio_cotton, 4))
draw_real = sp.linspace(0., 1.0, 50)
print np.polyfit(zone_position_polyester, ratio_poly, 4)
print np.polyfit(zone_position_real, ratio_cotton, 4)
raw_input("record the coefficients of polynomial equations")

pylab.figure()
pylab.subplot(121)
pylab.plot(zone_position_polyester, ratio_poly, 'o')
pylab.plot(draw_real, poly_polyester(draw_real), '--')
pylab.xlim(0., 1.05)
pylab.ylim(0., 0.8)
pylab.subplot(122)
pylab.plot(zone_position_real, ratio_cotton, '*')
pylab.plot(draw_real, poly_cotton(draw_real), '-')
pylab.xlim(0., 1.05)
pylab.ylim(0., 0.8)
pylab.axis()
pylab.show()
raw_input("check the figure")
#set up a yarn computation
from fiber.config import FiberConfigManager
from yarn.config import YarnConfigManager
from lib.utils.utils import set_outputdir
cfgf1 = FiberConfigManager.get_instance('tmpfiber1.ini', realdatastr=ini_fiber1)
cfgf2 = FiberConfigManager.get_instance('tmpfiber2.ini', realdatastr=ini_fiber2)
#cfgf3 = FiberConfigManager.get_instance('tmpfiber3.ini', realdatastr=ini_fiber3)
cfg = YarnConfigManager.get_instance('tmpyarn.ini', realdatastr=ini_yarn)
#create outputdir if not existing
if not os.path.isdir('temp'):
    os.mkdir('temp')
set_outputdir('temp')
#create 10 2D grids for statistics
from yarn2d.yarn2dgrid import Yarn2dGrid
grid = Yarn2dGrid(cfg)
ouroptions = {
                'x_central' : grid.x_central,
                'y_central' : grid.y_central,
                'number_fiber' : grid.number_fiber,
                'number_fiber_blend' : grid.number_fiber_blend,
                'radius_fiber' : grid.radius_fiber,
                'radius_yarn' : grid.radius_yarn,
                'theta_value' : grid.theta_value,
                'beta_value' : grid.beta_value,
                'mean_deviation': grid.mean_deviation,
                'prob_area': grid.prob_area,
                'radius_first_center': cfg.get(
                                    'domain.radius_first_center_virtloc'),
                'distribute_fiber': cfg.get('domain.distribute_fiber')
                }
from yarn2d.fiber_layout import virtlocoverlaplayout
#After generating n times of iteration, plot prob func result from the average ratio value
coefficient_function = [1.0, 0.0065]
iteration = 1
each_time_ratio = []
relative_error_each = []
relative_error_alpha = []
relative_error_alpha_1 = []
for i in range(iteration):
    x_position, y_position, x_position_alpha,y_position_alpha, x_position_alpha_1, y_position_alpha_1, all_radius_fibers, \
                    fiber_kind,type_fiber, ratio_each_kind, \
                    zone_position_each, ratio_each_alpha = virtlocoverlaplayout(ouroptions)
##    plot_yarn(x_position, y_position, all_radius_fibers, 
##              fiber_kind,) #title='Realization %d' % i)
    #raw_input("capture")
    ratio_each = [0] * type_fiber
    ratio_each_alpha = [0] * type_fiber
    ratio_each_alpha_1 = [0] * type_fiber
    zone_position = [0] * type_fiber
    #zone_position = [0]
    for i_type in sp.arange(type_fiber):
        x_position_cal = []
        y_position_cal = []
        x_position_alpha_cal = []
        y_position_alpha_cal = []
        x_position_alpha_1_cal = []
        y_position_alpha_1_cal = []
        radius_each = []
        for i_fiber in sp.arange(len(fiber_kind)):
            if fiber_kind[i_fiber] == i_type:
                x_position_cal.append(x_position[i_fiber])
                y_position_cal.append(y_position[i_fiber])
                x_position_alpha_cal.append(x_position_alpha[i_fiber])
                y_position_alpha_cal.append(y_position_alpha[i_fiber])
                x_position_alpha_1_cal.append(x_position_alpha_1[i_fiber])
                y_position_alpha_1_cal.append(y_position_alpha_1[i_fiber])
                radius_each.append(all_radius_fibers[i_fiber])
        x_position_cal = np.array(x_position_cal)
        y_position_cal = np.array(y_position_cal)
        x_position_alpha_cal = np.array(x_position_alpha_cal)
        y_position_alpha_cal = np.array(y_position_alpha_cal)
        x_position_alpha_1_cal = np.array(x_position_alpha_1_cal)
        y_position_alpha_1_cal = np.array(y_position_alpha_1_cal)
        radius_each = np.array(radius_each)
        probs = calculate_proportion(grid.radius_yarn, all_radius_fibers, 
                x_position_cal, y_position_cal, nrzones=5)
        probs_alpha = calculate_proportion(grid.radius_yarn, all_radius_fibers,
                    x_position_alpha_cal, y_position_alpha_cal, nrzones = 5)
        probs_alpha_1 = calculate_proportion(grid.radius_yarn, all_radius_fibers,
                    x_position_alpha_1_cal, y_position_alpha_1_cal, nrzones = 5)
##        print "when alpha = 1", probs_alpha_1
##        raw_input("check the value")
        nrzones = 5
        zone_position[i_type] = sp.zeros(len(zone_position_each) + 1, float)
        zone_position[i_type][0] = 0.
        zone_position[i_type][1:] = zone_position_each[:]
        zone_position[i_type][-1] = 1.
#calculate the relative error in the calculation
        value_from_function = grid.prob_area[i_type]
        #print 'the ratio value from the calculation', ratio_compare
        #raw_input('to generate the relative error')
        if i_type == 0:
            zone_position[i_type] = sp.zeros(len(probs[0]) +  2, float)
            #print 'check the length', len(zone_position[i_type][1:-1]), len(probs[0][:])
            zone_position[i_type][1:-1] = probs[0][:]#probs[0][:]
            zone_position[i_type][-2] = 1. - grid.radius_yarn / ((nrzones-1)*2 + 1)/2.
            zone_position[i_type][-1] = 1.
            ratio_each[i_type] = sp.zeros(len(zone_position[i_type]), float)
            ratio_each[i_type][1:-1] = probs[-1][:]
            ratio_each_alpha[i_type] = sp.zeros(len(zone_position[i_type]), float)
            ratio_each_alpha[i_type][1:-1] = probs_alpha[-1][:]
            ratio_each_alpha_1[i_type] = sp.zeros(len(zone_position[i_type]), float)
            ratio_each_alpha_1[i_type][1:-1] = probs_alpha_1[-1][:]
        else:
            zone_position[i_type] =sp.zeros(len(probs[0]) + 1, float)
##            zone_position[i_type][0] = 0.
            zone_position[i_type][1:] = probs[0][:] + + grid.radius_yarn / ((nrzones-1)*2 + 1) / 2.
            zone_position[i_type][-3] = zone_position[i_type][-3] + grid.radius_yarn / ((nrzones-1)*2 + 1)
            zone_position[i_type][-2] = 1. - grid.radius_yarn / ((nrzones-1)*2 + 1)
            zone_position[i_type][-1] = 1.
            ratio_each[i_type] = sp.ones(len(zone_position[i_type]), float)
            ratio_each[i_type][:-1] = probs[-1][:]
            ratio_each[i_type][-1] = 0.
            ratio_each_alpha[i_type] = sp.ones(len(zone_position[i_type]), float)
            ratio_each_alpha[i_type][:-1] = probs_alpha[-1][:]
            ratio_each_alpha[i_type][-1] = 0.
            ratio_each_alpha_1[i_type] = sp.ones(len(zone_position[i_type]), float)
            ratio_each_alpha_1[i_type][:-1] = probs_alpha_1[-1][:]
            ratio_each_alpha_1[i_type][-1] = 0.
        #the ratio value on the yarn's edge
        if i_type == 0:
            ratio_compare = sp.ones(len(zone_position[i_type]))
        else:
            ratio_compare = value_from_function(zone_position[i_type]) / \
                            coefficient_function[i_type]
        #print 'the ratio value for eacn kind of fiber', ratio_each_kind[i_type][:]
        if i_type == 0:
            relative_error_each_kind = sp.power(abs(ratio_compare[1:-1] - 
                                        ratio_each[i_type][1:-1])\
                                                / ratio_compare[1:-1], 2.)
            relative_error_kind_alpha = sp.power(abs(ratio_compare[1:-1] - 
                                        ratio_each_alpha[i_type][1:-1])\
                                                /ratio_compare[1:-1], 2.)
            relative_error_kind_alpha_1 = sp.power(abs(ratio_compare[1:-1] - 
                                        ratio_each_alpha_1[i_type][1:-1])\
                                                /ratio_compare[1:-1], 2.)
        else:
            relative_error_each_kind = sp.power(abs(ratio_compare[:-1] - 
                                        ratio_each[i_type][:-1])\
                                                / ratio_compare[:-1], 2.)
            relative_error_kind_alpha = sp.power(abs(ratio_compare[:-1] - 
                                        ratio_each_alpha[i_type][:-1])\
                                                /ratio_compare[:-1], 2.)
            relative_error_kind_alpha_1 = sp.power(abs(ratio_compare[:-1] - 
                                        ratio_each_alpha_1[i_type][:-1])\
                                                /ratio_compare[:-1], 2.)

        #raw_input('check the ratio value from the calculation')
        each_time_ratio.append(ratio_each[i_type])
        relative_error_each.append(relative_error_each_kind)
        relative_error_alpha.append(relative_error_kind_alpha)
        relative_error_alpha_1.append(relative_error_kind_alpha_1)
relative_error_each = sp.array(relative_error_each)
for i_error in relative_error_each:
    print "the error of two kinds of fibers:", i_error
relative_error_alpha = sp.array(relative_error_alpha)
relative_error_alpha_1 = sp.array(relative_error_alpha_1)
mean_value_each = [0] * type_fiber
mean_value_alpha = [0] * type_fiber
mean_value_alpha_1 = [0] * type_fiber
for i_type in sp.arange(type_fiber):
    each_kind_fiber = []
    each_kind_fiber_alpha = []
    each_kind_fiber_alpha_1 = []
    for i_iteration in sp.arange(iteration):
        index = i_type + 2 * i_iteration
        each_kind_fiber.append(relative_error_each[index])
        each_kind_fiber_alpha.append(relative_error_alpha[index])
        each_kind_fiber_alpha_1.append(relative_error_alpha_1[index])
    each_kind_fiber = sp.array(each_kind_fiber)
    each_kind_fiber_alpha = sp.array(each_kind_fiber_alpha)
    each_kind_fiber_alpha_1 = sp.array(each_kind_fiber_alpha_1)
    mean_value_each[i_type] = np.mean(each_kind_fiber, axis = 0)
    mean_value_alpha[i_type] = np.mean(each_kind_fiber_alpha, axis = 0)
    mean_value_alpha_1[i_type] = np.mean(each_kind_fiber_alpha_1, axis = 0)
##print 'relative error for each kind of fibers', mean_value_each, mean_value_alpha
##raw_input("enter")        
###sum_error_each_pos = sp.zeros(len(relative_error_each[0]))
###sum_error_alpha_pos = sp.zeros(len(relative_error_alpha[0]))
##for i_error in sp.arange(len(relative_error_each)):
##    for i_suberror in sp.arange(len(relative_error_each[i_error])):
##        sum_error_each_pos[i_suberror] += relative_error_each[i_error][i_suberror]
##        sum_error_alpha_pos[i_suberror] += relative_error_alpha[i_error][i_suberror]
##average_relative_error_each = sum_error_each_pos / iteration
##average_relative_error_alpha = np.mean(sum_error_alpha_pos) / iteration
##mean_value_each = np.mean(relative_error_each, axis = 0)
##mean_value_alpha = np.mean(relative_error_alpha, axis = 0)
##print 'relative error for proportionally moving', average_relative_error_each, \
##        mean_value_each
##print 'relative error for fixed value of moving', average_relative_error_alpha, \
##        mean_value_alpha
##==============================================================================
##each_time_ratio = []
##relative_error_each = []
##relative_error_alpha = []
##    
##for i in range(iteration):
##    x_position, y_position, all_radius_fibers, \
##                    fiber_kind,type_fiber, ratio_each_kind, \
##                    zone_position_each, ratio_each_alpha = virtlocoverlaplayout(ouroptions)
##    plot_yarn(x_position, y_position, all_radius_fibers, 
##              fiber_kind, title='Realization %d' % i)
##    #raw_input("capture")
##    ratio_each = [0] * type_fiber
##    zone_position = [0] * type_fiber
##    #zone_position = [0]
##    for i_type in sp.arange(type_fiber):
##        x_position_cal = []
##        y_position_cal = []
##        radius_each = []
##        for i_fiber in sp.arange(len(fiber_kind)):
##            if fiber_kind[i_fiber] == i_type:
##                x_position_cal.append(x_position[i_fiber])
##                y_position_cal.append(y_position[i_fiber])
##                radius_each.append(all_radius_fibers[i_fiber])
##        x_position_cal = np.array(x_position_cal)
##        y_position_cal = np.array(y_position_cal)
##        radius_each = np.array(radius_each)
##        probs = calculate_proportion(grid.radius_yarn, all_radius_fibers, 
##                x_position_cal, y_position_cal, nrzones=5)
##        nrzones = 5
##        zone_position[i_type] = sp.zeros(len(zone_position_each) + 1, float)
##        zone_position[i_type][0] = 0.
##        zone_position[i_type][1:] = zone_position_each[:]
##        zone_position[i_type][-1] = 1.
###calculate the relative error in the calculation
##        value_from_function = grid.prob_area[i_type]
##        ratio_compare = value_from_function(zone_position[i_type]) / coefficient_function[i_type]
##        print 'the ratio value from the calculation', ratio_compare
##        #raw_input('to generate the relative error')
##        
##        if i_type == 0:
##            zone_position[i_type] = sp.zeros(len(probs[0]) +  2, float)
##            print 'check the length', len(zone_position[i_type][1:-1]), len(probs[0][:])
##            zone_position[i_type][1:-1] = probs[0][:]
##            zone_position[i_type][-2] = 1. - grid.radius_yarn / ((nrzones-1)*2 + 1)
##            zone_position[i_type][-1] = 1.
##            ratio_each[i_type] = sp.zeros(len(probs[0]) + 1, float)
##            ratio_each[i_type][1:] = probs[-1][:]
##        else:
##            zone_position[i_type] =sp.zeros(len(probs[0]) + 1, float)
##            zone_position[i_type][0] = 0.
##            #zone_position[i_type][1:] = probs[0][:]
##            zone_position[i_type][-2] = 1. - grid.radius_yarn / ((nrzones-1)*2 + 1)
##            zone_position[i_type][-1] = 1.
##
##        #ratio_each[i_type] = sp.ones(len(probs[-1]), float)
##        #ratio_each[i_type][:] = probs[-1][:]
##        #the ratio value on the yarn's edge
##        
##        print 'the ratio value for eacn kind of fiber', ratio_each_kind[i_type][:]
##        relative_error_each_kind = sp.power(abs(ratio_compare[:-1] - ratio_each_kind[i_type][:])\
##                                            / ratio_compare[:-1], 2.)
##        relative_error_kind_alpha = sp.power(abs(ratio_compare[:-1] - ratio_each_alpha[i_type][:])\
##                                            /ratio_compare[:-1], 2.)
##        #raw_input('check the ratio value from the calculation')
##        each_time_ratio.append(ratio_each_kind[i_type])
##        relative_error_each.append(relative_error_each_kind)
##        relative_error_alpha.append(relative_error_kind_alpha)
##relative_error_each = sp.array(relative_error_each)
##relative_error_alpha = sp.array(relative_error_alpha)
##mean_value_each = [0] * type_fiber
##mean_value_alpha = [0] * type_fiber
##for i_type in sp.arange(type_fiber):
##    each_kind_fiber = []
##    each_kind_fiber_alpha = []
##    for i_iteration in sp.arange(iteration):
##        index = i_type + 2 * i_iteration
##        each_kind_fiber.append(relative_error_each[index])
##        each_kind_fiber_alpha.append(relative_error_alpha[index])
##    each_kind_fiber = sp.array(each_kind_fiber)
##    each_kind_fiber_alpha = sp.array(each_kind_fiber_alpha)
##    mean_value_each[i_type] = np.mean(each_kind_fiber, axis = 0)
##    mean_value_alpha[i_type] = np.mean(each_kind_fiber_alpha, axis = 0)
####print 'relative error for each kind of fibers', mean_value_each, mean_value_alpha
####raw_input("enter")        
#####sum_error_each_pos = sp.zeros(len(relative_error_each[0]))
#####sum_error_alpha_pos = sp.zeros(len(relative_error_alpha[0]))
####for i_error in sp.arange(len(relative_error_each)):
####    for i_suberror in sp.arange(len(relative_error_each[i_error])):
####        sum_error_each_pos[i_suberror] += relative_error_each[i_error][i_suberror]
####        sum_error_alpha_pos[i_suberror] += relative_error_alpha[i_error][i_suberror]
####average_relative_error_each = sum_error_each_pos / iteration
####average_relative_error_alpha = np.mean(sum_error_alpha_pos) / iteration
####mean_value_each = np.mean(relative_error_each, axis = 0)
####mean_value_alpha = np.mean(relative_error_alpha, axis = 0)
####print 'relative error for proportionally moving', average_relative_error_each, \
####        mean_value_each
####print 'relative error for fixed value of moving', average_relative_error_alpha, \
####        mean_value_alpha
        
plot_ratio_function(zone_position, each_time_ratio, type_fiber, grid.prob_area)
raw_input('show the figure')
compare_relative_error(mean_value_each, mean_value_alpha, mean_value_alpha_1,
                        nrzones, type_fiber)

#meanraw_input("finish one loop for one kind of fiber")
    #for prob in probs:
        #prob has zone_point and ratio_each, here we plot prob func of the fiber
        #TODO PEI
print 'layout created in directory temp'
raw_input('Press key to quit example')

    
