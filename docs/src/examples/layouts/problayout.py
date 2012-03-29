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


# We start with inifile settings in 
ini_yarn = u"""
[general]
method = 'FVM'
submethod = 'fipy'
read = False
verbose = False
[domain]
fiberlayout_method = 'virtlocoverlap'
distribute_fiber = ' integral'
cellsize_centre = 1.0e-1
cellsize_fiber = 2.50e-2
yarnradius = 1.0
[fiber]
number_type = 2
number_fiber = 400 
blend = [75., 25.]
eps_value = 0.001
fiber_config = ['tmpfiber1.ini', 'tmpfiber2.ini']
prob_area = '[lambda r: 0.00665 * ( 3.94231736 * r ** 4. - 8.00599612 * r ** 3. + 4.32590842 * r ** 2. - 0.67755321 * r + 0.41283374), lambda r: 0.02 * (4.44996934 * r ** 4. -9.05211446 * r ** 3. + 5.07757628 * r ** 2. -0.91234278 * r + 0.43632983)]'
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
radius_pure_fiber = 0.030
form = 'circle'
nrlayers = 2
mean_deviation = 0.00739
[fiberlayer_0]
thickness = 0.001
[fiberlayer_1]
thickness = 0.001
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
radius_pure_fiber = 0.05313
form = 'ellipse'
eccentricity = 0.7
nrlayers = 2
mean_deviation = 0.01201728333018
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

from stick.yarn2d.fiber_layout import plot_yarn
from stick.yarn2d.arearatioprobability import (calculate_proportion, 
            plot_ratio_function)

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
from stick.fiber.config import FiberConfigManager
from stick.yarn.config import YarnConfigManager
from stick.lib.utils.utils import set_outputdir
cfgf1 = FiberConfigManager.get_instance('tmpfiber1.ini', realdatastr=ini_fiber1)
cfgf2 = FiberConfigManager.get_instance('tmpfiber2.ini', realdatastr=ini_fiber2)
cfg = YarnConfigManager.get_instance('tmpyarn.ini', realdatastr=ini_yarn)
#create outputdir if not existing
if not os.path.isdir('temp'):
    os.mkdir('temp')
set_outputdir('temp')
#create 10 2D grids for statistics
from stick.yarn2d.yarn2dgrid import Yarn2dGrid
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
                }
from stick.yarn2d.fiber_layout import virtlocoverlaplayout
#After generating n times of iteration, plot prob func result from the average ratio value
coefficient_function = [0.00665, 0.02]
iteration = 5
each_time_ratio = []
relative_error_each = []
relative_error_alpha = []
for i in range(iteration):
    x_position, y_position, all_radius_fibers, \
                    fiber_kind,type_fiber, ratio_each_kind, \
                    zone_position_each, ratio_each_alpha = virtlocoverlaplayout(ouroptions)
    plot_yarn(x_position, y_position, all_radius_fibers, 
              fiber_kind, title='Realization %d' % i)
    ratio_each = [0] * type_fiber
    zone_position = [0] * type_fiber
    #zone_position = [0]
    for i_type in sp.arange(type_fiber):
        x_position_cal = []
        y_position_cal = []
        radius_each = []
        for i_fiber in sp.arange(len(fiber_kind)):
            if fiber_kind[i_fiber] == i_type:
                x_position_cal.append(x_position[i_fiber])
                y_position_cal.append(y_position[i_fiber])
                radius_each.append(all_radius_fibers[i_fiber])
        x_position_cal = np.array(x_position_cal)
        y_position_cal = np.array(y_position_cal)
        radius_each = np.array(radius_each)
##        probs = calculate_proportion(grid.radius_yarn, all_radius_fibers, 
##                x_position_cal, y_position_cal, nrzones=5)
        nrzones = 5
        zone_position[i_type] = sp.zeros(len(zone_position_each) + 1, float)
        zone_position[i_type][0] = 0.
        zone_position[i_type][1:] = zone_position_each[:]
        zone_position[i_type][-1] = 1.
#calculate the relative error in the calculation
        value_from_function = grid.prob_area[i_type]
        ratio_compare = value_from_function(zone_position[i_type]) / coefficient_function[i_type]
        print 'the ratio value from the calculation', ratio_compare
        raw_input('to generate the relative error')
##        if i_type == 0:
##            zone_position[i_type] = sp.zeros(len(probs[0]) +  2, float)
##            print 'check the length', len(zone_position[i_type][1:-1]), len(probs[0][:])
##            zone_position[i_type][1:-1] = probs[0][:]
##            zone_position[i_type][-2] = 1. - grid.radius_yarn / ((nrzones-1)*2 + 1)
##            zone_position[i_type][-1] = 1.
##            ratio_each[i_type] = sp.zeros(len(probs[0]) + 1, float)
##            ratio_each[i_type][1:] = probs[-1][:]
##        else:
##        zone_position[i_type] =sp.zeros(len(probs[0]) + 1, float)
##        zone_position[i_type][0] = 0.
##        zone_position[i_type][1:] = probs[0][:]
##        #zone_position[i_type][-2] = 1. - grid.radius_yarn / ((nrzones-1)*2 + 1)
##        zone_position[i_type][-1] = 1.

        #ratio_each[i_type] = sp.ones(len(probs[-1]), float)
        #ratio_each[i_type][:] = probs[-1][:]
        #the ratio value on the yarn's edge
        
        print 'the ratio value for eacn kind of fiber', ratio_each_kind[i_type][:]
        relative_error_each_kind = sp.power(abs(ratio_compare[:-1] - ratio_each_kind[i_type][:])\
                                            / ratio_compare[:-1], 2.)
        relative_error_kind_alpha = sp.power(abs(ratio_compare[:-1] - ratio_each_alpha[i_type][:])\
                                            /ratio_compare[:-1], 2.)
        raw_input('check the ratio value from the calculation')
        each_time_ratio.append(ratio_each_kind[i_type])
        relative_error_each.append(relative_error_each_kind)
        relative_error_alpha.append(relative_error_kind_alpha)
relative_error_each = sp.array(relative_error_each)
relative_error_alpha = sp.array(relative_error_alpha)
sum_error_each_pos = sp.zeros(len(relative_error_each[0]))
sum_error_alpha_pos = sp.zeros(len(relative_error_alpha[0]))
for i_error in sp.arange(len(relative_error_each)):
    for i_suberror in sp.arange(len(relative_error_each[i_error])):
        sum_error_each_pos[i_suberror] += relative_error_each[i_error][i_suberror]
        sum_error_alpha_pos[i_suberror] += relative_error_alpha[i_error][i_suberror]
average_relative_error_each = sum_error_each_pos / iteration
average_relative_error_alpha = np.mean(sum_error_alpha_pos) / iteration
mean_value_each = np.mean(relative_error_each, axis = 0)
mean_value_alpha = np.mean(relative_error_alpha, axis = 0)
print 'relative error for proportionally moving', average_relative_error_each, \
        mean_value_each
print 'relative error for fixed value of moving', average_relative_error_alpha, \
        mean_value_alpha
        
compare_relative_error(mean_value_each, mean_value_alpha, nrzones)
plot_ratio_function(zone_position, each_time_ratio, type_fiber, grid.prob_area)
#meanraw_input("finish one loop for one kind of fiber")
    #for prob in probs:
        #prob has zone_point and ratio_each, here we plot prob func of the fiber
        #TODO PEI
print 'layout created in directory temp'
raw_input('Press key to quit example')

def compare_relative_error(mean_value_each, mean_value_alpha, nrzones):
    ind = sp.arange(nrzones)
    width = 0.2
    
    plt.figure()
    propor_draw = plt.bar(ind, mean_value_each, width, color = 'b')
    alpha_draw = plt.bar(ind + width, mean_value_alpha, color = 'y')
    
    plt.ylabel(ur'$\Delta D = \frac{\sum\limits_{i=1}^{n}(p_{f}^{i} - p_{a}^{i})^2}{p_{f}^{i}}$')
    plt.xlabel(ind + width, ('Zone1', 'Zone2', 'Zone3', 'Zone4', 'Zone4'))
    plt.legend((propor_draw[0], alpha_draw[0]), ('$\alpha = \frac{R_{f}^{m}}{R_{f}^{n}}$', '$\alpha = \frac{1}{2}$'))
    
    plt.show()
    
