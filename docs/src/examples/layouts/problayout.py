"""
Example showing how a fiber-yarn layout can be constructed from a given
area probability function
"""

import os, sys
import numpy as np
import scipy as sp
import pylab
from matplotlib.patches import Circle, Wedge, Polygon, Ellipse
from matplotlib.collections import PatchCollection


# We start with inifile settings in 
ini_yarn = """
[general]
method = 'FVM'
submethod = 'fipy'
read = False
verbose = False
[domain]
fiberlayout_method = 'virtlocoverlap'
cellsize_centre = 1.0e-1
cellsize_fiber = 2.50e-2
yarnradius = 1.0
[fiber]
number_type = 2
number_fiber = 144
blend = [51.4, 48.6]
eps_value = 0.001
fiber_config = ['tmpfiber1.ini', 'tmpfiber2.ini']
prob_area = '[ lambda r: 0.0230 * (2.67999478 * r ** 4. -5.66339751 * r ** 3. + 3.07086875 * r ** 2. - 0.26131345 * r + 0.17273968), lambda r: 0.0275 * (1.13227791 * r **4 - 4.14156084 * r ** 3 + 4.39219104 * r ** 2 - 1.91494891 * r + 0.53275704)]'
[plot]
maxval = 0.0005
plotevery = 10
[writeout]
writeevery = 100
"""

ini_fiber1 = """
[general]
method = 'FVM'
submethod = 'odew'
read = False
verbose = True
[fiber]
radius_pure_fiber = 0.052
form = 'circle'
nrlayers = 2
mean_deviation = 0.0074
[fiberlayer_0]
thickness = 0.00085
[fiberlayer_1]
thickness = 0.00085
[plot]
plotevery = 1
"""
ini_fiber2 = """
[general]
method = 'FVM'
submethod = 'odew'
read = False
verbose = False
[fiber]
radius_pure_fiber = 0.055
form = 'ellipse'
eccentricity = 0.7
nrlayers = 2
mean_deviation = 0.012
[fiberlayer_0]
thickness = 0.001
[fiberlayer_1]
thickness = 0.001
[plot]
plotevery = 10
"""
#create the ini files
yarnf = open('tmpyarn.ini', 'w')
yarnf.write(ini_yarn)
yarnf.close()
fibf = open('tmpfiber1.ini', 'w')
fibf.write(ini_fiber1)
fibf.close()
fibf = open('tmpfiber2.ini', 'w')
fibf.write(ini_fiber2)
fibf.close()


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
radius_cotton = np.array(radius_cotton)
fiber_kind = np.zeros(len(radius_real_fiber), int)
fiber_kind[len(data_polyester):] = 1

from yarn2d.fiber_layout import plot_yarn
from yarn2d.arearatioprobability import (calculate_proportion, 
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
ratio_poly = sp.zeros(len(prob_poly[-1]) + 1, float)
ratio_cotton = sp.zeros(len(prob_cotton[-1])+1 , float)
ratio_poly[:-1] = prob_poly[-1][:]
ratio_cotton[:-1] = prob_cotton[-1][:]
print len(zone_position_real), len(ratio_poly), len(x_position_polyester), len(x_position_cotton)
poly_polyester = np.poly1d(np.polyfit(zone_position_real, ratio_poly, 4))
poly_cotton = np.poly1d(np.polyfit(zone_position_real, ratio_cotton, 4))
draw_real = sp.linspace(0., 1.0, 50)
print np.polyfit(zone_position_real, ratio_poly, 4)
print np.polyfit(zone_position_real, ratio_cotton, 4)
raw_input("record the coefficients of polynomial equations")

pylab.figure()
pylab.subplot(211)
pylab.plot(zone_position_real, ratio_poly, 'o')
pylab.plot(draw_real, poly_polyester(draw_real), '--')
pylab.xlim(0., 1.05)
pylab.ylim(0., 0.8)
pylab.subplot(212)
pylab.plot(zone_position_real, ratio_cotton, '*')
pylab.plot(draw_real, poly_cotton(draw_real), '-')
pylab.xlim(0., 1.05)
pylab.ylim(0., 0.8)
pylab.axis()
pylab.show()

#set up a yarn computation
from yarn.config import YarnConfigManager
from lib.utils.utils import set_outputdir
cfg = YarnConfigManager.get_instance('tmpyarn.ini')
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
                }
from yarn2d.fiber_layout import virtlocoverlaplayout
#After generating n times of iteration, plot prob func result from the average ratio value

iteration = 10
each_time_ratio = [] 
for i in range(iteration):
    x_position, y_position, all_radius_fibers, \
                    fiber_kind, type_fiber = virtlocoverlaplayout(ouroptions)
    plot_yarn(x_position, y_position, all_radius_fibers, 
              fiber_kind, title='Realization %d' % i)
    ratio_each = [0] * type_fiber
    #zone_position = [0]
    for i_type in sp.arange(type_fiber):
        x_position_cal = []
        y_position_cal = []
        for i_fiber in sp.arange(len(fiber_kind)):
            if fiber_kind[i_fiber] == i_type:
                x_position_cal.append(x_position[i_fiber])
                y_position_cal.append(y_position[i_fiber])
        x_position_cal = np.array(x_position_cal)
        y_position_cal = np.array(y_position_cal)
        probs = calculate_proportion(grid.radius_yarn, all_radius_fibers, 
                x_position_cal, y_position_cal, nrzones=5)
        zone_position =sp.zeros(len(probs[0]) + 1, float)
        zone_position[0] = 0.
        zone_position[1:] = probs[0][:]
        zone_position[-1] = 1.

        ratio_each[i_type] = sp.ones(len(probs[-1]), float)
        ratio_each[i_type][:] = probs[-1][:]
        #the ratio value on the yarn's edge
        each_time_ratio.append(ratio_each[i_type])

plot_ratio_function(zone_position, each_time_ratio, type_fiber, grid.prob_area)
#meanraw_input("finish one loop for one kind of fiber")
    #for prob in probs:
        #prob has zone_point and ratio_each, here we plot prob func of the fiber
        #TODO PEI
print 'layout created in directory temp'
raw_input('Press key to quit example')
