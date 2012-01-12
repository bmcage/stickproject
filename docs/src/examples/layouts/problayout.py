"""
Example showing how a fiber-yarn layout can be constructed from a given
area probability function
"""

import os, sys
import numpy as np
import scipy as sp
import pylab
import matplotlib
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
prob_area = '[lambda r: 0.0173 * (0.360 * r**4  -3.397 * r ** 3 + 4.531 * r ** 2 - 1.979 * r + 0.496), lambda r: 0.01975 * (-1.061 * r ** 4 + 0.397 * r ** 3 + 0.606 * r ** 2 - 0.093* r + 0.152)]'
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
mean_deviation = 0.0053
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
mean_deviation = 0.00541
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
data_polyester = np.loadtxt('fiber_polyester.csv')
data_cotton = np.loadtxt('fiber_cotton.csv')
x_position_real_fiber = []
y_position_real_fiber = []
radius_real_fiber = []
for i_polyester in sp.arange(len(data_polyester)):
  x_position_real_fiber.append(data_polyester[i_polyester][0])
  y_position_real_fiber.append(data_polyester[i_polyester][1])
  radius_real_fiber.append(data_polyester[i_polyester][2])
for i_cotton in sp.arange(len(data_cotton)):
  x_position_real_fiber.append(data_cotton[i_cotton][0])
  y_position_real_fiber.append(data_cotton[i_cotton][1])
  radius_real_fiber.append(data_cotton[i_cotton][2])
x_position_real_fiber = np.array(x_position_real_fiber)
y_position_real_fiber = np.array(y_position_real_fiber)
radius_real_fiber = np.array(radius_real_fiber)
fiber_kind = np.zeros(len(radius_real_fiber), int)
fiber_kind[len(data_polyester):] = 1

from yarn2d.fiber_layout import plot_yarn
plot_yarn(x_position_real_fiber, y_position_real_fiber, radius_real_fiber, 
          fiber_kind, title='Real fiber-yarn layout')

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
for i in range(2):
    x_position, y_position, all_radius_fibers, \
                    fiber_kind = virtlocoverlaplayout(ouroptions)
    plot_yarn(x_position, y_position, all_radius_fibers, 
              fiber_kind, title='Realization %d' % i)
    probs = calculate_proportion(grid.radius_yarn, all_radius_fibers, x_position, 
                         y_position, fiber_kind, nrzones=5)
    for prob in probs:
        #prob has zone_point and ratio_each, here we plot prob func of the fiber
        TODO PEI

print 'layout created in directory temp'
raw_input('Press key to quit example')

