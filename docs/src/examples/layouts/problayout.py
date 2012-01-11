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

#read the statistic data from a real yarn2d
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
fradius_real_fiber = np.array(radius_real_fiber)
radius_poly = radius_real_fiber[:(len(data_polyester) + 1)]
radius_cotton = radius_real_fiber[(len(data_polyester) + 1):]
x_polyester = x_position_real_fiber[:(len(data_polyester) + 1)]
x_cotton = x_position_real_fiber[(len(data_polyester) + 1):]
y_polyester = y_position_real_fiber[:(len(data_polyester) + 1)]
y_cotton = y_position_real_fiber[(len(data_polyester) + 1):]
area_cotton_fiber = sp.pi * sp.power(radius_cotton, 2.)
area_poly_fiber = sp.pi * sp.power(radius_poly, 2.)
fig = pylab.figure()
ax = fig.add_subplot(111, xlim = (-1.1, 1.1), ylim = (-1.1, 1.1))
patches_1 = []
patches_2 = []
    
for x_center, y_center, radii in zip(x_polyester[:], y_polyester[:], radius_poly[:]):
  circle = Circle((x_center, y_center), radii, facecolor = 'g', alpha = 0.4)
  patches_1.append(circle)
for x_center, y_center, radii in zip(x_cotton[:], y_cotton[:], radius_cotton[:]):
  circle = Circle((x_center, y_center), radii, facecolor = 'r', alpha = 0.4)
  patches_2.append(circle)
circle = Circle((0., 0.), 1.0)
patches_1.append(circle)
p_1 = PatchCollection(patches_1, facecolor = 'red', cmap = matplotlib.cm.jet, alpha = 0.4)
p_2 = PatchCollection(patches_2, facecolor = 'black', cmap = matplotlib.cm.jet, alpha = 0.4)
ax.add_collection(p_1)
ax.add_collection(p_2)
pylab.ion()
pylab.draw()
#pylab.ion()

#set up a yarn computation
from yarn2d.config import Yarn2dConfigManager
from lib.utils.utils import set_outputdir
cfg = Yarn2dConfigManager.get_instance('tmpyarn.ini')
#create outputdir if not existing
if not os.path.isdir('temp'):
    os.mkdir('temp')
set_outputdir('temp')
#create 10 2D grids for statistics
from yarn2d.yarn2dgrid import Yarn2dGrid
grid = Yarn2dGrid(cfg)
for i in range(1):
    mesh2d = grid.mesh_2d_generate(filename='yarn%0d.geo'%i,
                          regenerate=True)
print 'layout created in directory temp'
raw_input('Press key to quit example')

