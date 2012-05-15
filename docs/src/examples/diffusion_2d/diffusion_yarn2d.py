import os, sys
import numpy as np
import scipy as sp
import pylab
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon, Ellipse
from matplotlib.collections import PatchCollection

path = os.path.dirname(os.path.abspath(__file__))
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
number_fiber = 144
blend = [51.4, 48.6]
eps_value = 0.001
fiber_config = ['tmpfiber1.ini', 'tmpfiber2.ini']
prob_area = '[ lambda r: 0.0130* (-7.0349674 * r ** 4 + 13.2165906 * r ** 3. - 8.34242634 * r ** 2. + 2.1560694 * r + 0.00433641299), lambda r: 0.0185 * (1.13227791 * r **4 - 4.14156084 * r ** 3 + 4.39219104 * r ** 2 - 1.91494891 * r + 0.53275704)]'
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
radius_pure_fiber = 0.0503
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
radius_pure_fiber = 0.05441
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

from stick.fiber.config import FiberConfigManager
from stick.yarn.config import YarnConfigManager
from stick.lib.utils.utils import set_outputdir

cfgf1 = FiberConfigManager.get_instance('tmpfiber1.ini', realdatastr=ini_fiber1)
cfgf2 = FiberConfigManager.get_instance('tmpfiber2.ini', realdatastr=ini_fiber2)
cfg = YarnConfigManager.get_instance('tmpyarn.ini', realdatastr=ini_yarn)

if not os.path.isdir('temp'):
  os.mkdir('temp')
set_outputdir('temp')
#create a 2D grid for simulate DEET diffusion process
from stick.yarn2d.yarn2dmodel import Yarn2DModel
Yarn2DModel(cfg)

print 'calculate DEET diffusion'
raw_input('Press key to quit the calculation')
