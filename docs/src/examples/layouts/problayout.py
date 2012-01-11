"""
Example showing how a fiber-yarn layout can be constructed from a given
area probability function
"""

import os, sys

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
for i in range(10):
    mesh2d = grid.mesh_2d_generate(filename='temp'+os.sep+'yarn%0d.geo'%i,
                          regenerate=True)
print 'layout created in directory temp'
