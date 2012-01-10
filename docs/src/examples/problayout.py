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
[initial]
init_conc = 0.
[yarn]
tortuosity = 1.
[diffusion]
diffusion_conc = 64.8
[boundary]
boundary_exterior = 0.
boundary_interior = 1.
transfer_conc1 = 5.3e-9
evap_equilibrium = 0.362e-9
[time]
time_period = 10.
dt = 0.01
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
internal_diffusion = False
n_edge = 21
porosity_in = 0.2
mean_deviation = 0.0053
[fiberlayer_0]
n_edge = 21
thickness = 0.00085
diffusion_coef = 5.2e-7
porosity_layer = 1.
diffusion_polymer_exp_factor = 0.3
init_conc = 'lambda x: 0.961368'
porosity_layer = 1.0
[fiberlayer_1]
n_edge = 6
thickness = 0.00085
diffusion_coef = 7.2e-7
porosity_layer = 1.
diffusion_polymer_exp_factor = 0.3
init_conc = 'lambda x: 0.0'
[boundary]
type_left = 'flux'
type_right = 'evaporation'
boundary_fib_left = 0.0
boundary_fib_right = 0.0
transfer_right =  5.0e-10
evap_satconc = 'lambda T: 1.'
evap_transfer = 5.0e-10
evap_minbound = 0.
evap_Cf = 'lambda t: 0.'
[time]
time_period = 50000.
dt = 10.
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
internal_diffusion = True
n_edge = 21
diffusion_coef = 3.0e-8
porosity_in = 1.0
mean_deviation = 0.00541
[fiberlayer_0]
n_edge = 21
thickness = 0.001
diffusion_coef = 5.2e-7
porosity_layer = 1.
diffusion_polymer_exp_factor = 0.
init_conc = 'lambda x: 0.961368'
[fiberlayer_1]
n_edge = 21
thickness = 0.001
diffusion_coef = 7.2e-7
diffusion_polymer_exp_factor = 0.
init_conc = 'lambda x: 0.0'
porosity_layer = 1.
[boundary]
type_left = 'flux'
type_right = 'transfer'
boundary_fib_left = 0.0
boundary_fib_right = 0.0
transfer_right = 5.0e-14
[time]
time_period = 600.0
dt = 2.
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
#create a 2D grid
from yarn2d.yarn2dgrid import Yarn2dGrid
grid = Yarn2dGrid(cfg)
mesh2d = grid.mesh_2d_generate(filename='temp'+os.sep+'yarn01.geo',
                          regenerate=True)
print 'layout created in directory temp'
