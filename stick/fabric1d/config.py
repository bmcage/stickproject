#
# Copyright (C) 2005-2007  Donald N. Allingham
# Copyright (C) 2008-2009  Gary Burton 
# Copyright (C) 2009       Doug Blank <doug.blank@gmail.com>
# Copyright (C) 2009       Benny Malengier <bm@cage.ugent.be>
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
#


"""
This package implements config defaults for DiffusionIT
"""
#---------------------------------------------------------------
#
# local imports
#
#---------------------------------------------------------------
import os
import const
from lib.config import ConfigManager

#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'fiber' + os.sep + 'default.ini'

LONGOPTS = ["inifile", 'outputdir']
SHORTOPTS = "i:o" 

METHOD = {
    'FVM': ('Finite Volume Method discretization', ['cvode']),
    }

#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------

class Fabric2DConfigManager(ConfigManager):

    __instance = None
    
    def get_instance(inifile):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        if FiberConfigManager.__instance is None:
            FiberConfigManager.__instance = 1 # Set to 1 for __init__()
            FiberConfigManager.__instance = FiberConfigManager(inifile)
        return FiberConfigManager.__instance
    get_instance = staticmethod(get_instance)
    
    def __init__(self, filename = INIFILE_DEFAULT):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if FiberConfigManager.__instance is not 1:
            raise Exception("This class is a singleton. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename)

    def register_defaults(self):
        """default ini settings for a DiffusionIT problem"""
        self.register("general.read", False)
        self.register("general.verbose", False)
        self.register("general.method", 'FVM')
        self.register("general.submethod", 'fipy')
        
        self.register("domain.domain_size", [5.0e-2, 1.5e-1])
        self.register("domain.dx", 2.5e-3)
        self.register("domain.dy", 2.5e-3)
        
        self.register("sample.size_sample", [2e-2, 4.92e-3])
        self.register("sample.yarn_config", ['../yarn2d/defaultyarn.ini'])
        
        self.register("diffusion_DEET", 5.0e-8)
        self.register("domain.tortuosity_fab", 2.)
        self.register("domain.diff_DEET_void", 7.0e-8)
        
        self.register("initial.initial_conc_DEET", 0.0)
        self.register("initial.initial_DEET_void", 0.0)
        
        self.register("boundary.boundary_up", 0.0)
        self.register("boundary.boundary_bottom", 0.0)
        self.register("boundary_left", 0.0)
        self.register("boundary_right", 0.0)
        
        self.register("time.time_period", 5000)
        self.register("time.dt", 100.0)
        
        
        """
        self.register("general.method", 'FVM')
        self.register("general.submethod", METHOD['FVM'][1][0])
        self.register("general.components", ['DEET'])
        self.register("general.inverseproblem", False)
        
        self.register("init.initfromfile", False)
        self.register("init.initfile", "")
        self.register("init.initfunc", "lambda x: (0.01, 0.01)")
        
        self.register('time.type', 'fixstep')
        self.register('time.t0', 0.)
        self.register('time.dt', 5.0)
        self.register('time.nrsteps', 800)
        self.register('time.tend', 0.0)
        self.register('time.tout_every', 1.)
        
        self.register('fiber.material', 'polyester')
        self.register('fiber.diffcoef', 
                        [FIBER_MAT['polyester'][0][COMPONENTS['DEET'][0]]])
        self.register('fiber.radius', 1.)
        self.register('fiber.discrpoints', 101)
        self.register('fiber.nrlayer', 1)
        
        self.register('layer_0.thickness', 0.2)
        bind = BINDER.keys()[0]
        self.register('layer_0.binder', bind)
        self.register('layer_0.diffcoef', 
                        [BINDER[bind][0][COMPONENTS['DEET'][0]]])
        self.register('layer_0.discrpoints', 21)

        self.register("bcval.bcleft_type", ['homNeu', 'homNeu'])
        self.register("bcval.bcright_type", ['homNeu', 'homNeu'])
        self.register("bcval.bcleft_val", [0.0, 0.0])
        self.register("bcval.bcright_val", [0.0, 0.0])

        self.register("diffusion.type", 'piecewise')
        self.register("diffusion.defaultpar", True)
        self.register("diffusion.par", [[0.1,0.6,0.1]])

        self.register("experiments.nrexp", 0)
        """
