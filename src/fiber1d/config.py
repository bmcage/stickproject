#
# Copyright (C) 2005-2007  Donald N. Allingham
# Copyright (C) 2008-2009  Gary Burton 
# Copyright (C) 2009       Doug Blank <doug.blank@gmail.com>
# Copyright (C) 2009       Benny Malengier <bm@cage.ugent.be>
# Copyright (C) 2010       Pei.Li <pli@cage.ugent.be>
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
This package implements config defaults for diffusion in textile fibers
"""
#---------------------------------------------------------------
#
# local imports
#
#---------------------------------------------------------------
from __future__ import division
import os
import const
from lib.config import ConfigManager

#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'fiber' + os.sep + \
                     'defaultfiber2.ini'

LONGOPTS = ["inifile", 'outputdir', 'write-ini']
SHORTOPTS = "i:o" 

#solution methods and the possible submethods
METHOD = {
    'FVM': ('Finite Volume Method discretization', 
            ['odew', 'odeintw', 'odeu', 'odeintu', 'fipy']),
    }

FLUX = 0
TRANSFER = 1
BOUND_TYPE = {
    'flux': FLUX,
    'transfer': TRANSFER,
    }

#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------

class Fiber1dConfigManager(ConfigManager):

    __instance = {}
    
    def get_instance(inifile):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        if not (inifile in Fiber1dConfigManager.__instance):
            Fiber1dConfigManager.__instance[inifile] = None # Set for __init__()
            Fiber1dConfigManager.__instance[inifile] = Fiber1dConfigManager(inifile)
        return Fiber1dConfigManager.__instance[inifile]
    get_instance = staticmethod(get_instance)
    
    def __init__(self, filename = INIFILE_DEFAULT):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if filename not in Fiber1dConfigManager.__instance:
            raise Exception("This class is a singleton per filename. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename)

    def register_defaults(self):
        """default ini settings for a fiber1d problem"""
        self.register("general.read", False)
        self.register("general.verbose", False)
        self.register("general.method", 'FVM')
        self.register("general.submethod", 'odew')
        self.register("general.fiber_kind", 'polyester')

        self.register("fiber.radius_pure_fiber", 0.01)
        self.register("fiber.nrlayers", 1)
        self.register("fiber.internal_diffusion", False)
        self.register("fiber.diffusion_coef", 0.)
        self.register("fiber.diffusion_polymer_exp_factor", 0.)
        self.register("fiber.n_edge", 41)
        self.register("fiber.init_conc", 'lambda x: 0.')
        self.register("fiber.porosity_in", 0.2)
        self.register("fiber.percentage_active", 1.0)
        self.register("fiber.beta_value", 0.04)
        self.register("fiber.mean_deviation", 0.00355)

        self.register("fiberlayer_0.n_edge", 41)
        self.register("fiberlayer_0.thickness", 0.0017)
        self.register("fiberlayer_0.diffusion_coef", 5.2e-9)
        self.register("fiberlayer_0.init_conc", 'lambda x: 0.70')
        self.register("fiberlayer_0.porosity_layer", 1.0)

        self.register("boundary.type_left", 'flux')
        self.register("boundary.type_right", 'flux')
        self.register("boundary.boundary_fib_left", 0.0)
        self.register("boundary.boundary_fib_right", 0.0)
        self.register("boundary.transfer_right", 0.0)
        
        self.register("time.time_period", 500.)
        self.register("time.dt", 5.0)
