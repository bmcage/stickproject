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
This package implements config defaults for Diffusion in yarn
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
                     'defaultfiber.ini'

LONGOPTS = ["inifile", 'outputdir']
SHORTOPTS = "i:o" 

METHOD = {
    'FVM': ('Finite Volume Method discretization', ['odeint', 'ode', 'fipy']),
    }

#possible fiber materials, map to diff coeff of components
YARN_MAT = {
    'YARN_1': ([0.015], ),
    }

#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------


class FiberfipyConfigManager(ConfigManager):

    __instance = {}
    
    def get_instance(inifile):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        if not (inifile in FiberfipyConfigManager.__instance):
            FiberfipyConfigManager.__instance[inifile] = None # Set for __init__()
            FiberfipyConfigManager.__instance[inifile] = FiberfipyConfigManager(inifile)
        return FiberfipyConfigManager.__instance[inifile]
    get_instance = staticmethod(get_instance)
    
    def __init__(self, filename = INIFILE_DEFAULT):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if filename not in FiberfipyConfigManager.__instance:
            raise Exception("This class is a singleton per filename. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename)

    def register_defaults(self):
        """default ini settings for a DiffusionIT problem"""
        self.register("general.read", False)
        self.register("general.verbose", False)
        self.register("general.method", 'FVM')
        self.register("general.submethod", 'ode')

        self.register("fiber.radius_pure_fiber", 0.01)
        self.register("fiber.radius_fiber", 0.0117)
        self.register("fiber.n_edge", 41)
        self.register("fiber.nrlayers", 1)

        self.register("fiberlayer_1.n_edge", 41)
        self.register("fiberlayer_1.thickness", 0.0017)
        self.register("fiberlayer_1.diffusion_coef", 5.2e-9)
        self.register("fiberlayer_1.init_conc", 'lambda x: 0.70')

        self.register("diffusion.diffusion_polymer_exp_factor", 0.)

        self.register("boundary.boundary_fib_left", 0.0)
        self.register("boundary.boundary_fib_right", 0.0)
        
        self.register("transfer.transfer_conc1", 7.2e-6)
        
        self.register("time.time_period", 500.)
        self.register("time.dt", 5.0)
