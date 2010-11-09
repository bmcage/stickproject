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
INIFILE_DEFAULTFAB = const.INI_DIR + os.sep + 'yarn2d' + os.sep + \
                     'defaultyarn.ini'

LONGOPTS = ["inifile", 'outputdir']
SHORTOPTS = "i:o" 

METHOD = {
    'FVM': ('Finite Volume Method discretization', ['cvode']),
    }

#possible fiber materials, map to diff coeff of components
'''

''' 
YARN_MAT = {
    'YARN_1': ([0.015], ),
    }

#all components possible. Map to a unique number used in value maps, eg
#  DEET --> 0 ==> diffusion is BINDER[key][COMPONENT['DEET'][0]
COMPONENTS = {
    'DEET': [0], 
    }

#binders and typical diff coef for Deet and permethrin


#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------


class Yarn2dConfigManager(ConfigManager):

    __instance = None
    
    def get_instance(inifile):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        if Yarn2dConfigManager.__instance is None:
            Yarn2dConfigManager.__instance = 1 # Set to 1 for __init__()
            Yarn2dConfigManager.__instance = Yarn2dConfigManager(inifile)
        return Yarn2dConfigManager.__instance
    get_instance = staticmethod(get_instance)
    
    def __init__(self, filename = INIFILE_DEFAULTFAB):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if Yarn2dConfigManager.__instance is not 1:
            raise Exception("This class is a singleton. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename)

    def register_defaults(self):
        """default ini settings for a DiffusionIT problem"""
        self.register("general.read", False)
        self.register("domain.cellsize_centre", 5.0e-2)
        self.register("domain.cellsize_fiber", 5.0e-2)
        self.register("domain.radius", 1.)
        
        self.register("fiber.type", 'constant')
        self.register("fiber.number_fiber", 20)
        self.register("fiber.radius_fiber",0.1)
        self.register("fiber.n_point", 51)
        
        self.register("initial.init_conc1", 0.)
        self.register("initial.init_conc1_fiber", 'lambda x:(0.2,0.0)')
        
        self.register("diffusion.diffusion_conc1", 2e-5)
        self.register("diffusion.diffusion_co_l1", 5.0e-6)
        self.register("diffusion.diffusion_co_l2", 5.0e-6)
        
        self.register("boundary.boundary_fib_left", 0.0)
        self.register("boundary.boundary_fib_right", 1.0)
        
        self.register("transfer.transfer_conc1", 7.2e-11)
        
        self.register("time.time_period", 20.)
        self.register("time.dt", 0.1)
        
        