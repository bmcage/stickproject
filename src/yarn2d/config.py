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
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'yarn2d' + os.sep + \
                     'defaultyarn.ini'

LONGOPTS = ["inifile", 'outputdir']
SHORTOPTS = "i:o" 

METHOD = {
    'FVM': ('Finite Volume Method discretization', ['fipy']),
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
    
    def __init__(self, filename = INIFILE_DEFAULT):
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
        self.register("general.verbose", False)
        self.register("general.method", 'FVM')
        self.register("general.submethod", 'fipy')
        
        self.register("domain.cellsize_centre", 5.0e-2)
        self.register("domain.cellsize_fiber", 5.0e-2)
        self.register("domain.yarnradius", 1.)

        self.register("fiber.number_type", 2)
        self.register("fiber.number_fiber", 60)
        self.register("fiber.eps_value", 0.001)
        self.register("fiber.blend", [100.0])
        self.register("fiber.radius_fiber", [0.0117])
        self.register("fiber.fiber_config", ['../fiber/defaultfiber.ini'])
        
        self.register("initial.init_conc", 0.)
        
        self.register("diffusion.diffusion_conc", 2e-5)
        
        self.register("boundary.boundary_exterior", 0.0)
        self.register("boundary.boundary_interior", 1.0)

        self.register("size_hole.net_width", 1.0e-3)
        self.register("size_hole.net_length", 2.0e-3)
        self.register("size_hole.leng_yarn", 2.0)
        self.register("size_hole.domain_effect", 0.02)
        self.register("size_hole.dis_effect", 4)
        
        self.register("time.time_period", 4000.)
        self.register("time.dt", 5.0)
        
        