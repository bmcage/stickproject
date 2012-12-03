#
# Copyright (C) 2005-2007  Donald N. Allingham
# Copyright (C) 2008-2009  Gary Burton 
# Copyright (C) 2009       Doug Blank <doug.blank@gmail.com>
# Copyright (C) 2009-2012  Benny Malengier <bm@cage.ugent.be>
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
This package implements config defaults for FiberFabric, which models a 
fabric consisting of fibers/pcm (so neglecting yarn scale)
"""
#---------------------------------------------------------------
#
# local imports
#
#---------------------------------------------------------------
from __future__ import division
import os

import stick.const as const
from stick.lib.config import ConfigManager

#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'fiberfabric' + os.sep + 'defaultfiberfabric.ini'

LONGOPTS = ["inifile", 'outputdir', 'write-ini']
SHORTOPTS = "i:o" 

#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------

class FiberFabricConfigManager(ConfigManager):

    __instance = {}
    
    def get_instance(inifile, realdatastr=None):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        inifilebase = os.path.basename(inifile)
        if inifile in FiberFabricConfigManager.__instance:
            return FiberFabricConfigManager.__instance[inifile]
        elif inifilebase in FiberFabricConfigManager.__instance:
            return FiberFabricConfigManager.__instance[inifilebase]
        else:
            FiberFabricConfigManager.__instance[inifile] = None # Set for __init__()
            FiberFabricConfigManager.__instance[inifile] = FiberFabricConfigManager(inifile,
                                                                realdatastr)
            FiberFabricConfigManager.__instance[inifilebase] = FiberFabricConfigManager.__instance[inifile]
        return FiberFabricConfigManager.__instance[inifile]
    get_instance = staticmethod(get_instance)

    def delete(inifile):
        """remove the instance inifile from the loaded configurations"""
        del FiberFabricConfigManager.__instance[inifile]
        if inifile != os.path.basename(inifile):
            del FiberFabricConfigManager.__instance[os.path.basename(inifile)]
    delete = staticmethod(delete)

    def __init__(self, filename = INIFILE_DEFAULT, realdatastr=None):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if (filename not in FiberFabricConfigManager.__instance) or (
                FiberFabricConfigManager.__instance[filename] is not None):
            raise Exception("This class is a singleton per filename. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename, realdatastr)

    def register_defaults(self):
        """default ini settings for a DiffusionIT problem"""
        self.register("general.read", False)
        self.register("general.verbose", False)
        
        self.register("fabric.type", 'square',
            "type of domain. options are: square")
        self.register("fabric.length", 100.,
            "If fabric type has length, it's length in mm")
        self.register("fabric.width", 100.,
            "If fabric type has width, it's width in mm")
        self.register("fabric.height", 3.,
            "If fabric type has height, it's height in mm")

        self.register("fiber.fiber_config", ['../fiber/defaultfiber.ini'],
            "List of ini files describing the fibers used in this fabric")
        self.register("fiber.volfrac", [0.3],
            'Volume fraction of the fibers. % of volume of fabric occupied by fiber')

        self.register("pcm.pcm_config", ['../pcm/defaultpcm.ini'],
            "List of ini files describing the pcms used in this fabric")
        self.register("pcm.volfrac", [0.03],
            'Volume fraction of the pcms. % of volume of fabric occupied by pcm')

        self.register("component.present", False,
            "Indicate if an active component is present. If not, it will not be"
            " tracked")

        self.register("initial.init_conc", 'lambda x,y,z: 0.0', 
            'Initial concentration of an active component we track in kg/m^3.')
        self.register("initial.init_concvap", 'lambda x,y,z: 1.0', 
            'Initial vapour concentration in kg/m^3.')
        self.register("initial.init_concair", 'lambda x,y,z: 1.0', 
            'Initial air concentration in kg/m^3.')
        self.register("initial.init_temp", 'lambda x,y,z: 20.0', 
            'Initial temperature in degrees celcius')
    
        self.register("fabriccoeff.diff_coef", 25.,
            "Diffusion coefficient in the room in mm^2/s")
        self.register("fabriccoeff.therm_cond_K", 0.025,
            "Thermal conductivity of air in W / (m K)")
        self.register("fabriccoeff.spec_heat_c", 1.21,
            "Volumetric Specific heat of air J/(mm^3 K). References: Air, 1.21;"
            " cotton, 1925.5; ABS, 1647.; Polyester, 1411.5")

        self.register("boundary.T_type", 'heatingplate', 
            "Type of boundary condition for temperature. Possibilies: "
            "  1. heatingplate = bottom on a heating plate"
            "  2. insulated    = bottom on insulated plate")
        self.register("boundary.T_dir", 55.,
            "Dirichilet boundary condition for Temperature in degree Celcius"
            " as appropriate for the set T_type")

        self.register("discretization.el_length", 8,
            "Number of elements along the length")
        self.register("discretization.el_width", 8,
            "Number of elements along the width")
        self.register("discretization.el_height", 4,
            "Number of elements along the height")

        #plotting output
        self.register("plot.plotevery", 10,
            "When plotting over time, indicate how many steps dt to skip "
            "before plotting again")
        self.register("plot.writeevery", 100,
            "When writing data out over time, indicate how many steps dt to skip")
        
        #time info: how long to run and timestep
        self.register("time.time_period", 5000.)
        self.register("time.dt", 0.1)
