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
This package implements config defaults for RoomFabric
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
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'room' + os.sep + 'defaultroom.ini'

LONGOPTS = ["inifile", 'outputdir', 'write-ini']
SHORTOPTS = "i:o" 

NONE = 0
BOTCENT = 1
PLACEMENT = {
    'none': NONE,
    'bottomcenter': BOTCENT
    }
#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------

class RoomConfigManager(ConfigManager):

    __instance = {}
    
    def get_instance(inifile, realdatastr=None):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        inifilebase = os.path.basename(inifile)
        if inifile in RoomConfigManager.__instance:
            return RoomConfigManager.__instance[inifile]
        elif inifilebase in RoomConfigManager.__instance:
            return RoomConfigManager.__instance[inifilebase]
        else:
            RoomConfigManager.__instance[inifile] = None # Set for __init__()
            RoomConfigManager.__instance[inifile] = RoomConfigManager(inifile,
                                                                realdatastr)
            RoomConfigManager.__instance[inifilebase] = RoomConfigManager.__instance[inifile]
        return RoomConfigManager.__instance[inifile]
    get_instance = staticmethod(get_instance)

    def delete(inifile):
        """remove the instance inifile from the loaded configurations"""
        del RoomConfigManager.__instance[inifile]
        if inifile != os.path.basename(inifile):
            del RoomConfigManager.__instance[os.path.basename(inifile)]
    delete = staticmethod(delete)

    def __init__(self, filename = INIFILE_DEFAULT, realdatastr=None):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if (filename not in RoomConfigManager.__instance) or (
                RoomConfigManager.__instance[filename] is not None):
            raise Exception("This class is a singleton per filename. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename, realdatastr)

    def register_defaults(self):
        """default ini settings for a DiffusionIT problem"""
        self.register("general.read", False)
        self.register("general.verbose", False)
        
        self.register("domain.type", 'beam',
            "type of domain. options are: beam")
        
        self.register("domain.length", 200.,
            "If domain type has length, it's length in mm")
        self.register("domain.width", 200.,
            "If domain type has width, it's width in mm")
        self.register("domain.height", 100.,
            "If domain type has height, it's height in mm")
        self.register("domain.load_msh", False,
            "If True, msh_file is loaded instead of constructing the domain again."
            " The msh_file should be consistent with domain and fabric values !!!")
        self.register("domain.msh_file", 'room.geo',
            "The msh file to load")
        self.register("fabric.fabric_config", '../fiberfabric/defaultfiberfabric',
            "The ini file describing the fabric used in this room")
        self.register("fabric.fabricposition", 'bottomcenter',
            'Placement of the fabric in the room. Choose from ' + 
                ",".join(PLACEMENT.keys()))
        self.register("fabric.type", 'square',
            'The form of the fabric. Choose from: square. Should be consistent '
            'with the data in the fabric config file !!' )
        self.register("fabric.length", 100.,
            "If fabric form has length, it's length in mm. Should be consistent "
            'with the data in the fabric config file !!' )
        self.register("fabric.width", 100.,
            "If fabric form has width, it's width in mm. Should be consistent "
            'with the data in the fabric config file !!' )
        self.register("fabric.height", 3.,
            "If fabric form has height, it's height in mm. Should be consistent "
            'with the data in the fabric config file !!' )
        
        self.register("initial.init_concvap", 'lambda x,y,z: 1.0', 
            'Initial vapour concentration in kg/m^3.')
        self.register("initial.init_concair", 'lambda x,y,z: 1.0', 
            'Initial vapour concentration in kg/m^3.')
        self.register("initial.init_temp", 'lambda x,y,z: 20.0', 
            'Initial temperature in degrees celcius')
    
        self.register("roomcoeff.diff_coef", 25.,
            "Diffusion coefficient in the room in mm^2/s")
        self.register("roomcoeff.therm_cond_K", 10.,
            "Thermal conductivity of air in W / (m K)")
        self.register("roomcoeff.spec_heat_c", 10.,
            "Volumetric Specific heat of air J/m^3")

        self.register("boundary.dirichletval_T_BC", 21., 
            "Dirichilet boundary condition for Temperature in degree Celcius")

        self.register("discretization.el_length", 10,
            "Number of elements along the length")
        self.register("discretization.el_width", 10,
            "Number of elements along the width")
        self.register("discretization.el_height", 6,
            "Number of elements along the height")

        #info about the diffusion outside fabric
        self.register("diffusion.diff_coef", 5.0e-8,
            "Diffusion in air of the active component, in mm^2/s")

        #plotting output
        self.register("plot.plotevery", 10,
            "When plotting over time, indicate how many steps dt to skip "
            "before plotting again")
        self.register("plot.writeevery", 100,
            "When writing data out over time, indicate how many steps dt to skip")
        
        #time info: how long to run and timestep
        self.register("time.time_period", 5000.)
        self.register("time.dt", 0.1)
