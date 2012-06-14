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
from __future__ import division
import os
import const
from lib.config import ConfigManager

#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'fabric' + os.sep + 'defaultfabric.ini'

LONGOPTS = ["inifile", 'outputdir']
SHORTOPTS = "i:o" 

#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------

class BednetConfigManager(ConfigManager):

    __instance = {}
    
    def get_instance(inifile, realdatastr=None):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        inifilebase = os.path.basename(inifile)
        if inifile in BednetConfigManager.__instance:
            return BednetConfigManager.__instance[inifile]
        elif inifilebase in BednetConfigManager.__instance:
            return BednetConfigManager.__instance[inifilebase]
        else:
            BednetConfigManager.__instance[inifile] = None # Set for __init__()
            BednetConfigManager.__instance[inifile] = BednetConfigManager(inifile,
                                                                realdatastr)
            BednetConfigManager.__instance[inifilebase] = BednetConfigManager.__instance[inifile]
        return BednetConfigManager.__instance[inifile]
    get_instance = staticmethod(get_instance)

    def delete(inifile):
        """remove the instance inifile from the loaded configurations"""
        del BednetConfigManager.__instance[inifile]
        if inifile != os.path.basename(inifile):
            del BednetConfigManager.__instance[os.path.basename(inifile)]
    delete = staticmethod(delete)

    def __init__(self, filename = INIFILE_DEFAULT, realdatastr=None):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if (filename not in BednetConfigManager.__instance) or (
                BednetConfigManager.__instance[filename] is not None):
            raise Exception("This class is a singleton per filename. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename, realdatastr)

    def register_defaults(self):
        """default ini settings for a DiffusionIT problem"""
        self.register("general.read", False)
        self.register("general.verbose", False)
        
        self.register("observer.x0", [1e-3,])
        
        self.register("domain.nr_vert_yarns",50)
        self.register("domain.nr_hor_yarns",100)
        self.register("domain.domain_size", [5.0e-2, 1.5e-1], 'size fabric in m')
        self.register("domain.dx", 2.5, 'size hole in mm')
        self.register("domain.dy", 2.5, 'size hole in mm')
        
        self.register("sample.size_sample", [2e-2, 4.92e-3])
        self.register("sample.yarn_config", ['../yarn/defaultyarn.ini'])

        #size_hole section describing the square holes in the net
        self.register("size_hole.net_width", 1.0e-3,
            "the width of the hole of the bed net (mm)")
        self.register("size_hole.net_length", 2.0e-3,
            "the length of the hole of the bed net (mm)")
        self.register("size_hole.length_yarn", 2.0,
            "the length of the yarn of a bed net")
        self.register("size_hole.domain_effect", 0.02,
            "the domain of repelling the mosquito")
        self.register("size_hole.dis_effect", 4)

        self.register("diffusion.diff_coef", 5.0e-8)
        self.register("diffusion.tortuosity_fab", 2.)
        self.register("diffusion.diff_DEET_void", 7.0e-8)
        
        self.register("saturation.saturation_conc", 5.0)
        
        self.register("initial.init_conc", 0.0)
        self.register("initial.init_void", 0.0)
        
        self.register("boundary.boundary_up", 0.0)
        self.register("boundary.boundary_bottom", 0.0)
        self.register("boundary.boundary_left", 0.0)
        self.register("boundary.boundary_right", 0.0)
        
        self.register("plot.plotevery", 10,
            "When plotting over time, indicate how many steps dt to skip before plotting again")
        
        self.register("time.time_period", 5000.)
        self.register("time.dt", 100.0)
