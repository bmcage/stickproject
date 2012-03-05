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
This package implements config defaults for phase change material capsules
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
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'pcm' + os.sep + \
                     'defaultpcm.ini'

LONGOPTS = ["inifile", 'outputdir', 'write-ini']
SHORTOPTS = "i:o" 

#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------

class PCMConfigManager(ConfigManager):

    __instance = {}
    
    def get_instance(inifile, realdatastr=None):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        The configuration can be obtained with a filename or with the original
        string given to get_instance, i.e., file '/home/me/myinifile.ini' can 
        be obtained with '/home/me/myinifile.ini' and 'myinifile.ini'
        """
        inifilebase = os.path.basename(inifile)
        if inifile in PCMConfigManager.__instance:
            return PCMConfigManager.__instance[inifile]
        elif inifilebase  in PCMConfigManager.__instance:
            return PCMConfigManager.__instance[inifilebase]
        else:
            PCMConfigManager.__instance[inifile] = None # Set for __init__()
            PCMConfigManager.__instance[inifile] = PCMConfigManager(inifile,
                                                                realdatastr)
            PCMConfigManager.__instance[inifilebase] = PCMConfigManager.__instance[inifile]
        return PCMConfigManager.__instance[inifile]
    get_instance = staticmethod(get_instance)

    def delete(inifile):
        """remove the instance inifile from the loaded configurations"""
        del PCMConfigManager.__instance[inifile]
        if inifile != os.path.basename(inifile):
            del PCMConfigManager.__instance[os.path.basename(inifile)]
    delete = staticmethod(delete)

    def __init__(self, filename = INIFILE_DEFAULT, realdatastr=None):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if (filename not in PCMConfigManager.__instance) or (
                PCMConfigManager.__instance[filename] is not None):
            raise Exception("This class is a singleton per filename. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename, realdatastr)

    def register_defaults(self):
        """default ini settings for a fiber1d problem"""
        self.register("general.read", False)
        self.register("general.verbose", False)

        self.register("pcm.radius", 0.01,
            "radius of the pcm in mm")
        self.register("pcm.melting_point", 28.24,
            "Melting point of the pcm in degree Celcius")
        self.register("pcm.latent_heat_fusion", 238.76,
            "Amount of heat that can be absorbed during melting, in kJ/kg.")
        self.register("pcm.density_solid", 779,
            "Density of the solid for normal conditions (20 degree Celcius), in kg/m^3")
        self.register("pcm.specific_heat_solid", 1.9,
            "Specific heat of the solid for normal conditions (25 degree Celcius), in kJ/(kg K)")
        self.register("pcm.specific_heat_liquid", 2.1,
            "Specific heat of the liquid for normal conditions (35 degree Celcius), in kJ/(kg K)")
        self.register("pcm.thermal_cond_solid", 0.4,
            "Thermal conductivity of the solid, in W/(m K)")
        self.register("pcm.thermal_cond_liquid", 0.3,
            "Thermal conductivity of the liquid, in W/(m K)")

        self.register("discretization.n_edge", 41,
            "Numerical method needs a grid discritization, here number of edges")
    
        self.register("init.init_temp", 'lambda x: 26.',
            "Initial temperature of the pcm, in degree Celcius")

        self.register("boundary.heat_transfer_coeff", 8.5,
            "heat transfer coefficient through pcm boundary, in W/(m^2 K)")
        self.register("boundary.T_out", 31.5,
            "Outside temperature, in degree Celcius")
        
        self.register("fabric.simulate", False,
            "simulate effect of PCM on a 1D depth profile of a fabric, True or False")
        self.register("fabric.thickness", 2.,
            "Thickness of the fabric in mm")
        self.register("fabric.PCM_fraction_distribution", "lambda x: 10.",
            "In % of a REV, how many PCM volume there is, as a function of "
            "position over the thickness")
        self.register("fabric.porosity", 0.88,
            "Porosity of the fabric without PCM, in % of REV volume")
        self.register("fabric.vol_heat_capacity", 1847.6,
            "Volumetric heat capaity of fabric, in kJ/(m^3 K)")
        self.register("fabric.therm_cond", 0.0404,
            "Thermal conductivity fabric (K_{mix}), in W/(m K)")
        
        self.register("time.time_period", 600.)
        self.register("time.dt", 5.0)

        #plot section
        self.register("plot.plotavgtemp", True, 
            'plot average temperature of the PCM over time')
        self.register("plot.plotinterface", True,
            'plot position interface liquid/solid in PCM over time')
        self.register("plot.plotevery", 10,
            "When plotting over time, indicate how many steps dt to skip before plotting again."
            " If 0, no plot occurs")
