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

from __future__ import division
import os

#---------------------------------------------------------------
#
# local imports
#
#---------------------------------------------------------------
import stick.const as const
from stick.lib.config import ConfigManager

#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'yarn' + os.sep + \
                     'defaultyarn.ini'

LONGOPTS = ["inifile", 'outputdir', 'mesh', 'write-ini']
SHORTOPTS = "i:o" 

METHOD = {
    'FVM': ('Finite Volume Method discretization', ['fipy']),
    }

FIBERLAYOUTS = {
    'uniform': ('the porosity is the same everywhere'),
    'virtlocoverlap': ('adapted virtual locations (different fiber size, overlap), \
                    area distribution is required',),
    'random': ('do not use: random position of fibers in the yarn',),
    'virtloc': ('do not use: default virtual locations (start in center, no overlap)',),
    }

DIFF_FLUX = 0
TRANSFER = 1
BOUND_TYPE = {
    'diff_flux': DIFF_FLUX,
    'transfer': TRANSFER,
    }
    
#---------------------------------------------------------------
#
# DiffitConfigManager class
#
#---------------------------------------------------------------

class YarnConfigManager(ConfigManager):

    __instance = {}
    
    def get_instance(inifile, realdatastr=None):
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        The configuration can be obtained with a filename or with the original
        string given to get_instance, i.e., file '/home/me/myinifile.ini' can 
        be obtained with '/home/me/myinifile.ini' and 'myinifile.ini'
        """
        inifilebase = os.path.basename(inifile)
        if inifile in YarnConfigManager.__instance:
            return YarnConfigManager.__instance[inifile]
        elif inifilebase  in YarnConfigManager.__instance:
            return YarnConfigManager.__instance[inifilebase]
        else:
            YarnConfigManager.__instance[inifile] = None # Set to 1 for __init__()
            YarnConfigManager.__instance[inifile] = YarnConfigManager(inifile,
                                                        realdatastr)
            YarnConfigManager.__instance[inifilebase] = YarnConfigManager.__instance[inifile]
        return YarnConfigManager.__instance[inifile]
    get_instance = staticmethod(get_instance)

    def delete(inifile):
        """remove the instance inifile from the loaded configurations"""
        del YarnConfigManager.__instance[inifile]
        if inifile != os.path.basename(inifile):
            del YarnConfigManager.__instance[os.path.basename(inifile)]
    delete = staticmethod(delete)

    def __init__(self, filename = INIFILE_DEFAULT, realdatastr=None):
        """ 
        A singleton implementation of config.ConfigManager
        """
        if (filename not in YarnConfigManager.__instance) or (
                YarnConfigManager.__instance[filename] is not None):
            raise Exception("This class is a singleton per filename. "
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename, realdatastr)

    def register_defaults(self):
        """default ini settings for a DiffusionIT problem"""
        #general section
        self.register("general.read", False,
            "Set True to read fiber-yarn layout data from the output"
            " directory instead of regenerating")
        self.register("general.verbose", False)
        self.register("general.method", 'FVM', 
            "Solution method to use, one of " + ",".join(METHOD.keys()))
        self.register("general.submethod", 'fipy', 
            "Solution submethod to use")

        #domain section
        self.register("domain.distribute_fiber", 'integral',
            "distribute fiber to the yarn layout")
        self.register("domain.cellsize_centre", 5.0e-2, 
            "preferred edge length of each mesh for yarn")
        self.register("domain.cellsize_fiber", 5.0e-2,
            "preferred edge length of each mesh for fiber")
        self.register("domain.yarnradius", 1.,
            "radius of yarn domain in mm")
        self.register("domain.useextension", False,
            "extend the yarn domain with an extension outside the yarn")
        self.register("domain.extensionfraction", 1.,
            "extend the domain with this fraction*yarnradius")
        self.register("domain.fiberlayout_method", 'uniform',
            "Method for the fiber layout, one of " 
            + ",".join(FIBERLAYOUTS.keys()))
        self.register("domain.n_edge", 10, 'For 1D modelling, how many edges '
            'the domain mesh has.')
        self.register("domain.theta_value", 0.05)
        self.register("domain.beta_value", 0.04)
        self.register("domain.radius_first_center_virtloc", 0.)
        #fiber section
        self.register("fiber.number_type", 1, 
            "Number of fiber types present")
        self.register("fiber.number_fiber", 60, 
            "Total number of fibers in the yarn")
        self.register("fiber.eps_value", 0.001)
        self.register("fiber.blend", [100.0], 
            "Blend distribution over the fiber types")
        self.register("fiber.radius_fiber", [100.])
        self.register("fiber.fiber_config", ['../fiber/defaultfiber.ini'])
        self.register("fiber.prob_area", '[lambda r: 0.2]')

        #initial section
        self.register("initial.init_conc1d", 'lambda r: 0.', 
            "Initial concentration of tracked compound in the yarn in terms of"
            " radius in 1D")
        self.register("initial.init_conc2d", 'lambda r, theta: 0.', 
            "Initial concentration of tracked compound in the yarn in terms of"
            " radius, angle in 2D")
        self.register('initial.init_extension', 0., 
            'Initial concentration of tracked compound in the extension region')
        
        #diffusion
        self.register("diffusion.diffusion_coeff", 25.,
            "Diffusion coefficient of tracked compound in the yarn in mm2/s")
        
        #boundary section
        self.register("boundary.type_right", 'diff_flux', 
            "boundary type at surface yarn: transfer or diff_flux")
        self.register("boundary.conc_out", 0.,
            "outside concentration, so if type diff_flux, "
            "flux = - D_out * (conc_out - yarn_edge_conc)/dist_conc_out")
        self.register("boundary.D_out", 7.77,
            "outside diffusion coef [unit mm2/s], so if type diff_flux, "
            "flux = - D_out * (conc_out - yarn_edge_conc)/dist_conc_out")
        self.register("boundary.dist_conc_out", 0.1,
            "distance of conc_out to yarn edge in mm, so if type diff_flux, "
            "flux = - D_out * (conc_out - yarn_edge_conc)/dist_conc_out")
        self.register("boundary.transfer_coef", 5.3e-3,
            "tracked compound exterior transfer coef, so if type transfer, "
            "flux - D dC/dx = + transfer_coef C")

        #time section
        self.register("time.time_period", 4000.,
            "the time domain for the simulation (s)")
        self.register("time.dt", 1.0,
            "the time step to use in the simulation")

        #plot section
        self.register("plot.maxval", 0.0005,
            "When plotting tracked compound, set a max value for vertical axis")
        self.register("plot.plotevery", 10,
            "When plotting over time, indicate how many steps dt to skip before plotting again")
        self.register("plot.writeevery", 0,
            "When writing data out over time, indicate how many steps dt to skip")
