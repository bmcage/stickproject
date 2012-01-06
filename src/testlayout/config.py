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
This package implements config defaults for test layout
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
INIFILE_DEFAULT = const.INI_DIR + os.sep + 'test' + os.sep + \
                    'defaulttest.ini'

LONGOPTS = ["inifile", 'outputdir', 'mesh', 'write-ini']
SHORTOPTS = "i:o" 

METHOD = {
    'FVM': ('Finite Volume Method discretization', ['fipy']),
    }

FIBERLAYOUTS = {
    'random': ('random position of fibers in the yarn',),
    'virtloc': ('default virtual locations (start in center, no overlap)',),
    'virtlocoverlap': ('adapted virtual locations (different fiber size, overlap)',),
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
class TestFiberLayoutManager(ConfigManager):
    
    __instance = None
    
    def get_instance(inifile):
        """
        Use this function to get the instance of the ConfigManager
        that will work on inifile.
        """
        if TestFiberLayoutManager.__instance is None:
            TestFiberLayoutManager.__instance = 1
            TestFiberLayoutManager.__instance = TestFiberLayoutManager(inifile)
        return TestFiberLayoutManager.__instance
    get_instance = staticmethod(get_instance)
    
    def __int__(self, filename = INIFILE_DEFAULT):
        """
        A singleleton implementation of config.ConfigManager
        """
        if TestFiberLayoutManager.__instance is not 1:
            raise Exception("This class is a singleton."
                            "Use the get_instance() method")
        ConfigManager.__init__(self, filename)
        
    def register_defaults(self):
        """
        default ini settings for a test problem
        """
        #general section
        self.register("general.read", False, 
                    "only set to False")
        self.register("general.verbose", False)
        
        self.register("information.domain_radius", 0.6)
        print 'read the size of domain'
        self.register("information.each_circle_size", 0.15)
        
        self.register("position.domain_central", [0., 0.])
        self.register("position.circle_1", [0., 0.])
        self.register("position.circle_2", [0., 0.])
        self.register("position.circle_3", [0., 0.])
        #self.register("position.circle_4", [0., 0.])
        
        