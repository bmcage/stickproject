#
# Copyright (C) 2010  B. Malengier
# Copyright (C) 2010  P.Li
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

""" 
Module for functions of determine the fiber kind. 
"""
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
from __future__ import division
import os.path
import sys
import const
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
from fipy import Gmsh2D

#-------------------------------------------------------------------------
#
# Yarn2dFiber class
#
#-------------------------------------------------------------------------

class Yarn2dFiber(object):
    def __init__(self, cfg):
        self.cfg = cfg
        
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        self.blend = self.cfg.get('fiber.blend')
        self.verbose = self.cfg.get('general.verbose')
        
    def create_fiber_kinds(self, filename = 'determine_kinds.dat', regenerate = True):
        """
        create the file to determine the fiber kinds in the yarn and the fibers in it
        returns string with defenition, path to file
        """
        filepath = utils.OUTPUTDIR + os.sep + filename
        if regenerate:
            self.determine_file = open(filepath, 'w')
            index_1 = 0.0 #index for the loop to generate the number of ifber
            index_2 = 0.0 #index for the loop to generate the cotton
            index_3 = 0.0 #index for the loop to generate the polyester
            fiber_kind1 = self.blend[0] * self.number_fiber / 100.
            fiber_kind2 = self.blend[1] * self.number_fiber / 100.
            self.kind1 = sp.empty(fiber_kind1)
            self.kind2 = sp.empty(fiber_kind2)
            while index_2 < fiber_kind1 or index_3 < fiber_kind2:
                determine_a = np.random.uniform(0, 1)
                if determine_a <= 0.5 and index_2 < fiber_kind1:
                    self.kind1 = 0.0
                    index_2 += 1.0
                    #self.determine_file.write(self.kind1)
                    self.determine_file.write("%s\n" %(repr(self.kind1)))
                elif determine_a > 0.5 and index_3 < fiber_kind2:
                    self.kind2 = 1.0
                    index_3 += 1.0
                    #self.determine_file.write(self.kind2)
                    self.determine_file.write("%s\n" %(repr(self.kind2)))
            
            self.determine_file.close()
            """
            deterfile = open(filepath, 'r')
            self.deterfile_value = eval(deterfile.read())
            #print deterfile_value
            deterfile.close()
            """
        else:
            deterfile = open(filepath, 'r')
            self.deterfile_value = eval(deterfile.read())
            deterfile.close()
            