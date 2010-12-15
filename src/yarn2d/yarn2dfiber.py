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
        self.type_fiber = self.cfg.get('fiber.number_type')
        self.blend = self.cfg.get('fiber.blend')
        self.verbose = self.cfg.get('general.verbose')
        
    def create_fiber_kinds(self, filename = 'determine_kinds.dat', filename1 = 'index_fiber.dat',
                    regenerate = True):
        """
        create the file to determine the fiber kinds in the yarn and the fibers in it
        returns string with defenition, path to file
        """
        filepath = utils.OUTPUTDIR + os.sep + filename
        filepath1 = utils.OUTPUTDIR + os.sep + filename1
        if regenerate:
            self.determine_file = open(filepath, 'w')
            index_1 = 0.0
            index_2 = 0.0 
            index_3 = 0.0
            
            for type_fiber_i in sp.arange(self.type_fiber):
                fiber_kind = self.blend[type_fiber_i] * self.number_fiber / 100.
                index_1 = 0
                while index_1 < fiber_kind:
                    determine_a = np.random.uniform(0, 1)
                    if determine_a > (type_fiber_i ) / self.type_fiber and determine_a <= (type_fiber_i + 1) / self.type_fiber:
                        self.kind = float(type_fiber_i)
                        index_1 += 1.0
                        self.determine_file.write("%s\n" %(repr(self.kind)))
            index_fiber = sp.empty(self.number_fiber, int)
            index_i = 0
            self.index_file = open(filepath1, 'w')
            print 'the number of fiber', self.number_fiber
            while index_i < self.number_fiber:
                init_index = np.random.uniform(1, self.number_fiber + 1)
                if index_i == 0:
                    index_fiber[index_i] = int(init_index)
                    print 'this is index number', index_fiber[index_i]
                    self.index_file.write("%s\n" %(repr(index_fiber[index_i])))
                    index_i += 1
                else:
                    determine_index = (int(init_index) == index_fiber)
                    if determine_index.any() != True:
                        index_fiber[index_i] = int(init_index)
                        self.index_file.write("%s\n" %(repr(index_fiber[index_i])))
                        index_i += 1
            self.determine_file.close()
            self.index_file.close()
        else:
            deterfile = open(filepath, 'r')
            #self.deterfile_value = eval(deterfile.read())
            deterfile.close()