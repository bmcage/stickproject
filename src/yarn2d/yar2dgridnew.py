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
Module for functions of a yarn 2D grid. 
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
# Yarn2dGrid new class
#
#-------------------------------------------------------------------------

class Yarn2dNewGrid(object):
    def __init__(self, cfg):
        self.cfg = cfg
        
        self.Ry = self.cfg.get('domain.yarnradius')
        self.Rf = self.cfg.get('fiber.radius_fiber')
        self.number_circle = self.cfg.get('yarn.circle_num')
        self.angle_virtual = self.cfg.get('yarn.angel_virtual')
        
        self.radius_yarn = self.scaleL * self.Ry
        self.radius_fiber =  [self.scaleL * rad for rad in self.Rf]
        self.radius_boundlayer = max(self.radius_fiber)/2.
        self.radius_domain = self.radius_yarn + self.radius_boundlayer
        self.cellsize_centre = self.cfg.get('domain.cellsize_centre')
        self.cellSize = self.cfg.get('domain.cellsize_fiber')
        self.number_fiber = self.cfg.get('fiber.number_fiber')
        
        self.verbose = self.cfg.get('general.verbose')
        
    def create_circle_domain_gmesh(self, filename='yarn_new.geo', filename1 = 'fib_centers_x_new.geo',
        filename2 = 'fib_centers_y_new.geo' ,regenerate=True):
        filepath = utils.OUTPUTDIR + os.sep + filename
        filepath1 = utils.OUTPUTDIR + os.sep + filename1
        filepath2 = utils.OUTPUTDIR + os.sep + filename2
        self.number_for_circles = int((self.radius_domain - self.radius_fiber) / (2 * \
                                    self.radius_fiber))
        self.total_circles = self.number_for_circles + 1
        
        if regenerate:
            
            
            
        
        
        