"""
To test whether the module: arearatioprobability.py within the function working
well
"""
from __future__ import division
import os.path
import sys
import const
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import pylab
import matplotlib

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
import yarn2d.arearatioprobability as arearp
import testlayout.config as conf


class TestLayoutCalculation(object):
    """
    TestLayoutCalculation is the special class to confirm whether the layout 
    calculation is right
    """
    def __init__(self, config):
        self.datatime = []
        self.cfg = config
        
        self.verbose = self.cfg.get('general.verbose')
        self.domain_radius = self.cfg.get('information.domain_radius')
        self.each_circle_size = self.cfg.get('information.each_circle_size')
        print 'value of self.each_circle_size', self.each_circle_size
        self.domain_central = self.cfg.get('position.domain_central')
        self.circle_1 = self.cfg.get('position.circle_1')
        self.circle_2 = self.cfg.get('position.circle_2')
        self.circle_3 = self.cfg.get('position.circle_3')
        #self.circle_4 = self.cfg.get('position.circle_4')
        
        
    def check_ratio_calculation(self):
        check_yarn = self.domain_radius
        check_fiber = sp.zeros(3, float)
        check_x = sp.zeros(3, float)
        check_y = sp.zeros(3, float)
        check_fiber[:] = self.each_circle_size
        
        check_central_x = self.domain_central[0]
        check_central_y = self.domain_central[1]
        
        check_x[0] = self.circle_1[0]
        check_x[1] = self.circle_2[0]
        check_x[2] = self.circle_3[0]
        #check_x[3] = self.circle_4[0]
        
        check_y[0] = self.circle_1[1]
        check_y[1] = self.circle_2[1]
        check_y[2] = self.circle_3[1]
        #check_y[3] = self.circle_4[1]
        
        check_position, ratio_check = arearp.calculate_proportion(check_yarn, 
                                    check_fiber, check_x, check_y)
        print 'ratio_check value is', ratio_check
        
        return (ratio_check)
        
    def run(self):
        self.check_ratio_calculation()