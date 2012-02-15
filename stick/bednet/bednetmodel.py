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
    upscaling to fabric1d domain
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
from scipy import integrate
import matplotlib.pyplot as plt
import sets
import time
import math
from mpmath import jtheta

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import lib.utils.utils as utils
import lib.utils.gridutils as GridUtils
import bednet.config as conf
from yarn.config import YarnConfigManager
from yarn1d.yarn1dmodel import Yarn1DModel

#-------------------------------------------------------------------------
#
# DiffusionModel-fabric1d class 
#
#-------------------------------------------------------------------------
class Bednet(object):
    """
    Fabirc2dUpscaling is a special diffusion model to calculate the concentration of 
    a volatile outside the fabric. This upscaling method is simple with the equation:
    C_outside = C_yarn * R_yarn/d_distance
    """
    def __init__(self, config):
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.timesteps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.timesteps + 1)
        #set correct delta_t
        self.delta_t = self.times[1]-self.times[0]
        if self.verbose:
            print "Timestep used in bednet model:", self.delta_t

        #self.n = self.cfg.get('domain.nr_vert_yarns')
        #self.m = self.cfg.get('domain.nr_hor_yarns')   
        self.domain_size = self.cfg.get('domain.domain_size')
        self.dx = self.cfg.get('domain.dx')
        self.dy  = self.cfg.get('domain.dy')
        self.size_sample = self.cfg.get('sample.size_sample')
        #self.thickness_sample = self.cfg.get('sample.thickness_sample')
        self.boundary_up = self.cfg.get('boundary.boundary_up')
        self.boundary_bottom = self.cfg.get('boundary.boundary_bottom')
        self.boundary_left = self.cfg.get('boundary.boundary_left')
        self.boundary_right = self.cfg.get('boundary.boundary_right')
        self.diff_coef = self.cfg.get('diffusion.diff_coef')
        self.saturation_conc = self.cfg.get('saturation.saturation_conc')
        x0 = self.cfg.get('observer.x0')
        self.x0 = sp.empty(len(x0)+1, float)
        self.x0[1:] = x0[:]
        """
        self.distance_yarn = self.cfg.get('domain.distance_yarn')
        self.grid_space_vertical = self.cfg.get('domain.grid_space_vertical')
        self.number_nodes = (self.domain_size[0] / self.distance_yarn + 1) * \
                            (self.domain_size[1] / self.grid_space_vertical + 1)
        """
        #we set a distance for the yarn bc
        self.x0[0] = 0.1
        self.cfg_yarn = []
        for filename in self.cfg.get('sample.yarn_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                        os.path.dirname(self.cfg.filename), filename))
            self.cfg_yarn.append(YarnConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            print 'time', self.time_period
            self.cfg_yarn[-1].set("time.time_period", self.time_period)
            self.cfg_yarn[-1].set("boundary.dist_conc_out", float(self.x0[0]))
            
        #create yarn models
        self.yarn_models = []
        for cfg in self.cfg_yarn:
            self.yarn_models.append(Yarn1DModel(cfg))
        self.nr_models = len(self.yarn_models)

        #some memory
        self.source_mass = np.empty((self.nr_models, self.timesteps + 1), float)
        
        #plot the result every few seconds so outcome becomes visible during calculations
        self.plotevery = self.cfg.get("plot.plotevery")

    def initial_void_conc(self):
        """ initial concentration over the domain"""
        self.init_void = self.cfg.get('initial.init_void')
        if self.plotevery:
            self.viewerplotcount = 1
        
    def init_yarn(self):
        self.yarn_mass = [0] * len(self.yarn_models)
        self.tstep = 0
        for ind, model in enumerate(self.yarn_models):
            model.do_yarn_init()
            self.yarn_mass[ind] = model.calc_mass()
            # no mass released at start time
            self.source_mass[ind, self.tstep] = 0
    
    def initial_boundary_conc(self):
        self.initconc = self.cfg.get('initial.init_conc')
        self.initvoidconc = self.cfg.get('initial.init_void')

    def solve_timestep(self, t):
        print "solve timestep", t
        self.tstep += 1
        # 1. step one, solve the yarn model
        for ttype, model in enumerate(self.yarn_models):
            rt, rety = model.do_yarn_step(t)
            tmp = model.calc_mass()
            self.source_mass[ttype, self.tstep] = self.yarn_mass[ttype] - tmp
            self.yarn_mass[ttype] = tmp
        
        # 2. step two, solve the bednet model
        #    to obtain value near yarn, and at observer
        #    we know that self.source_mass[ttype] has been released since
        #    last step
        x0 = self.x0
        # concentration is a consequence of all previous releases, so sum 
        # over all times, and compute contribution of that moment.
        self.sol[self.tstep, :] = 0.
        for tt in np.arange(self.tstep):
            release = tt * self.delta_t
            timesincerelease = t-release
            factor = 1 / (4 * self.diff_coef * timesincerelease)
            for ttype in np.arange(self.nr_models):
                # TODO: this should be dx/dy per yarn ttype!!
                AH = math.pow(self.dx, 2) * factor
                AW = math.pow(self.dy, 2) * factor
                res = (
                         self.source_mass[ttype, tt+1] *factor / (2 * sp.pi)
                         * sp.exp(-sp.power(x0[:], 2) * factor) * 
                         (jtheta(3, 0, sp.exp(-AH)) 
                          + jtheta(3, 0, sp.exp(-AW)) + 2)
                        )
                for ind, val in enumerate(res):
                    self.sol[self.tstep, ind] += float(val)
        # 3. for next timestep, we need to set correct boundary condition
        #    on the yarn level
        for ind, model in enumerate(self.yarn_models):
            model.boundary_conc_out = self.sol[self.tstep, 0]

    def view_sol(self):
        #maxv = np.max(self.sol)
        #minv = np.min(self.sol)
        #print 'max', maxv, minv
        plt.ion()
        for ind, pos in enumerate(self.x0[1:]):
            plt.figure(ind)
            plt.plot(self.times, self.sol[:, ind+1])
            #plt.ylim(0, maxv*1.1)
            plt.title('Concentration at position %g' % pos)
            plt.show()

    def init_bednet(self):
        self.sol = sp.empty((self.timesteps+1, len(self.x0)), float)
        self.sol[0, :] = 0
        self.initial_boundary_conc()
        self.init_yarn()

    def run(self, wait=False):
        self.init_bednet()
        for t in self.times[1:]:
            self.solve_timestep(t)

        self.view_sol()

        if wait:
            raw_input("Finished bednet run")
