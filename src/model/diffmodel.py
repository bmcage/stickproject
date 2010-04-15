#!/usr/bin env python

# Copyright (C) 2009-2010  B. Malengier
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

""" module holding a generic diffusion model. 
"""
#-------------------------------------------------------------------------
#
# Imports
#
#-------------------------------------------------------------------------

import sys
import os.path
import const
import lib.utils.utils as utils
from lib.utils.utilstime import TimeParameters

#-------------------------------------------------------------------------
#
# DiffusionModel class
#
#-------------------------------------------------------------------------
class DiffusionModel(object):
    """
    DiffusionModel is the generic framework to run a diffusion model. 
    It has knowledge of the
    implemented models, parses the ini file and passes responsibility to the
    correct diffusionmodel implementation.
    Inherit this method for real functionality
    """
    def __init__(self, config):
        """
        a config class must be passed in that contains the required settings
        """
        self.cfg = config
        self.model = None      # the actual model
        self._grid = None      # the grid on which to work
        self.datatime = []     # time on which there is data, [0] = init val
        self.datamesh = []     # mesh over which there is data, [0] = init val
        self.dataval  = []     # value of the data, [0] = init val
        self.diff = None       # the diffusion used
        self.component = []    # several components can be present, which ones
                               #   to track

    def _create_compgrid(self):
        """
        Construct the computational grid we need to work on 
        """
        raise NotImplementedError

    def _read_init_cond(self):
        """
        read the initial condition in
        """
        raise NotImplementedError

    def _read_experiments(self):
        """
        read in experiments
        """
        raise NotImplementedError

    def _setup_data_on_compgrid(self):
        """
        grid and data should be constructed/read. We interpolate data on the
        computational grid. These are numpy arrays
        """
        raise NotImplementedError

    def _setup_model(self):
        """
        set self.model as a diffusion problem
        """
        raise NotImplementedError

    def setup_for_solve(self):
        """
        setup all data as required before solving
        """
        self.__create_compgrid()
        self.__read_init_cond()
        self.__read_experiments()
        self.__setup_data_on_compgrid()
        self.__setup_model()
        self.model.set_boundary_cond(self.cfg)
        timepar = TimeParameters(self.cfg, exptime=self.datatime)
        self.model.set_time(timepar)
        self.model.set_experiment_data(self.datatime, self.datamesh,
                            self.data, self.component)

    def solve(self):
        """
        Solve the underlying model. It is required that the model is in a 
        state that allows it to be solved
        """
        self.model.solve()
        
    def run(self):
        """
        starts the diffusion model
        """
        self.setup_for_solve()
        self.solve()
        self.model.writeinfo(utils.OUTPUTDIR)
        self.model.showmass(utils.OUTPUTDIR)
        self.model.plot_sol(utils.OUTPUTDIR)
