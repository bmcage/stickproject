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
import os.path
import sys
import const
import lib.utils.utils as utils
from model.diffmodel import DiffusionModel
import lib.utils.gridutils as GridUtils
import config as conf
import lib.diff.diffusion as diffusion

#-------------------------------------------------------------------------
#
# DiffusionModel class
#
#-------------------------------------------------------------------------
class FiberModel(DiffusionModel):
    """
    FiberModel is the generic diffusion model for fibers. It has knowledge of 
    the    implemented models, parses the ini file and passes responsibility to the
    correct diffusionmodel implementation
    """
    def __init__(self, config):
        """
        a config class must be passed in that contains the required settings
        """
        DiffusionModel.__init__(self)
        
        self.datatime = []

    def _create_compgrid(self):
        """
        Construct the computational grid we need to work on as a 
        grid over [0,1]
        """
        ## TODO
        # make mesh on [0,1]
        self.__grid = GridUtils.create_grid(
                        self.cfg.get('discr.points'), 0., 1., "borderrefined", 
                        1)

    def _read_init_cond(self):
        """
        Read the initial condition, determine cutoff and length of sample 
        from it, and project on the grid to obtain the init cond over the
        grid
        """
        if len(self.datatime)> 0:
            print 'ERROR: initial condition already read in, do not do it twice'
            sys.exit()
        comps = self.cfg.get('general.components')
        if self.cfg.get('init.initfromfile'):
            starttime, meshinitpos, components, conc_components = \
                    GridUtils.read_exp_file(utils.add_root_to_file(
                                                self.cfg.get('init.initfile')))
            self.datatime.append(starttime)
            self.datamesh.append(meshinitpos)
            dataval = []
            for comp in comps:
                #select this from component from the experiment
                try:
                    dataval.append(conc_components[components[comp]])
                except KeyError:
                    print 'ERROR: component %s not present in initfile' % comp
                    sys.exit()
            self.data.append(dataval)
        else :
            #read from function
            initfunc = eval(self.cfg.get('init.initfunc'))
            startgrid = self.cfg.get('init.initfuncstart')
            endgrid = self.cfg.get('init.initfuncend')
            starttime = 0.
            nrc = len(initfunc(pp))
            data = []
            for i in range(nrc):
                data.append([])
            data = tuple(data)
            self.datamesh.append(self.__grid * (endgrid-startgrid) + startgrid)
            for p in self.__grid : 
                pp = p * (endgrid-startgrid) + startgrid
                conc = initfunc(pp)
                ind = 0
                for val in conc:
                    data[ind].append(val)
                    ind += 1
            self.datatime.append(starttime)
            self.data.append(data)
            self.left_cutoff.append(startgrid)
            self.lenght.append(endgrid-startgrid)

    def _read_experiments(self):
        """
        read in experiments and add to data
        """
        if not len(self.data) > 0:
            print 'ERROR: Set Initial condition first'
            sys.exit()
        nrexp = self.cfg.get('experiments.nrexp')
        comps = self.cfg.get('general.components')
        for exp in range(nrexp):
            file = utils.add_root_to_file(
                        self.cfg.get('experiment_%i.file' % exp))
            time_exp, meshexp, components, conc_components = \
                    GridUtils.read_exp_file(file)
            if (meshexp[0] < self.datamesh[0][0] \
                or meshexp[-1] > self.datamesh[0][-1]):
                print ('ERROR: Experiment files contain data that conflicts'
                            ' with data in init file')
                sys.exit()
            if time_exp <= self.datatime[-1] : 
                print ('ERROR: Time in experiment file conflicts earlier '
                            'given times')
                sys.exit()
        
            self.datatime.append(time_exp)
            self.datamesh.append(meshexp)
            dataval = []
            for comp in comps:
                #select this from component from the experiment
                try:
                    dataval.append(conc_components[components[comp]])
                except KeyError:
                    print 'ERROR: component %s not present in initfile' % comp
                    sys.exit()
            self.data.append(dataval)

    def _setup_data_on_compgrid(self):
        """
        grid and data should be constructed/read. We interpolate data on the
        computational grid. These are numpy arrays
        """
        ## TODO
        for mesh, data, coff, leng in zip(self.datamesh, self.data, 
                                          self.left_cutoff, self.lenght):
            datagrids = []
            for component in data:
                datagrids.append(GridUtils.project_exp_to_grid(
                                        mesh, component,
                                        self.__grid, coff, leng))
            self.datagrid.append(datagrids)

    def _setup_model(self):
        ## TODO
        dtype =  self.cfg.get('diffusion.type')
        method = self.cfg.get('general.method')
        if not (method == diffitconfig.FINITE_VOLUME_ML):
            print 'Only method = 0 (=FINITE_VOLUME_ML) is supported'
            sys.exit(1)
        submethod = self.cfg.get('general.submethod')
        if submethod == diffitconfig.CVODE:
            subm = 'cvode'
        elif submethod == diffitconfig.ODEINT:
            subm = 'odeint'
        else:
            print 'submethod given not supported, give CVODE (0) or ODEINT (1)'
            sys.exit(1)
        if dtype == 'const':
            #linear diffusion model, constant diffusion
            self.diff = diffusion.Diffusion_const(self.cfg.get(
                            'diffusion.par')[0])
            nrcomp = self.cfg.get('general.nrcomponents')
            if nrcomp == 1:
                #one comp linear diffusion
                self.component = self.cfg.get('data.component')[0]
                from lindiff import Lin1DDiffusion
                self.model = Lin1DDiffusion(
                                    self.__grid,
                                    self.datagrid[0][self.component],
                                    self.diff,
                                    length_grid=self.lenght[0],
                                    left_cutoff=self.left_cutoff[0],
                                    solmeth='conservative_fv')
            else:
                raise NotImplementedError
        elif dtype == 'D(u)_lin_interp':
            vals = self.cfg.get('diffusion.par')
            self.diff = diffusion.Diffusion_u(vals[0], vals[1], '1D2pint')
            nrcomp = self.cfg.get('general.nrcomponents')
            if nrcomp == 1:
                #one comp nonlinear diffusion
                self.component = self.cfg.get('data.component')[0]
                from nonlindiff import Nonlin1DDiffusion
                self.model = Nonlin1DDiffusion(
                                    self.__grid,
                                    self.datagrid[0][self.component],
                                    self.diff,
                                    length_grid=self.lenght[0],
                                    left_cutoff=self.left_cutoff[0],
                                    solmeth='conservative_fv')
            else:
                raise NotImplementedError
        elif dtype == 'D(u,v)_lin_interp':
            vals = self.cfg.get('diffusion.par')
            self.diff = diffusion.Diffusion_u_MC(vals[0], vals[1], 
                                                 vals[2], vals[3])
            nrcomp = self.cfg.get('general.nrcomponents')
            if nrcomp == 2:
                #multi comp nonlinear diffusion
                from nonlindiffmc import Nonlin1DDiffusion_MC
                self.model = Nonlin1DDiffusion_MC(
                                    self.__grid,
                                    self.datagrid[0][0], self.datagrid[0][1],
                                    self.diff,
                                    length_grid=self.lenght[0],
                                    left_cutoff=self.left_cutoff[0],
                                    solmeth='conservative_fv', 
                                    submeth=subm,
                                    jac=False)
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError
