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
from scipy.integrate import quad
import matplotlib.pyplot as plt
import pylab
import sets
import time
import math
import fipy

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
        self.times = np.linspace(0, self.time_period, self.timesteps + 1)
        #set correct delta_t
        self.delta_t = self.times[1]-self.times[0]
        if self.verbose:
            print "Timestep used in bednet model:", self.delta_t
        self.initconc = self.cfg.get('initial.init_conc')
        #self.n = self.cfg.get('domain.nr_vert_yarns')
        #self.m = self.cfg.get('domain.nr_hor_yarns')   
        self.domain_size = self.cfg.get('domain.domain_size')
        self.dx = self.cfg.get('domain.dx')
        self.dy  = self.cfg.get('domain.dy')
        self.nvertyarns = self.cfg.get('domain.nr_vert_yarns')
        self.nhoryarns = self.cfg.get('domain.nr_hor_yarns')
        #self.thickness_sample = self.cfg.get('sample.thickness_sample')
        self.boundary_up = self.cfg.get('boundary.boundary_up')
        self.boundary_bottom = self.cfg.get('boundary.boundary_bottom')
        self.boundary_left = self.cfg.get('boundary.boundary_left')
        self.boundary_right = self.cfg.get('boundary.boundary_right')
        self.diff_coef = self.cfg.get('diffusion.diff_coef')
        self.saturation_conc = self.cfg.get('active_component.saturation_conc')
        self.treshold = self.cfg.get('active_component.treshold_effect')
        x0 = self.cfg.get('observer.x0')
        self.x0 = np.empty(len(x0)+1, float)
        self.x0[1:] = x0[:]
        
        #we set a distance for the yarn bc
        EXTFRAC = 1.
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
            self.cfg_yarn[-1].set("boundary.D_out", self.diff_coef)
            self.cfg_yarn[-1].set("boundary.conc_out", self.initconc)
            self.cfg_yarn[-1].set("domain.useextension", True)
            ##TODO How much overlap region? Take one yarn radius for now
            self.cfg_yarn[-1].set("domain.extensionfraction", EXTFRAC)
        radyarn = self.cfg_yarn[-1].get("domain.yarnradius")
        for conf in self.cfg_yarn:
            assert radyarn == self.cfg_yarn[-1].get("domain.yarnradius")
            
        # we need concentration at half of the extra zone
        self.x0[0] = radyarn + EXTFRAC*radyarn/2
        self.overlapsize = EXTFRAC*radyarn
        self.overlaparea = np.pi * (np.power(radyarn+self.overlapsize,2)
                                    - np.power(radyarn,2))
        
        #create yarn models
        self.yarn_models = []
        for cfg in self.cfg_yarn:
            self.yarn_models.append(Yarn1DModel(cfg))
        self.nr_models = len(self.yarn_models)
        
        #some memory
        self.source_mass = np.empty((self.nr_models, self.timesteps + 1), float)
        
        #plot the result every few seconds so outcome becomes visible during calculations
        self.plotevery = self.cfg.get("plot.plotevery")

    def init_yarn(self):
        self.yarn_mass = [0] * len(self.yarn_models)
        self.tstep = 0
        for ind, model in enumerate(self.yarn_models):
            model.do_yarn_init()
            if model.bound_type != 0 : 
                print ' ***********************************************'
                print ' ******  WARNING: Boundary condition not diffusion flux,'\
                      '\n        so yarn does not consider the fabric !!'
                print ' ***********************************************'
            self.yarn_mass[ind] = model.calc_mass(model.init_conc)
            # no mass released at start time
            self.source_mass[ind, self.tstep] = 0

    def __loop_over_yarn(self, nryarns, dx):
        x0 = self.x0
        x02 = np.power(x0[:], 2)
        factor = 4 * self.diff_coef
        n = 0
        termV = np.zeros(len(x0),float)
        expn1V = np.empty(len(x0),float)
        expn2V = np.empty(len(x0),float)
        sol = 0 
        #determine max timestep to take into account in the summation
        # After 100 sec we assume no longer influence
        maxtimestep = int(100 / self.delta_t)
        while n <= nryarns: 
            RV = x02 + n*n*dx*dx
            #determine first timestep to take into account
            firsttimestep = int(np.max(RV*1e-6/(factor*self.delta_t)))
            firsttimestep = max(firsttimestep, self.tstep-maxtimestep)
            for ttt in np.arange(int(firsttimestep), self.tstep):
                tt = ttt+1
                for ttype in np.arange(self.nr_models):
                    # TODO: this should be dx/dy per yarn ttype!!
                    #integralV = np.empty(len(x0),float)
                    #integralH = np.empty(len(x0),float)
                    k=1
                    #for ind, xstart in enumerate(self.x0):
                    #print 'x0', x0, 'delta_t', self.delta_t, 'tt', tt
                    while k <= tt:
                        z2V = -RV/(factor*(self.times[tt]-(k-1)*self.delta_t))
                        for i, z2 in enumerate(z2V):  
                            expn2V[i] = sp.special.expi(z2)
                            if k == tt:
                                expn1V[i] = 0.
                            else:
                                expn1V[i] = sp.special.expi(\
                                    z2 * (self.times[tt]-(k-1)*self.delta_t)\
                                       /(self.times[tt]-k*self.delta_t))
                        termV += self.source_mass[ttype,tt] \
                                    * (expn1V-expn2V) / (factor*np.pi)
                        k += 1

                    #print 'termV', termV 
                        #raw_input()
                    solstep = (
                        self.source_mass[ttype,0] /(factor*np.pi*self.times[tt]) 
                        * np.exp(-RV/(factor*self.times[tt])) + termV 
                                )
                    sol += solstep        
            n += 1
        return sol

    def solve_timestep(self, t):
        print "solve timestep", t
        self.tstep += 1
        print 'tstep', self.tstep
        #raw_input()
        # 1. step one, solve the yarn model, calculate the mass coming out of one yarn and calculate 
        # the corresponding concentration by dividing by the volume of a yarn pi Ry^2
        for ttype, model in enumerate(self.yarn_models):
            rt, rety = model.do_yarn_step(t)
            #print 'rety',rety
            #raw_input()
            tmp = model.calc_mass(rety)
            #V = np.pi * np.power(model.end_point, 2)    
            self.source_mass[ttype, self.tstep] = self.yarn_mass[ttype] - tmp
            #self.source_mass[ttype, self.tstep] /= V
            print 'mass yarn now', tmp, 'prev', self.yarn_mass[ttype], 'release', self.source_mass[ttype, self.tstep]
            self.yarn_mass[ttype] = tmp
            if self.source_mass[ttype, self.tstep] < 0.:
                if abs(self.source_mass[ttype, self.tstep]) < 1e-7:
                    self.source_mass[ttype, self.tstep] = 0.
                    print 'WARNING: small negative release, set to 0'
                else:
                    raise NotImplementedError, 'source must be positive, negative not supported'
        ##raw_input('Continue press ENTER')

        # 2. step two, solve the bednet model
        #    to obtain value near yarn, and at observer
        #    we know that self.source_mass[ttype] has been released since
        #    last step
        
        # concentration is a consequence of all previous releases, so sum 
        # over all times, and compute contribution of that moment.
        self.sol[self.tstep, :] = 0.
        solV = self.__loop_over_yarn(self.nvertyarns, self.dx)
        print 'solV', solV
        solH = self.__loop_over_yarn(self.nhoryarns, self.dy)
        print 'solH', solH        
        self.sol[self.tstep, :] = solV + solH 
        print 'solution', self.sol[self.tstep,:]
        # 3. for next timestep, we need to set correct boundary condition
        #    on the yarn level
        for ind, model in enumerate(self.yarn_models):
            #the source mass is what was present in the overlap
            massoverlapold = self.source_mass[ttype, self.tstep]
            #the new mass there we approximate
            massoverlapnew = self.sol[self.tstep, 0] * self.overlaparea
            massremoved = massoverlapold - massoverlapnew
            print 'prev mass overlap', massoverlapold, 'new', massoverlapnew, 'removed:', massremoved
            #based on removed, we set a source term in the overlap zone of 
            # of the yarn
            model.source_overlap = -massremoved / self.tstep / self.overlaparea
            #print 'prev boundary conc', model.boundary_conc_out, self.sol[self.tstep-1, 0]
            #model.boundary_conc_out = self.sol[self.tstep, 0]
            #print 'boundary conc yarn', model.boundary_conc_out

        fipy.dump.write({
                        'time':self.tstep,
                        'concentration': self.sol[self.tstep,0] },
                        filename=utils.OUTPUTDIR + os.sep + 'bednet_sol_%08d.gz'%(self.tstep)   ,
                        extension='.gz')

    def view_sol(self, times, sol):
        #maxv = np.max(self.sol)
        #minv = np.min(self.sol)
        #print 'max', maxv, minv
        #self.plottimes = np.arange(self.times[0],self.times[-1]+1,self.plotevery)
        plt.ion()
        for ind, pos in enumerate(self.x0[1:]):
            plt.rc("font", family="serif")
            plt.rc("font", size=10)
            width = 4.5  #width in inches
            height = 1.4 #height in inches
            plt.rc("figure.subplot", left=(50/72.27)/width)
            plt.rc("figure.subplot", right=(width-10/72.27)/width)
            plt.rc("figure.subplot", bottom=(14/72.27)/height)
            plt.rc("figure.subplot", top=(height-7/72.27)/height)
            plt.figure(ind)
            plt.gca().set_xlabel('Time [s]')
            plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
            #plt.gca().yaxis.set_major_formatter(pylab.FormatStrFormatter('%e'))
            plt.title('Concentration at position %g mm' % pos)
            plt.plot(self.times, self.sol[:, ind+1])
            #plt.ylim(0, maxv*1.1)
            plt.plot(self.times, np.ones(len(self.times)) * self.treshold, 'b--')
            plt.show()
            plt.savefig(utils.OUTPUTDIR + os.sep 
                        + 'AIconc_%03.1f_mm' % pos + const.FIGFILEEXT)
                
        #fipy.dump.write({plt.plot},filename=utils.OUTPUTDIR + os.sep + 'bednetconc%08.4f.png' % t)

    def init_bednet(self):
        self.sol = np.empty((self.timesteps+1, len(self.x0)), float)
        self.sol[0, :] = 0
        self.init_yarn()

    def run(self, wait=False):
        self.init_bednet()
        for t in self.times[1:]:
            self.solve_timestep(t)
          
        self.view_sol(self.times,self.sol)
        
        #self.view_sol(self.times,self.sol)    

        if wait:
            raw_input("Finished bednet run")
