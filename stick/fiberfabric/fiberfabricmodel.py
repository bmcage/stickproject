#
# Copyright (C) 2012  B. Malengier
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
Module holding a generic diffusion model for a fiberfabric for textile exp. 
"""
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
from __future__ import division
import os.path
import shutil
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import stick.const as const
import stick.fiberfabric.config as conf
import stick.lib.utils.utils as utils

#-------------------------------------------------------------------------
#
#Fipy Imports
#-------------------------------------------------------------------------
from fipy import *

#-------------------------------------------------------------------------
#
# Functions
#
#-------------------------------------------------------------------------

rho_w =  10**(-3) # density water g/mm^3 
M_Molw = 18.01528 # molar mass water g/mol
Kw = 0.6 / 10**3  # heat conductivity water in W / (mm K)
Ka = 0.025 / 10**3# heat conductivity air in W / (mm K) at 298 K
Cvw = 4.1796 / 10**3 # volumetric heat capacity water in J/ (mm^3 K) at 298 K

def porosity(por_no_water, water_content_surface, outside):
    """
    Current porisity based on the porosity with no water, and the water content
    on the surface of the fibers
    """
    por = por_no_water - water_content_surface*(1-por_no_water)
    por[outside] = 1.
    return por

def eff_heat_cond(por_no_water, heatcond_gas, heatcond_fibs, density_fibers, 
        volfracfibers, water_content_surface,
        water_contents_fiber_relative_fiber_weight, outside):
        """ The effective heat conductivity of the fiberfabric, excluding the
            PCM
        """
        Wtilde = water_content_surface
        Wfibers = water_contents_fiber_relative_fiber_weight
        Kg = heatcond_gas # should be 0.025 W/ (m K) for 298 K (=25 C)
        Kfs = heatcond_fibs
        rho_fs = density_fibers
        # first the solid
        Ks = Wtilde*Kw
        Ksnom = Wtilde
        for Kf, Wfiber, rho_f, volfrac in zip(Kfs, Wfibers, rho_fs, volfracfibers):
            Ks += Kf * volfrac / (1-por_no_water) + rho_f/rho_w*Wfiber*Kf
            Ksnom += volfrac / (1-por_no_water) + rho_f/rho_w*Wfiber
        Ks = Ks / Ksnom
        por = porosity(por_no_water, water_content_surface, outside)
        K = por * Kg + (1-por) * Ks
        K[outside] = Kg
        return K

def eff_heat_capacity(por_no_water, heatcap_gas, heatcap_fibs, density_fibers, 
        volfracfibers, water_content_surface,
        water_contents_fiber_relative_fiber_weight, outside):
        """ The effective heat capacity of the fiberfabric, excluding the PCM
        """
        Wtilde = water_content_surface
        Wfibers = water_contents_fiber_relative_fiber_weight
        Cvg = heatcap_gas # should be 0.025 W/ (m K) for 298 K (=25 C)
        Cvfs = heatcap_fibs
        rho_fs = density_fibers
        # first the solid
        Cvs = Wtilde*Cvw
        Cvsnom = Wtilde
        for Cvf, Wfiber, rho_f, volfrac in zip(Cvfs, Wfibers, rho_fs, volfracfibers):
            Cvs += Cvf * volfrac / (1-por_no_water) + rho_f/rho_w*Wfiber*Cvf
            Cvsnom += volfrac / (1-por_no_water) + rho_f/rho_w*Wfiber
        Cvs = Cvs / Cvsnom
        por = porosity(por_no_water, water_content_surface, outside)
        Cv = por * Cvg + (1-por) *Cvs
        Cv[outside] = Cvg
        return Cv

#-------------------------------------------------------------------------
#
# DiffusionModel class 
#
#-------------------------------------------------------------------------
class FiberFabricModel(object):
    """
    FiberFabricModel is a special diffusion model for a fabric modeled as a
    collection of fibers and other materials (PCM, microcap, ...),
    in which we will do textile experiments.
    The Fabric is solved with fipy, it connects to the outside room via 
    direct domain splitting. 
    The fibers and PCM are incorporated via volume averaging upscaling
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        """
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')
        #
        self.modelcomp = self.cfg.get("component.present")
        #time data
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.delta_t = self.times[1] - self.times[0]
        if self.verbose:
            print "Timestep used in fiberfabric model:", self.delta_t
        
        #construct cfg for fiber
        self.cfg_fiber = []
        for filename in self.cfg.get('fiber.fiber_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            from stick.fiber.config import FiberConfigManager
            self.cfg_fiber += [FiberConfigManager.get_instance(filename)]
            #set values from the fiberfabric on this inifile
            self.cfg_fiber[-1].set("time.time_period", self.time_period)
            if self.cfg_fiber[-1].get("time.dt") > self.cfg.get("time.dt"):
                self.cfg_fiber[-1].set("time.dt", self.cfg.get("time.dt"))

        #construct cfg for pcm
        self.cfg_pcm = []
        for filename in self.cfg.get('pcm.pcm_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                            os.path.dirname(self.cfg.filename), filename))
            from stick.pcm.config import PCMConfigManager
            self.cfg_pcm += [PCMConfigManager.get_instance(filename)]
            #set values from the fiberfabric on this inifile
            self.cfg_pcm[-1].set("time.time_period", self.time_period)
            self.cfg_pcm[-1].set("time.dt", self.cfg.get("time.dt"))

        self.plotevery = self.cfg.get("plot.plotevery")
        self.writeevery = self.cfg.get("plot.writeevery")
        self.writeoutcount = 0
        
        #allow a multiscale model to work with a source in overlap zone
        self.source_overlap = 0.
        
        self.initialized = False

    def create_mesh(self):
        """
        Create a mesh for use in the model
        """
        self.length = self.cfg.get('fabric.length')
        self.width = self.cfg.get('fabric.width')
        self.height = self.cfg.get('fabric.height')
        el_length = self.cfg.get('discretization.el_length')
        el_width = self.cfg.get('discretization.el_width')
        el_height = self.cfg.get('discretization.el_height')
        
        dxe = self.length / el_length # mesh size in x direction 
        dye = self.width / el_width # mesh size in x direction
        dze = self.height / el_height # mesh size in x direction
        outsidesize = self.height
        dz = [dze] * el_height
        dz = dz + dz
        dz = np.array(dz)
        
        dx = [dxe] * el_length
        dx = np.array(dx)
        dy = [dye] * el_width
        dy = np.array(dy)

        # construct the fipy mesh
        self.mesh = Grid3D(dx=dx, dy=dy, dz=dz)
        
        xc, yc, zc = self.mesh.cellCenters
        self.outsidecells = (zc > self.height)
        xfc, yfc, zfc = self.mesh.faceCenters
        self.textilesurface = (self.height - dze/10 * self.height < zfc ) & \
                            ( zfc < self.height + dze/10 * self.height)
        self.toptextlayer = (self.height - 2*dze/3 < zc) & \
                            (zc < self.height )
        self.textilebody = (zc < self.height )
        
        # we construct fibermodels and pcmmodels in every fipy cell
        self.fiber_model = [0] * len(self.cfg_fiber)
        for ind in range(len(self.fiber_model)):
            self.fiber_model[ind] = []
        self.pcm_model = [0] * len(self.cfg_pcm)
        for ind in range(len(self.pcm_model)):
            self.pcm_model[ind] = []
        

    def initial_fabric(self):
        """
        Do initial setup to solve the problem
        """
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        
        self.init_concvap = eval(self.cfg.get('initial.init_concvap'))
        self.init_concair = eval(self.cfg.get('initial.init_concair'))
        self.init_temp = eval(self.cfg.get('initial.init_temp'))
        
        cellCenter_x, cellCenter_y, cellCenter_z = self.mesh.cellCenters
        initialConcVap = []
        initialConcAir = []
        initialTemp = []
        for i_x, i_y, i_z in zip(cellCenter_x, cellCenter_y, cellCenter_z):
            initialConcVap.append(self.init_concvap(i_x, i_y, i_z))
            initialConcAir.append(self.init_concair(i_x, i_y, i_z))
            initialTemp.append(self.init_temp(i_x, i_y, i_z))
        initialConcVap = sp.array(initialConcVap)
        initialConcAir = sp.array(initialConcAir)
        initialTemp = sp.array(initialTemp)

        from stick.fiber1d.fibermodel import FiberModel
        from stick.pcm.pcmmodel import PCMModel

        self.nrcells = len(self.mesh.cellCenters[0])
        self.nrcellsfabric = len(self.mesh.cellCenters[0][self.textilebody])
        self.nr_pcmmodels = len(self.cfg_pcm)
        self.nr_fibermodels = len(self.cfg_fiber)
        if self.nr_pcmmodels:
            self.pcm_E = np.empty((self.nrcellsfabric, self.nr_pcmmodels), float)
            self.tmp_E = np.empty((self.nrcellsfabric, self.nr_pcmmodels), float)
            self.source_energy = np.empty((self.nrcells, self.nr_pcmmodels), float)
        if self.nr_fibermodels and self.modelcomp:
            self.fiber_mass = np.empty((self.nrcellsfabric, self.nr_pcmmodels), float)
        
        for x, y, z in zip(self.mesh.cellCenters[0][self.textilebody], 
                self.mesh.cellCenters[1][self.textilebody],
                self.mesh.cellCenters[2][self.textilebody]):
            if self.modelcomp:
                for ind, cfg in enumerate(self.cfg_fiber):
                    #create fiber model
                    self.fiber_model[ind] += [FiberModel(cfg)]
            for ind, cfg in enumerate(self.cfg_pcm):
                #set correct init temp of the PCM
                cfg.set("init.init_temp", 'lambda x: %g' % 
                                self.init_temp(x, y, z))
                #create fiber model
                self.pcm_model[ind] += [PCMModel(cfg)]

        self.step_old_sol_vap = initialConcVap
        self.step_old_sol_air = initialConcAir
        self.step_old_sol_temp = initialTemp
        self.step_old_sol_temp_body = self.step_old_sol_temp[self.textilebody]

        self.conVap = CellVariable(name = "Conc. water vapour", 
                                    mesh = self.mesh,
                                    value = initialConcVap)
        self.concAir = CellVariable(name = "Conc. air", 
                                    mesh = self.mesh,
                                    value = initialConcAir)
        self.Temp = CellVariable(name = "Temperature", 
                                    mesh = self.mesh,
                                    value = initialTemp)

        self.heat_cond = CellVariable(name='Effective Heat Conductivity',
                                      mesh=self.mesh)
        self.heat_cap = CellVariable(name="Effective Heat Capacity",
                                     mesh=self.mesh)
        
        # constans in the equations
        self.porosity = self.cfg.get("fabric.porosity")
        self.volfrac = self.cfg.get("fiber.volfrac")
        # test
        totfracs = self.porosity
        for frac in self.volfrac:
            totfracs += frac
        assert totfracs == 1., "fraction of voids and fibers does not sum to 1"
        self.Da = self.cfg.get('fabriccoeff.diff_coef')
        self.Ka = self.cfg.get('fabriccoeff.therm_cond_K') * (10**(-3))
        self.ca = self.cfg.get('fabriccoeff.spec_heat_c')
        self.Kf = []
        self.cf = []
        self.rhof = []
        self.Wfinit = [] #original water content absorbed relative to fiber density
        for cfg in self.cfg_fiber:
            self.Kf += [cfg.get("fiber.therm_cond_K") * 10**(-3)] #in  W / (mm K)
            self.cf += [cfg.get("fiber.spec_heat_c")]
            self.rhof += [cfg.get("fiber.density") * 10**(-3)] # in  g/mm^3
            Wfinit = np.empty(self.nrcells, float)
            Wfinit[:] = cfg.get("fiber.water_absorbed_rel_dens")
            Wfinit[self.outsidecells] = 0.
            self.Wfinit += [Wfinit]
        #we store the amount of water around the fibers (percentage of voids)
        self.W_v = np.empty(self.nrcells, float)
        self.W_v[:] = self.cfg.get("fabric.water_content_voids")
        self.W_v[self.outsidecells] = 0.
        #derive water content surface as in articles
        self.W_s = self.W_v * self.porosity / (1-self.porosity)

        self.BC_T_type = self.cfg.get('boundary.T_type')
        self.BC_T_dir = None
        if self.BC_T_type == 'heatingplate':
            self.BC_T_dir = self.cfg.get('boundary.T_dir')
        elif self.BC_T_type == 'insulated':
            pass
        else:
            print "ERROR: unknown Temperature boundary type"
            sys.exit()
        self.with_overlap = self.cfg.get("boundary.overlapzone")
        self.outTemp = self.cfg.get("boundary.outtemp")
        self.Temp[self.outsidecells] = self.outTemp
        self.minTempup = np.empty(len(self.times), float)
        self.maxTempup = np.empty(len(self.times), float)
        self.avgTempup = np.empty(len(self.times), float)
        self.minTempup[0] = np.min(self.Temp[self.toptextlayer].value)
        self.maxTempup[0] = np.max(self.Temp[self.toptextlayer].value)
        self.avgTempup[0] = (self.minTempup[0] + self.maxTempup[0]) / 2
        
        self.viewer = None
        self.viewer = Viewer(vars = self.Temp, title = 'Temperature Distribution', 
                            datamin = 15., datamax = 45.)
        self.viewer.plot()
        #raw_input("take the example of the initial condition")
        self.viewerplotcount = 1
        if self.plotevery:
            self.viewerplotcount = self.viewerplotcount % self.plotevery

        #now we initialize our submodels
        self.solve_fiber_init()
        self.solve_pcm_init()
        #initialize upscale data
        self.volfrac_pcm = self.cfg.get("pcm.volfrac")
        self.vol_pcm = []
        self.nr_pcm_permm3 = []
        for cfg in self.cfg_pcm:
            self.vol_pcm += [4/3*np.pi * cfg.get("pcm.radius")]
        for vol, volfrac in zip(self.vol_pcm, self.volfrac_pcm):
            self.nr_pcm_permm3 += [volfrac/vol]
        
        #now print out some data for user to verify
        print 'sample size in mm is', self.length, 'x', self.width, 'x', self.height
        self.volume = self.length * self.width * self.height
        print '  ==> volume', self.volume, 'mm3'
        self.weight = 0.
        for cfg, volfracf in zip(self.cfg_fiber, self.volfrac):
            self.weight += self.volume * volfracf * cfg.get("fiber.density") * 1e-3
        print '  ==> weight', self.weight, 'g/mm3'
        print 'Last due to fiber volume fractions of', self.volfrac
        if self.volfrac_pcm:
            print '\nPCM'
            print 'volume fractions of', self.volfrac_pcm
            print 'weight added in g/mm3',
            weightpcms = 0.
            for volfrac, nrpcm, cfg in zip(self.volfrac_pcm, self.nr_pcm_permm3, 
                                            self.cfg_pcm):
                weight = volfrac * self.volume * cfg.get("pcm.density") * 1e-6
                print weight,
                weightpcms += weight
            print ' '
            print 'weight pickup pcms is', weightpcms/self.weight * 100, '%'
            print 'total weight fabric + pcm', self.weight + weightpcms
        #raw_input('\nPress key to start computation\n')

    def solve_fiber_init(self):
        """
        Initialize the solvers that do the fiber simulations
        """
        for type, models in enumerate(self.fiber_model):
            #first index is type of fiber
            for ind, model in enumerate(models):
                #ind is position of the cell in the datastructure
                model.run_init()
                model.solve_init()
                print 'ERROR: not implemented yet to run fibers!'
                sys.exit()
                #rebind the out_conc method to a call to fiberfabric
                model.set_userdata(self.get_data(ind))
                model.out_conc = lambda t, data: self.out_conc(data, t)
                self.fiber_mass[ind, type] = model.calc_mass(model.initial_c1)

    def solve_pcm_init(self):
        """
        Initialize the solvers that do the pcm simulations
        """
        for type, models in enumerate(self.pcm_model):
            #first index is type of pcm
            for ind, model in enumerate(models):
                #ind is position of the cell in the datastructure
                model.run_init()
                model.solve_init()
                #rebind the Tout method to a call to fiberfabric
                model.set_userdata(self.get_data(ind))
                model.Tout = lambda t, data: self.out_temp(data, t)
                self.pcm_E[ind, type] = model.calc_energy(model.initial_T_in, 
                            model.initial_T_out)

    def get_data(self, cellnr):
        return cellnr

    def out_temp(self, data, t):
        """
        return the temperature at cellnr data at time t
        """
        timenowyarn = self.step_old_time
        if t >= timenowyarn:
            #return data
            return self.step_old_sol_temp_body[data]
        raise ValueError, 'out temperature should only be requested at a later time'

    def update_heat_cond(self):
        """
        Update heat conductivity in the fabric based on the previous solution
        """
        ## we assume for now no water change
        heat_cond = eff_heat_cond(self.porosity, self.Ka, self.Kf,
                    self.rhof, self.volfrac, self.W_s, self.Wfinit,
                    self.outsidecells)
        self.heat_cond.value = heat_cond

    def update_heat_cap(self):
        """
        Update heat capacity in the fabric based on the previous solution
        """
        ## we assume for now no water change
        heat_cap = eff_heat_capacity(self.porosity, self.ca, self.cf,
                    self.rhof, self.volfrac, self.W_s, self.Wfinit,
                    self.outsidecells)
        self.heat_cap.value = heat_cap


    def do_pcm_step(self, stopt):
        """
        Solve the diffusion process on the fiber up to stoptime, starting
        from where we where last. 
        The flux is the BC: S*h(C_equi - C_yarn(t))*H(C-C_b,C_equi-C_yarn(t))
        """
        #Outside temp is current temp at that position via func out_temp!!
        for type, models in enumerate(self.pcm_model):
            #first index is type of pcm
            for ind, model in enumerate(models):
                model.solve_step()
                self.tmp_E[ind, type] = model.calc_last_energy()
        self.source_energy[self.textilebody, :] = self.pcm_E[:, :] - self.tmp_E[:,:]
        self.pcm_E[:, :] = self.tmp_E[:,:]

    def solve_fabric(self):
        """
        Solve the unknowns in the fiberfabric
        Solve fabric in steps of delta t. This does:
           1. solve the pcm for delta t
           2. set correct source term for the fabric
           3. solve the fabric for delta t
        """
        #input the transient equation
        self.update_heat_cap()
        self.update_heat_cond()
        print 'Min - Max heat cond', min(self.heat_cond.value), max(self.heat_cond.value)
        print 'Min - Max heat cap ', min(self.heat_cap.value), max(self.heat_cap.value)
        
        rhs = DiffusionTerm(coeff=self.heat_cond)
        self.source_var = []
        for pcmtype, Np in enumerate(self.nr_pcm_permm3):
            self.source_var += [CellVariable(name="Source pcm %d" % pcmtype,
                                     mesh=self.mesh)]
            rhs = rhs + self.source_var[-1]

        self.eqTmp = TransientTerm(coeff=self.heat_cap) == rhs
                            

        #boundary conditions
        if self.BC_T_type == 'heatingplate':
            self.Temp.constrain(self.BC_T_dir, where=self.mesh.facesFront)
        elif self.BC_T_type == 'insulated':
            pass
        else:
            pass
        if not self.with_overlap:
            #we set on top boundary the outside temperature
            self.Temp.constrain(self.outTemp, where=self.mesh.facesBack)

        #all other boundaries are automatically Neumann BC
        
        #now loop in time to solve
        t = self.initial_t
        stop_time = self.times[-1]
        compute = True
        post = 0
        while compute:
            t += self.delta_t
            post += 1
            print 'computing time', t
            if t >= stop_time -self.delta_t / 100:
                t = stop_time
                compute = False

            # step 1 solve the pcm models. 
            self.do_pcm_step(t)
            
            # step 2 update source by setting value of the source variables
            for type, var in enumerate(self.source_var):
                # upscale the energy of a single PCM
                var.value = self.source_energy[:, type] /self.delta_t \
                            * self.nr_pcm_permm3[type]
            
            # step 3 solve fabric model
            self.eqTmp.solve(var = self.Temp,
                          dt = self.delta_t)
            
            #store solution
            self.step_old_time += self.delta_t
            self.step_old_sol_temp[:] = self.Temp.value[:]
            self.step_old_sol_temp_body = self.step_old_sol_temp[self.textilebody]
            
            #store data for plot
            self.minTempup[post] = np.min(self.Temp[self.toptextlayer].value)
            self.maxTempup[post] = np.max(self.Temp[self.toptextlayer].value)
            self.avgTempup[post] = (self.minTempup[post] + self.maxTempup[post]) / 2
        
            if self.viewer is not None and self.viewerplotcount == 0:
                    self.viewer.plot()
                    outvtk =  utils.OUTPUTDIR + os.sep + \
                                        'fiberfabricTemp%08.3f.vtk' % t
                    invtk = self.viewer.vtkcellfname
                    shutil.copy(invtk, outvtk)
            self.viewerplotcount += 1
            self.viewerplotcount = self.viewerplotcount % self.plotevery

            if self.writeoutcount == 0:
                pass
            self.writeoutcount += 1
            self.writeoutcount = self.writeoutcount % self.writeevery

        #now we plot the min and max temp in upper layer
        plt.figure(num=None)
        plt.plot(self.times, self.minTempup, 'b', 
                    self.times, self.maxTempup, 'r',
                    self.times, self.avgTempup, 'k--')
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (C)')
        plt.show()
        
        raw_input("Finished <press return>.....")

    def run(self):
        """
        Method that is called to do a full model run
        """
        self.create_mesh()
        self.initial_fabric()
        self.solve_fabric()
