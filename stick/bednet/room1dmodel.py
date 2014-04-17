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
    Here a room of LxWxH is modeled with a net spanning the room at L/2
    This can be reduced to a 1D model.
    The yarn in the net is modeled as 1D cilindrical over 2 x Ry. The mass 
    outside the yarn (Ry, 2Ry) (the overlap zone) is put in the room model
    which runs from Ry to L/2. Note that this is a cube domain!
    The room allows for some ventilation L/s which is modeled as removing
    the corresponding volume of AI from everywhere
"""
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
from __future__ import division, print_function
import os.path
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import pi

MAX_STORE_LENGTH = 1000
MAX_PLOT_LENGTH = 500000
INSPECT_MEM = False

HAVE_ODES = False
try:
    from scikits.odes import ode as sc_ode
    HAVE_ODES = True
except:
    print ('Could not load scikits.odes, odes solver not available')

#consider absorption problem possible or not?
ABSORPTION = True

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import stick.const as const
import stick.lib.utils.utils as utils
from stick.yarn.config import YarnConfigManager
from stick.yarn1d.yarn1dmodel import Yarn1DModel

#-------------------------------------------------------------------------
#
# DiffusionModel-fabric1d class 
#
#-------------------------------------------------------------------------
class Room1DModel(object):
    """
    upscaling to fabric1d domain
    Here a room of LxWxH is modeled with a net spanning the room at L/2
    This can be reduced to a 1D model.
    The yarn in the net is modeled as 1D cilindrical over 2 x Ry. The mass 
    outside the yarn (Ry, 2Ry) (the overlap zone) is put in the room model
    which runs from Ry to L/2. Note that this is a cube domain!
    The room allows for some ventilation L/s which is modeled as removing
    the corresponding volume of AI from everywhere
    """
    def __init__(self, config):
        #indicate if only one side of bednet is free, or both are free
        self.singleside = False
        #other settings from config
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')
        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.timesteps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        #set correct delta_t
        self.delta_t = self.time_period / self.timesteps
        if self.verbose:
            print ("Timestep used in bednet model:", self.delta_t)
        self.initconc = eval(self.cfg.get('initial.init_conc'))

        self.dx = self.cfg.get('domain.dx')
        self.dy  = self.cfg.get('domain.dy')
        #size room in mm
        self.room_L = self.cfg.get("domain.room_L") * 1000
        self.room_W = self.cfg.get("domain.room_W") * 1000
        self.room_H = self.cfg.get("domain.room_H") * 1000
        self.nvertyarns = self.room_W/self.dx
        self.nhoryarns = self.room_H/self.dy
        self.diff_coef = self.cfg.get('diffusion.diff_coef')
        self.saturation_conc = self.cfg.get('active_component.saturation_conc')
        self.treshold = self.cfg.get('active_component.treshold_effect')
        self.x0 = self.cfg.get('observer.x0')
        #define whether there is the ventilation existing
        self.ventilation = self.cfg.get('domain.ventilation')
        self.vel_ventilation = self.cfg.get('domain.vel_ventilation')
        #we set a distance for the yarn bc
        EXTFRAC = 1.
        self.cfg_yarn = []
        self.radius_yarn = []
        for filename in self.cfg.get('sample.yarn_config'):
            if not os.path.isabs(filename):
                filename = os.path.normpath(os.path.join(
                        os.path.dirname(self.cfg.filename), filename))
            self.cfg_yarn.append(YarnConfigManager.get_instance(filename))
            #set values from the yarn on this inifile
            print ('time', self.time_period)
            self.cfg_yarn[-1].set("time.time_period", self.time_period)
            self.cfg_yarn[-1].set("boundary.dist_conc_out", float(self.x0[0]))
            self.cfg_yarn[-1].set("boundary.D_out", self.diff_coef)
            self.cfg_yarn[-1].set("boundary.conc_out", 
                    float(self.initconc(self.cfg_yarn[-1].get("domain.yarnradius"))))
            self.cfg_yarn[-1].set("domain.useextension", True)
            ## How much overlap region? Take one yarn radius for now
            self.cfg_yarn[-1].set("domain.extensionfraction", EXTFRAC)
            self.radius_yarn.append(self.cfg_yarn[-1].get("domain.yarnradius"))
            assert self.radius_yarn[-1] == self.radius_yarn[0], 'ERROR, yarns'\
                ' must have equal radius for now, as massperyarn is equally '\
                'distributed over the yarns'

        #we want the overlap zone of the yarns to end at 2* maximum radius:
        self.maxyarnrad = max(self.radius_yarn)
        self.minyarnrad = min(self.radius_yarn)
        voloverlapyarn = (np.pi*((self.maxyarnrad * (1 + EXTFRAC))**2 -
                                (self.maxyarnrad**2)) * (self.nhoryarns 
                                    * self.room_W + self.nvertyarns 
                                    * self.room_H)
                         )
        #self.endoverlap = self.maxyarnrad * (1 + EXTFRAC)
        if self.singleside:
            self.endoverlap = self.minyarnrad + voloverlapyarn /  self.room_W /  self.room_H
        else:
            self.endoverlap = self.minyarnrad + voloverlapyarn /  self.room_W /  self.room_H / 2
        for config, rad in zip(self.cfg_yarn, self.radius_yarn):
            config.set("domain.extensionfraction", EXTFRAC)
        
        #create yarn models
        self.yarn_models = []
        for cfg in self.cfg_yarn:
            self.yarn_models.append(Yarn1DModel(cfg))
        self.nr_models = len(self.yarn_models)
        #some memory
        self.source_mass = np.empty(self.nr_models, float)
        #self.mass_build =  [0.,]
        
        #plot the result every few seconds so outcome becomes visible during calculations
        self.plotevery = self.cfg.get("plot.plotevery")
        self.viewerwritecount = 0
        self.writeevery = self.cfg.get("plot.writeevery")
        
        #now some output on density
        self.volbednet = 0.
        self.surfbednet = self.room_H * self.room_W
        for rad in self.radius_yarn:
            print ('vert vol yarns', self.nvertyarns * pi * rad**2 * self.room_H, 'mm3')
            print ('horz vol yarns', self.nhoryarns * pi * rad**2 * self.room_W, 'mm3')
            self.volbednet += self.nvertyarns * pi * rad**2 * self.room_H
            self.volbednet += self.nhoryarns * pi * rad**2 * self.room_W
        print ('volume_bednet space =', (2 * self.maxyarnrad * self.room_H
                                                * self.room_W))
        # The total volume of the bednet incl void space is the area of the
        # net * the tickness of the net.
        # This thickness is taken to be twice a yarndiameter.
        self.totalvolume_net = self.room_H * self.room_W * 4 * self.maxyarnrad
        self.voidvolume = self.totalvolume_net - self.volbednet
        self.densitybednet =  self.volbednet / (2 * self.maxyarnrad * self.room_H
                                                * self.room_W)
        self.fabporosity = self.voidvolume / self.totalvolume_net
        
        self.initialized = False
    
        self.yarnconc_center = np.empty((self.timesteps, 2),float)
        self.yarnconc_surface = np.empty((self.timesteps, 2),float)

    def times(self, timestep, end=None):
        """ Compute the time at one of our steps
        If end is given, all times between step timestep and step end are 
        returned as a list, with end included
        """
        if end is None:
            return timestep * self.delta_t
        else:
            begin = timestep * self.delta_t
            end = end * self.delta_t
            return np.linspace(begin, end, end-begin + 1)

    def create_mesh(self):
        """ create a mesh for the room model
        """
        self.begin_point = self.minyarnrad
        self.end_point = self.room_L / 2
        self.nr_edge = self.cfg.get('domain.n_edge')
        self.nr_cell = self.nr_edge - 1
        self.grid_edge = np.empty(self.nr_edge, float)
        self.grid_edge[0] = self.begin_point
        self.grid_edge[1:] = np.linspace(self.endoverlap, self.end_point,
                                self.nr_edge-1)
        if self.singleside:
            self.overlapvolume = (self.endoverlap - self.begin_point) \
                                * self.room_H * self.room_W
        else:
            self.overlapvolume = (self.endoverlap - self.begin_point) \
                                * self.room_H * self.room_W * 2
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_x = self.grid_edge[1:] - self.grid_edge[:-1]
        
        self.plotdata = []
        for xplot in self.x0:
            assert self.grid[0] < xplot < self.grid[-1], "%f < %f < %f "\
                "Not satisfied, observer out of domain" % (self.grid[0], 
                                                    xplot, self.grid[-1])
            for ind, xcell in enumerate(self.grid):
                if xcell >= xplot:
                    interpol_start = (xcell-xplot)/(self.grid[ind]-self.grid[ind-1])
                    self.plotdata.append((ind-1, interpol_start))
                    break

    def upscale_yarnmass(self, mass):
        """
        Upscale the mass in one yarn, to the mass in the entire bednet,
        returns the mass
        """
        return (self.nhoryarns * self.room_W + self.nvertyarns * self.room_H) * mass

    def calc_mass(self):
        """
        Calculate the mass in the room and the bednet of Active Component
        at this specific state the model is in.
        """
        #First, the mass in the yarns
        yarnmass = [None] * len(self.yarn_models)
        yarnmassoverlap = [None] * len(self.yarn_models)
        for ttype, model in enumerate(self.yarn_models):
            yarnmass[ttype] = model.calc_mass(model.step_old_sol) 
            yarnmassoverlap[ttype] = model.calc_mass_overlap(model.step_old_sol)
        #Next, the upscaled mass in the yarns
        totyarnmass = [None] * len(self.yarn_models)
        totyarnmassoverlap = [None] * len(self.yarn_models)
        for ttype, (massy, massyo) in enumerate(zip(yarnmass, yarnmassoverlap)):
            totyarnmass[ttype] = self.upscale_yarnmass(massy)
            totyarnmassoverlap[ttype] = self.upscale_yarnmass(massyo)

        #Next, the mass in the room, divided in overlapzone and rest.
        if self.singleside:
            roomoverlapmass = self.overlapvolume * self.step_old_sol[0]
            roommass = np.sum(self.delta_x[1:] * self.room_H * self.room_W * self.step_old_sol[1:])
        else:
            # factor 2 because we only model half of the room
            roomoverlapmass = self.overlapvolume * self.step_old_sol[0]
            roommass = 2*np.sum(self.delta_x[1:] * self.room_H * self.room_W * self.step_old_sol[1:])

        return (yarnmass, yarnmassoverlap, totyarnmass, totyarnmassoverlap,
                roommass, roomoverlapmass)

    def initial_room(self):
        """ initial concentration in the room domain
        """
        self.init_conc = np.empty(self.nr_cell, float)
        self.init_conc[:] = self.initconc(self.grid[:])

    def init_yarn(self):
        self.yarn_mass = [0] * len(self.yarn_models)
        self.yarn_mass_overlap = [0] * len(self.yarn_models)
        self.yarn_mass_overlap_old = [0] * len(self.yarn_models)
        #self.fibermass = [0] * len(self.yarn_models)
        self.tstep = 0
        for ind, model in enumerate(self.yarn_models):
            model.do_yarn_init()
            if model.bound_type != 0 : 
                print (' ***********************************************')
                print (' ******  WARNING: Boundary condition not diffusion flux,'\
                      '\n        so yarn does not consider the fabric !!')
                print (' ***********************************************')
            self.yarn_mass[ind] = model.calc_mass(model.init_conc)
            self.yarn_mass_overlap[ind] = model.calc_mass_overlap(model.init_conc)
            # no mass released at start time
            self.source_mass[ind] = 0.
            #self.fibermass = model.get_fiber_mass()
            #print(self.fibermass)

    def f_conc_ode(self, t, conc_x, diff_u_t):
        """
        Solving the room 1D diffusion equation: 
        
          \partial_t (C) =  \partial_x (D \partial_x C) + Source
        
        with Source the concentration amount per time unit added/removed at x. 
        Solution is obtained by integration over a cell, so
        
           \delta x d_t (C) = flux_right - flux_left + Source (\delta x)
        
        so 
        
          d_t C = 1 / (\delta x) * (flux_right - flux_left) + Source 
        
        We have homogeneous Neumann BC
        """
        grid = self.grid
        n_cellcenters = len(grid)
        #Initialize the flux rate on the edges
        flux_edge = self.__tmp_flux_edge
        #set flux on edge 0, self.nr_edge-1
        flux_edge[0] = 0.
        flux_edge[-1] = 0.

        #calculate flux rate in each edge of the domain
        flux_edge[1:self.nr_edge-1] = (2 * self.diff_coef *
            (conc_x[1:]-conc_x[:-1]) / (self.delta_x[:-1]+self.delta_x[1:])
            )
        if self.ventilation == 'advection':
            raise NotImplementedError, 'This piece needs testing before use!' # needs testing before activating this!
            flux_edge[-1] = - self.vel_ventilation * conc_x[-1] * self.delta_t
            flux_edge[1:self.nr_edge-1] += - 2 * self.vel_ventilation \
                            * (conc_x[1:] + conc_x[:-1]) / 2.
        elif self.ventilation == "zero_on_edge":
            #we assume always 0 outside the edge, this means outside is refreshed
            flux_edge[self.nr_edge-1] = (self.diff_coef *
                                        (0-conc_x[-1]) / self.delta_x[-1] )
            ##print ('flux edge room', flux_edge[self.nr_edge-1])

        diff_u_t[:] = ((flux_edge[1:]-flux_edge[:-1])
                            / self.delta_x[:]
                      )
        ## we add a source term in the first cell where the overlap is
        diff_u_t[0] += self.source_room_from_yarn

##    def f_conc_ode_vel(self, t, conc_x, diff_u_t, vel_ventilation):
##        """
##        Solving the room 1D diffusion equation: 
##        
##          \partial_t (C) =  \partial_x (D \partial_x C) - v\partial_x C + Source
##        
##        with Source the concentration amount per time unit added/removed at x. 
##        Solution is obtained by integration over a cell, so
##        
##           \delta x d_t (C) = flux_right - flux_left + Source (\delta x)
##        
##        so 
##        
##          d_t C = 1 / (\delta x) * (flux_right - flux_left) + Source 
##        
##        We have homogeneous Neumann BC
##        """		
##        grid = self.grid
##	n_cellcenter = len(grid)
##        flux_edge = self.__tmp_flux_edge
##        #set flux on edge 0, self.nr_edge-1
##        flux_edge[0] = 0.
##	flux_edge[-1] = -vel_ventilation * conc_x[-1] 
##        #flux_edge[-1] = 0.
##        #calculate flux rate in each edge of the domain
##        flux_edge[1:self.nr_edge-1] = (2 * self.diff_coef *
##            (conc_x[1:]-conc_x[:-1]) / (self.delta_x[:-1]+self.delta_x[1:])
##            ) - vel_ventilation * (conc_x[1:] + conc_x[:-1]) / (self.delta_x[:-1]
##	    +self.delta_x[1:])
##        diff_u_t[:] = ((flux_edge[1:]-flux_edge[:-1])
##                            / self.delta_x[:]
##                      )
##        ## we add a source term in the first cell where the overlap is
##        diff_u_t[0] += self.source_room_from_yarn

    def solve_ode_init(self):
        """
        Initialize the ode solver
        """
        self.initial_t = self.times(0)
        self.step_old_time = self.initial_t
        #storage for solution
        self.solstoreind = 0
        if self.timesteps+1 > MAX_STORE_LENGTH:
            self.solpart = np.empty((MAX_STORE_LENGTH, self.nr_cell), float)
        else:
            self.sol = np.empty((self.timesteps+1, self.nr_cell), float)
            self.solpart = self.sol
        self.solpart[0, :] = self.init_conc[:]
        self.ret_y = np.empty(self.nr_cell, float)
        self.__tmp_flux_edge = np.zeros(self.nr_cell+1, float)
        self.tstep = 0
        self.step_old_sol = self.solpart[0]
        
        self.solver = sc_ode('cvode', self.f_conc_ode,
                             min_step_size=1e-8, rtol=1e-6, atol=1e-6, 
                             max_steps=50000, lband=1, uband=1)
        print (self.step_old_time)
        self.solver.init_step(self.step_old_time, self.init_conc)
        self.initialized = True

    def do_ode_step(self, stoptime):
        """Solve the roommodel up to stoptime, continuing from the present
           state, return the time, concentration after step
        """
        self.solver.init_step(self.step_old_time, self.step_old_sol)
        self.solver.set_options(tstop=stoptime)
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'

        flag, realtime = self.solver.step(stoptime, self.ret_y)
        if flag < 0:
            raise Exception, 'could not find solution, flag %d' % flag
        assert np.allclose(realtime, stoptime), "%f %f" % (realtime, stoptime)
        return stoptime, self.ret_y

    def solve_timestep(self, t):
        print ("solve up to time", t, "s")
        self.tstep += 1
        # 1. step one, solve the yarn model, calculate the mass coming out of one yarn and calculate 
        # the corresponding concentration by dividing by the volume of a yarn pi Ry^2
        for ttype, model in enumerate(self.yarn_models):
            rt, rety = model.do_yarn_step(t)
            self.yarnconc_center[self.tstep-1,0] = t
            self.yarnconc_center[self.tstep-1,1] = rety[0]
            self.yarnconc_surface[self.tstep-1,0] = t
            self.yarnconc_surface[self.tstep-1,1] = rety[-1]
            #filedata= open(utils.OUTPUTDIR + os.sep + "yarnconc_%05d" %t + ".txt",'w')
            #filedata.write("conc on %.10f is %s" % (t,rety))
            #filedata.close()
            #filedata= open(utils.OUTPUTDIR + os.sep + "yarnconc_center%05d" %t + ".txt",'w')
            #filedata.write("conc on %.10f is %s" % (t,rety[0]))
            #filedata.close()
            #filedata= open(utils.OUTPUTDIR + os.sep + "yarnconc_surface_%05d" %t + ".txt",'w')
            #filedata.write("conc on %.10f is %s" % (t,rety[-1]))
            #filedata.close()
            tmp = model.calc_mass(rety)
            tmp_overlap = model.calc_mass_overlap(rety)
            # mass that goes into overlap is the mass that disappeared.
            self.source_mass[ttype] = tmp_overlap \
                - (self.yarn_mass_overlap[ttype] + model.source_overlap
                                            * self.delta_t * model.areaextend)
            self.yarn_mass[ttype] = tmp
            self.yarn_mass_overlap_old[ttype] = self.yarn_mass_overlap[ttype]
            self.yarn_mass_overlap[ttype] = tmp_overlap
            if (ABSORPTION != True):
                # we check on absorption, and give error if too big
                if self.source_mass[ttype] < 0.:
                    print ("source mass", self.source_mass[ttype])
                    if abs(self.source_mass[ttype]) < 100:
                        #self.source_mass[ttype, self.tstep] = 0.
                        print ('WARNING: small negative release, reduce timestep fiber/yarn if needed')
                    else:
                        raise NotImplementedError, 'source must be positive, negative not supported'

        # 2. step two, solve the room model
        #    to obtain new concentration value near yarn.
        #    We know that self.source_mass[ttype] has been released in the 
        #    overlap region since last step
        # 2.a upscale source_mass (mass in ring zone area) to a source per second per mm^3
        # concentration is a consequence of all yarn types, so sum 
        # over yarn types, and compute contribution of that moment.
        concreleased = (self.nhoryarns * self.room_W + self.nvertyarns * self.room_H) \
                        * np.sum(self.source_mass[:]) / self.overlapvolume

        self.source_room_from_yarn = concreleased / self.delta_t
        
        # 2.b solve the room model
        self.step_old_time, self.step_old_sol = self.do_ode_step(t)
        if self.tstep % MAX_STORE_LENGTH == 0 :
            #dump to file, and restart
            self.dump_sol(self.solstoreind)
            self.solstoreind += 1
        self.solpart[self.tstep % MAX_STORE_LENGTH, :] = self.step_old_sol[:]
        # 3. for next timestep, we need to set correct boundary condition
        #    on the yarn level, so downscale the mass to keep mass balance
        ##massdiff = massoverlapnew - massoverlapold
        massperyarn = (self.step_old_sol[0] * self.overlapvolume 
            / (self.nhoryarns * self.room_W + self.nvertyarns * self.room_H)
            / len(self.cfg_yarn)
            )
        for ind, model in enumerate(self.yarn_models):
            #the source mass is what was present in the overlap before doing room model
            massyarnoverlapold = self.yarn_mass_overlap[ind]
            #the new mass there we approximate from concentration
            massyarnoverlapnew = massperyarn
            massyarndiff = massyarnoverlapnew - massyarnoverlapold
            #based on removed, we set a source term in the overlap zone of 
            # of the yarn
            model.source_overlap = massyarndiff / self.delta_t / model.areaextend
        #store masses
        self.__store_masses(self.tstep)

    def view_sol(self):
        #maxv = np.max(self.sol)
        #minv = np.min(self.sol)
        #print 'max', maxv, minv
        #self.plottimes = np.arange(self.times[0],self.times[-1]+1,self.plotevery)
        plotextra = False
        times = self.times(self.timesteps+1-((self.tstep - 1) % MAX_STORE_LENGTH+1),
                           self.timesteps)
        sol = self.solpart[:(self.tstep - 1) % MAX_STORE_LENGTH+1]
        extravals = self.cfg.get("plot.extra_time_room")
        if extravals:
            extravals = extravals.split("|")
            if len(extravals) == 3 and not eval(extravals[1]) == []:
                plotextra = True
        plt.ion()
        ind = 0
        for ind, interpdat in enumerate(self.plotdata):
            xval = self.x0[ind]
            cellstart, interpval = interpdat
            conc_in_point = interpval * sol[:, ind] + (1-interpval) * sol[:, ind+1]
            print ('conc in end point', conc_in_point[-1])
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
            plt.title('Concentration at position %g mm' % xval)
            if plotextra:
                plt.plot(eval(extravals[1]), eval(extravals[2]), extravals[0])
            plt.plot(times, conc_in_point)
            #plt.ylim(0, maxv*1.1)
            plt.plot(times, np.ones(len(times)) * self.saturation_conc, 'k--')
            plt.plot(times, np.ones(len(times)) * self.treshold, 'b--')
            plt.show()
            plt.savefig(utils.OUTPUTDIR + os.sep
                        + 'AIconc_%03.1f_mm' % xval + const.FIGFILEEXT)
        return ind

    def plot_room_sol(self, ind):
        print ('Generating fig of solution over the room domain')
        self.viewerplotcount = 0
        times = self.times(self.timesteps+1-((self.tstep - 1) % MAX_STORE_LENGTH+1),
                           self.timesteps)
        sol = self.solpart[:(self.tstep - 1) % MAX_STORE_LENGTH+1]
        minval = np.min(sol)
        maxv = np.max(sol)
        try:
            maxval = np.power(10., int(math.log10(maxv))+1)
        except ValueError:
            maxval = minval + 10
        plt.ion()
        if self.plotevery:
            for time, ssol in zip(times, sol):
                if self.viewerplotcount == 0:
                    print ('plotting for time', time)
                    plt.rc("font", family="serif")
                    plt.rc("font", size=10)
                    width = 4.5  #width in inches
                    height = 1.4 #height in inches
                    plt.rc("figure.subplot", left=(50/72.27)/width)
                    plt.rc("figure.subplot", right=(width-10/72.27)/width)
                    plt.rc("figure.subplot", bottom=(14/72.27)/height)
                    plt.rc("figure.subplot", top=(height-7/72.27)/height)
                    fig = plt.figure(ind)
                    plt.gca().set_xlabel('Position [mm]')
                    plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
                    plt.gca().set_ylim(minval, maxval)
                    #plt.gca().yaxis.set_major_formatter(pylab.FormatStrFormatter('%e'))
                    plt.title('Concentration in the room at t = %g s' % time)
                    plt.ioff()
                    lines = plt.plot(self.grid, ssol, 'r')
                    plt.draw()
                    try:
                        fig.canvas.flush_events()
                    except NotImplementedError:
                        pass
                    plt.ion()
                    plt.savefig(utils.OUTPUTDIR + os.sep 
                                    + 'AIconc_%08.1f_sec' % time + const.FIGFILEEXT)
                    #remove the line again
                    lines.pop(0).remove()
                self.viewerplotcount += 1 
                self.viewerplotcount = self.viewerplotcount % self.plotevery
        else:
            #plot last
            time = times[-1]
            ssol = sol[-1]
            print ('plotting for time', time)
            plt.rc("font", family="serif")
            plt.rc("font", size=10)
            width = 4.5  #width in inches
            height = 1.4 #height in inches
            plt.rc("figure.subplot", left=(50/72.27)/width)
            plt.rc("figure.subplot", right=(width-10/72.27)/width)
            plt.rc("figure.subplot", bottom=(14/72.27)/height)
            plt.rc("figure.subplot", top=(height-7/72.27)/height)
            fig = plt.figure(ind)
            plt.gca().set_xlabel('Position [mm]')
            plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
            plt.gca().set_ylim(minval, maxval)
            #plt.gca().yaxis.set_major_formatter(pylab.FormatStrFormatter('%e'))
            plt.title('Concentration in the room at t = %g s' % time)
            plt.ioff()
            lines = plt.plot(self.grid, ssol, 'r')
            plt.draw()
            try:
                fig.canvas.flush_events()
            except NotImplementedError:
                pass
            plt.ion()
            plt.savefig(utils.OUTPUTDIR + os.sep 
                            + 'AIconc_%08.1f_sec' % time + const.FIGFILEEXT)

    def view_sol_mass(self, ind):
        """
        Plot the evolution of the mass at current state of the solution
        """
        times = self.times(self.timesteps+1-((self.tstep - 1) % MAX_STORE_LENGTH+1),
                           self.timesteps)
        fignr = ind
        plt.ion()
        for ind, ymass in enumerate(self.yarnmass):
            plt.rc("font", family="serif")
            plt.rc("font", size=10)
            width = 4.5  #width in inches
            height = 1.4 #height in inches
            plt.rc("figure.subplot", left=(50/72.27)/width)
            plt.rc("figure.subplot", right=(width-10/72.27)/width)
            plt.rc("figure.subplot", bottom=(14/72.27)/height)
            plt.rc("figure.subplot", top=(height-7/72.27)/height)
            plt.figure(fignr)
            plt.gca().set_xlabel('Time [s]')
            plt.gca().set_ylabel('Mass [$\mu$g]')
            plt.title('Mass AC in yarn type %d' % ind)
            plt.plot(times, ymass[:(self.tstep-1) % MAX_STORE_LENGTH+1])
            plt.savefig(utils.OUTPUTDIR + os.sep 
                        + 'AImass_yarn_%d' % ind + const.FIGFILEEXT)
            fignr += 1

        plt.figure(fignr)
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Mass [$\mu$g]')
        plt.title('Mass AC in the bednet')
        #plt.plot([0,],[28.935,], 'r*')
        plt.plot(times, self.totyarnmass[:(self.tstep-1) % MAX_STORE_LENGTH+1])
        plt.savefig(utils.OUTPUTDIR + os.sep 
                        + 'AImass_bednet' + const.FIGFILEEXT)
        fignr += 1
        plt.figure(fignr)
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Mass [$\mu$g]')
        plt.title('Mass AC in the room')
        plt.plot(times, self.totroommass[:(self.tstep-1) % MAX_STORE_LENGTH+1])
        plt.savefig(utils.OUTPUTDIR + os.sep 
                        + 'AImass_in_room' + const.FIGFILEEXT)
        fignr += 1
        #plot to check mass conservation
        plt.figure(fignr)
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Mass [$\mu$g]')
        plt.title('Total Mass AC')
        plt.plot(times, self.totroommass[:(self.tstep-1) % MAX_STORE_LENGTH+1] \
                    + self.totyarnmass[:(self.tstep-1) % MAX_STORE_LENGTH+1])
        plt.savefig(utils.OUTPUTDIR + os.sep 
                        + 'AImass_total' + const.FIGFILEEXT)
        return fignr

    def init_room(self):
        #first we create the room mesh
        self.create_mesh()
        self.initial_room()
        self.init_yarn()
        if not self.initialized:
            self.solve_ode_init()
        #storage for mass values
        self.yarnmass = [None] * len(self.yarn_models)
        for ind in range(len(self.yarn_models)):
            self.yarnmass[ind] = np.empty(MAX_STORE_LENGTH, float)
        self.totyarnmass = np.empty(MAX_STORE_LENGTH, float)
        self.totroommass = np.empty(MAX_STORE_LENGTH, float)
        #compute the initial masses
        self.__store_masses(0)

    def __store_masses(self, ind):
        (yarnmass, yarnmassoverlap, totyarnmass, totyarnmassoverlap,
                roommass, roomoverlapmass) = self.calc_mass()
        #compute initial values
        for mind, val in enumerate(yarnmass):
            self.yarnmass[mind][ind%MAX_STORE_LENGTH] = val
        self.totyarnmass[ind%MAX_STORE_LENGTH] = np.sum(totyarnmass) #+ np.sum(totyarnmassoverlap)
        self.totroommass[ind%MAX_STORE_LENGTH] = roommass + roomoverlapmass

    def write_info(self):
        """
        Write generic info on the bednet
        """
        print ("\n\nINFO ON BEDNET")
        print ("**************")
        print  ("volume bednet = %g m^3, which means calculated porosity"\
                " %f mm^3 fabric/mm^3" \
                % (self.volbednet/1e9, self.fabporosity))
        print  ("surface bednet = %g m^2, which means calculated surface"\
                "mass %f gram/m^2" \
                % (self.surfbednet/1e6, (self.totyarnmass[0]/1e6)/(self.surfbednet/1e6)))
        print (" initial mass in bednet", self.totyarnmass[0]/1e6, "gram, room",\
                self.totroommass[0]/1e6, "gram")
        print (" number of yarns in fabric", "vertical", self.nvertyarns, \
                "horizontal", self.nhoryarns)
        print (" masses in the yarns ")
        
        for mind, val in enumerate(self.yarn_mass):
            print ("Yarn %d has initial mass AC %f" % (mind, val))
            for ind, models in enumerate(self.yarn_models[mind].fiber_models):
                print ('yarn cell', ind)
                for type, model in enumerate(models):
                    print ("fibertype %d: fibermass %f ; " % (type,
                            self.yarn_models[mind].fiber_mass[ind, type]))
                print (' ')
            print ("Blend in yarn is", self.yarn_models[mind].blend)
        print ("**************\n\n")
            #raw_input("Press key to start")

    def dump_sol(self, index):
        """ Dump solpart to file with extension index """
        for mind, val in enumerate(self.yarn_mass):
            for ind, models in enumerate(self.yarn_models[mind].fiber_models):
                for type, model in enumerate(models):
                   self.fibermass = self.yarn_models[mind].get_fiber_mass()
        times = self.times(index*MAX_STORE_LENGTH, self.tstep-1)
        timestxt = self.times(index*MAX_PLOT_LENGTH,self.tstep-1)
        newyarnmass = [0] * len(self.yarnmass)
        for ind in range(len(self.yarnmass)):
            newyarnmass[ind] = self.yarnmass[ind][:len(times)]
        np.savez(utils.OUTPUTDIR + os.sep + 'bednetroom1d_solpart_%05d.npz' % index,
                 times=times,
                 sol = self.solpart[:len(times)],
                 tresh_sat = [self.saturation_conc, self.treshold],
                 grid_cellcenters = self.grid,
                 fibermass = self.fibermass,
                 yarnmass = newyarnmass,
                 totyarnmass = self.totyarnmass[:len(times)],
                 totroommass = self.totroommass[:len(times)]
                 )
        #roomconc over time to textfile
        filedata= open(utils.OUTPUTDIR + os.sep + "roomconc" + ".txt",'w')
        filedata.write("conc in the room is %s" %(self.solpart[:len(timestxt)]) )
        filedata.close()
        #roomconc at the outermost left position over time to textfile
        self.roomconcleft = np.empty((len(timestxt),2),float)
        self.roomconcleft[:,0] = timestxt
        self.roomconcleft[:,1] = self.solpart[:len(timestxt),0]
        
        self.roomconcmiddle  = np.empty((len(timestxt),2),float)
        self.roomconcmiddle[:,0] = timestxt
        self.roomconcmiddle[:,1] = self.solpart[:len(timestxt),int(self.nr_cell/2)]
        
        self.roomconcright = np.empty((len(timestxt),2),float)
        self.roomconcright[:,0] = timestxt
        self.roomconcright[:,1] = self.solpart[:len(timestxt),-1]
        
        filedata= open(utils.OUTPUTDIR + os.sep + "roomconcLEFT" + ".txt",'w')
        for i in range(0,len(self.roomconcleft)):
            filedata.write("%.5f %.5f\n" % (self.roomconcleft[i,0],self.roomconcleft[i,1]))
        #filedata.write("conc at outermost LEFT in the room is %s" %(self.solpart[:len(times),0]) )
        filedata.close()
        #roomconc at the middle of the room over time to textfile
        filedata= open(utils.OUTPUTDIR + os.sep + "roomconcMIDDLE" + ".txt",'w')
        for i in range(0,len(self.roomconcmiddle)):
            filedata.write("%.5f %.5f\n" % (self.roomconcmiddle[i,0],self.roomconcmiddle[i,1]))
        #filedata.write("conc at outermost LEFT in the room is %s" %(self.solpart[:len(times),0]) )
        filedata.close()
        #roomconc at the outermost right position over time to textfile
        filedata= open(utils.OUTPUTDIR + os.sep + "roomconcRIGHT" + ".txt",'w')
        for i in range(0,len(self.roomconcright)):
            filedata.write("%.5f %.5f\n" % (self.roomconcright[i,0],self.roomconcright[i,1]))
        filedata.close()

    def run(self, wait=False):
        self.init_room()
        self.write_info()
        t = self.times(self.tstep+1)
        while t <= self.time_period+self.delta_t/10:
            self.solve_timestep(t)
            t = self.times(self.tstep+1)
        #we set tstep one after last value
        self.tstep += 1

        if INSPECT_MEM:
            import gc
            notr =gc.collect()
            print ('unreachable objects:')
            print (notr)
            notr =gc.collect()
            print ('unreachable objects:')
            print (notr)
            raw_input('press key')
            print ("Remaining Garbage")
            print (gc.garbage)
            raw_input('press key')
            
        #save solution to output file
        self.dump_sol(self.solstoreind)
        filedata= open(utils.OUTPUTDIR + os.sep + "yarnconccenter" + ".txt",'w')
        for i in range(len(self.yarnconc_center)):
            filedata.write("%.8f %.8f\n" % (self.yarnconc_center[i,0],self.yarnconc_center[i,1]))
        filedata.close()
        filedata= open(utils.OUTPUTDIR + os.sep + "yarnconcsurface" + ".txt",'w')
        for i in range(len(self.yarnconc_surface)):
            filedata.write("%.8f %.8f\n" % (self.yarnconc_surface[i,0],self.yarnconc_surface[i,1]))
        filedata.close()
        fignr = self.view_sol()
        fignr = self.view_sol_mass(fignr+1)
        self.plot_room_sol(fignr+1)

        for ymod in self.yarn_models:
            ymod.view_sol([ymod.step_old_time], [ymod.step_old_sol])
            for ind_cell, models in enumerate(ymod.fiber_models):
                for ftype, model in enumerate(models):
                    model.view_last_sol(" cell %d, type %d" % (ind_cell, ftype))
                break

        if wait:
            raw_input("Finished bednet run")