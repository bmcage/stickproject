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
Module holding a generic diffusion model for a yarn. 
"""
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
from __future__ import division
import os.path
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

HAVE_ODES = False
try:
    from scikits.odes import ode as sc_ode
    HAVE_ODES = True
except:
    print 'Could not load scikits.odes, odes solver not available'

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------
import stick.lib.utils.utils as utils
from stick.fiber.config import (METHOD, FLUX, TRANSFER, EVAP, EVAPYARN,
                    BOUND_TYPE, FIBER_FORM,
                    CIRCLE, ELLIPSE)

#-------------------------------------------------------------------------
#
#Fipy Imports
#
#-------------------------------------------------------------------------
from fipy import CylindricalGrid1D, Matplotlib1DViewer, CellVariable
from fipy.tools.dump import write as fipywrite

def Heaviside_oneside(val, control):
    """
    a Heaviside function of val, if control is positive, otherwise identity
    """
    if control < 0.:
        return 1.
    if val < 0.:
        return 0.
    return 1.
    #return val

def Heaviside_oneside_smoothed(val, control):
    """
    a smoothed Heaviside function of val, if control is positive, otherwise identity
    Smoothed in the sense:
    If val<0: 0
    if val-smoothwidth<0: val/smoothwidth
    if val> smoothwidth: 1
    """
    smoothwidth = 1e-6  #what would be good value here??
    if control < 0.:
        return 1.
    if val < 0.:
        return 0.
    if val < smoothwidth:
        #smoothing to avoid problems in ode solver
        return val/smoothwidth
    return 1.

#-------------------------------------------------------------------------
#
# DiffusionModel class 
#
#-------------------------------------------------------------------------
class FiberModel(object):
    """
    fiber1d.FiberModel is a special diffusion model for a single radial 
    symmetric fiber which is composed of a specific material (cotton, polyester,...)
    and has a number of coatings.
    
    The equation solved is \partial_t C = 1/r \partial_r (r D \partial_r C). 
    So the flux of rC over an interface is r D \partial_r C, or in other words,
    in radial components, the flux is r times the flux in x component
    
    If a flux boundary condition is given as F, it is interpreted as 
    requiring \partial_r C . n = F.
    This is done by setting the normal flux of r D \partial_r C = r D F
    
    If an evaporation boundary condition is given, we assume an outward flux
    at the surface of S h_lg (C_sat(T) - C_free) H(C - C_bo), so again, 
    radially, this is multiplied with r.
    """
    def __init__(self, config):
        """ 
        a config class must be passed in that contains the required settings
        
        boundary_left: the Neumann boundary condition for left side of 1D domain;
        boundary_right: the Neumann boundary condition for right side of 1D domain;
        transfer: the Robin BC for the right side of 1D domain
        initial_c: the initial concentration of DEET in the whole domain
        """
        self.cfg = config
        self.verbose = self.cfg.get('general.verbose')
        self.temp = 273.15 + 21 #temperature in Kelvin
        self.name = '0'

        self.solver = None
        self.method = self.cfg.get('general.method')
        if not (self.method in METHOD):
            print 'ERROR: unkown solution method %s' % self.method
            sys.exit(0)
        self.submethod = self.cfg.get('general.submethod')
        if not (self.submethod in METHOD[self.method][1]):
            print 'ERROR: unkown solution submethod %s' % self.submethod
            sys.exit(0)
        elif self.submethod in ['cvode', 'cvode_step'] and not HAVE_ODES:
            print 'ERROR: %s method not available, could not load scikits.odes' % self.submethod
            sys.exit(0)

        self.fiber_diff = self.cfg.get('fiber.internal_diffusion')
        if self.fiber_diff and (self.cfg.get('fiber.diffusion_coef') == 0.):
            raise Exception, 'Internal diffusion requires diff coeff > 0.'

        self.time_period = self.cfg.get('time.time_period')
        self.delta_t = self.cfg.get('time.dt')
        self.steps = int((self.time_period*(1.+self.delta_t*1e-6)) // self.delta_t)
        self.times = sp.linspace(0, self.time_period, self.steps + 1)
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        #set correct delta_t
        self.delta_t = self.times[1]-self.times[0]
        if self.verbose:
            print "Timestep used in fiber model:", self.delta_t
        #storage for output
        self.fiber_surface = sp.empty(len(self.times), float)
        #the radial flux at surface. Do times 2 pi radius() to obtain all flux
        self.flux_at_surface = sp.empty(len(self.times), float)
        
        #print 'the times', self.times
        #self.delta_t = 0.1#self.times[1] - self.times[0]
        #read the initial and boundary information for fiber
        self.n_edge = self.cfg.get('fiber.n_edge') #discretize the fiber radius
        self.bound_left = BOUND_TYPE[self.cfg.get('boundary.type_left')]
        self.bound_right = BOUND_TYPE[self.cfg.get('boundary.type_right')]
        self.boundary_fib_left = self.cfg.get('boundary.boundary_fib_left')
        self.boundary_fib_right = self.cfg.get('boundary.boundary_fib_right')
        self.boundary_transf_right = self.cfg.get('boundary.transfer_right')
        if self.bound_right == EVAP:
            self.evap_satconc = eval(self.cfg.get('boundary.evap_satconc'))
            self.evap_transfer = self.cfg.get('boundary.evap_transfer')
            self.evap_minbound = self.cfg.get('boundary.evap_minbound')
            self.out_conc = eval(self.cfg.get('boundary.out_conc'))
            self.only_compound = self.cfg.get('boundary.only_compound')
        #density in g/mm^3
        self.density_compound = self.cfg.get('compound.density') * 1e-3

        self.area_extend = self.cfg.get("fiber.extendarea")
        #data for stepwise operation
        self.initialized = False
        self.__userdata = None
        
        self.__Rf_pure = None
        self.__Rf = None

        self.plotevery = self.cfg.get("plot.plotevery")

    def set_userdata(self, data):
        """
        Read the data of the concentration in the void space and overwrite the 
        default value 
        """
        self.__userdata = data

    def get_userdata(self):
        return self.__userdata

    def set_areaextend(self, area):
        """
        Areaextend is based on calculation in the yarn model. With this method
        yarn model can set the area extend of this fiber model
        """
        self.area_extend = area

    def radius_pure(self):
        """
        method that returns the radius of the fiber seen as a circle,
        without coatings
        """
        if self.__Rf_pure:
            return self.__Rf_pure
        rad = self.cfg.get('fiber.radius_pure_fiber')
        form = FIBER_FORM[self.cfg.get('fiber.form')]
        if form == CIRCLE:
            pass
        elif form == ELLIPSE:
            #radius_pure_fiber is the long axis
            ecc = self.cfg.get('fiber.eccentricity')
            if not (ecc>0. and ecc <= 1.):
                raise Exception, 'Eccentricity must be between 0 and 1'
            #we determine radius of the fiber with same surface area
            #pi a b = pi r**2
            a = rad
            b = sp.sqrt(a**2 * (1-ecc**2))
            rad = sp.sqrt(a*b)
        else:
            raise Exception, 'Fiber form is not supported for a 1D fiber model'
        self.__Rf_pure = rad
        return self.__Rf_pure
    
    def radius(self):
        """ method that returns the total radius of the fiber
        """
        if self.__Rf:
            return self.__Rf
        rad = self.radius_pure()
        for i in range(self.cfg.get('fiber.nrlayers')):
            section = 'fiberlayer_%i' % i
            rad += self.cfg.get(section + '.thickness')
        self.__Rf = rad
        return self.__Rf

    def create_mesh(self):
        """
        Create a mesh for use in the model.
        We use an equidistant mesh
        grid: the space position of each central point for every cell element;
        """
        if not self.fiber_diff:
            n_edge = [0]
            self.diff_coef = [0.]
            self.diff_exp_fact = [0.]
            self.init_conc = [lambda x: 0.0 ]
            self.ind_first_zone = 1
            self.porosity = [0.]
        else:
            n_edge = [self.cfg.get('fiber.n_edge')]
            self.diff_coef = [self.cfg.get('fiber.diffusion_coef')]
            self.diff_exp_fact = [self.cfg.get('fiber.diffusion_polymer_exp_factor')]
            self.init_conc = [eval(self.cfg.get('fiber.init_conc')) ]
            self.ind_first_zone = 0
            self.porosity = [self.cfg.get('fiber.porosity_in')]
        self.nrlayers = self.cfg.get('fiber.nrlayers')
        self.surf_begin = [0.]
        self.surf = [self.radius_pure()]
        for i in range(self.nrlayers):
            section = 'fiberlayer_%i' % i
            n_edge += [self.cfg.get(section + '.n_edge')]
            self.surf_begin += [self.surf[-1]]
            self.surf += [self.surf[-1] + self.cfg.get(section + '.thickness')]
            self.diff_coef += [self.cfg.get(section + '.diffusion_coef')]
            self.diff_exp_fact += [self.cfg.get(section + '.diffusion_polymer_exp_factor')]
            self.init_conc += [eval(self.cfg.get(section + '.init_conc')) ]
            self.porosity += [self.cfg.get(section + '.porosity_layer')]
        #now take extend into account
        self.use_extend = self.cfg.get("fiber.useextension")

        #we now construct the full edge grid
        self.tot_edges = 0
        first = True
        for nr in n_edge:
            if nr == 0 and self.tot_edges == 0 and not first:
                print 'ERROR, no discretization points given'
                sys.exit(0)
            if not (nr == 0):
                first = False
                if self.tot_edges:
                    self.tot_edges += nr - 1
                else:
                    self.tot_edges = nr
        if self.use_extend:
            self.end_extend = np.sqrt(self.area_extend/np.pi + self.surf[-1]**2)
            sizegrid = (self.surf[-1]-self.surf_begin[-1])/n_edge[-1]
            self.nr_edge_extend = max(3,
                    int((self.end_extend - self.surf[-1])/sizegrid))
            self.nr_edge_extend = min(20, self.nr_edge_extend)
            n_edge += [self.nr_edge_extend]
            self.surf_begin += [self.surf[-1]]
            self.surf += [self.end_extend]
            self.diff_coef += [self.cfg.get('fiber.extenddiff')]
            self.diff_exp_fact += [0.]
            exconval = self.cfg.get('fiber.extendinit_conc')
            self.init_conc += [lambda x: exconval]
            self.porosity  += [1.]
            self.tot_edges_no_extend = self.tot_edges
            self.tot_edges += self.nr_edge_extend - 1
        else:
            self.tot_edges_no_extend = self.tot_edges
        self.grid_edge = sp.empty(self.tot_edges , float)
        left = 0.
        totnr = 0
        first = True
        for nr, right in zip(n_edge, self.surf):
            if nr:
                if first:
                    #left = 0.
                    self.grid_edge[totnr:totnr+nr] = sp.linspace(left, right, nr)
                    totnr += nr
                else:
                    self.grid_edge[totnr:totnr+nr-1] = sp.linspace(left, right, nr)[1:]
                    totnr += nr-1
                first = False
            left = right
        #construct cell centers from this
        self.grid = (self.grid_edge[:-1] + self.grid_edge[1:])/2.
        #obtain cell sizes
        self.delta_r = self.grid_edge[1:] - self.grid_edge[:-1]
        if self.verbose:
            pass
            #print 'mesh center', self.grid
            #print 'mesh edges ', self.grid_edge

        #create a fipy mesh for visualization and fipy computation
        self.mesh_fiber = CylindricalGrid1D(dx=tuple(self.delta_r))
        self.mesh_fiber_plot = CylindricalGrid1D(dx=tuple(self.delta_r[:self.tot_edges_no_extend-1]))
        self.mesh_fiber.periodicBC = False
        self.mesh_fiber = self.mesh_fiber + \
                        ((self.surf_begin[self.ind_first_zone],),)
        self.mesh_fiber_plot.periodicBC = False
        self.mesh_fiber_plot = self.mesh_fiber_plot + \
                        ((self.surf_begin[self.ind_first_zone],),)

    def initial_fiber(self):
        """ initial concentration over the domain"""
        if self.fiber_diff:
            #percentage of active component
            self.percentage_active = self.cfg.get('fiber.percentage_active')
        self.initial_c1 = sp.empty(self.tot_edges-1, float)
        self.diffusion_coeff = sp.empty(self.tot_edges-1, float)
        self.diffusion_exp_fact = sp.empty(self.tot_edges-1, float)
        self.porosity_domain = sp.empty(self.tot_edges - 1, float)
        st = 0
        surf = self.surf[st]
        for i, pos in enumerate(self.grid):
            while pos > surf:
                st += 1
                try:
                    surf = self.surf[st]
                except IndexError:
                    print st, self.surf, pos, surf
                    raise Exception, 'Index error'
            if st == 0:
                self.initial_c1[i] = self.init_conc[st](pos) * \
                                            (self.percentage_active / 100.)
            else:
                self.initial_c1[i] = self.init_conc[st](pos) 
            self.diffusion_coeff[i] = self.diff_coef[st]
            self.diffusion_exp_fact[i] = self.diff_exp_fact[st]
            self.porosity_domain[i] = self.porosity[st]
        self.initial_w1 = self.initial_c1 * self.grid
        self.volume = self.calc_volume()
        if self.verbose:
            print 'initial concentration', self.initial_c1
            print 'initial mass = ', self.calc_mass(self.initial_c1)

    def calc_mass(self, conc_r):
        """calculate the mass of component present given value in cell center
        This is given by 2 \pi int_r1^r2 n C(r)r dr
        
        conc_r: concentration in self.grid
        """
        return sp.sum(conc_r[:self.tot_edges_no_extend-1] * (
            sp.power(self.grid_edge[1:self.tot_edges_no_extend], 2) -
                sp.power(self.grid_edge[:self.tot_edges_no_extend-1], 2))
                * self.porosity_domain[:self.tot_edges_no_extend-1]) * sp.pi

    def calc_mass_surface(self, surface_conc_r):
        return surface_conc_r * self.porosity_domain[self.tot_edges_no_extend-2] * sp.pi * \
                (np.power(self.grid_edge[self.tot_edges_no_extend-1], 2)
                    - np.power(self.grid_edge[self.tot_edges_no_extend-2], 2))

    def calc_volume(self):
        """calculate the volume over which the compound can move. We have
           Cavg = mass/volume
        """
        return sp.sum((sp.power(self.grid_edge[1:self.tot_edges_no_extend],2) - 
                        sp.power(self.grid_edge[:self.tot_edges_no_extend-1],2)) 
                * self.porosity_domain[:self.tot_edges_no_extend-1]) * sp.pi 

    def set_outconc(self, conc):
        """
        set concentration in the extend area to a specific value
        Conc is C, solution stored is w = C*r
        """
        self.step_old_sol[self.tot_edges_no_extend-1:] = conc * self.grid[self.tot_edges_no_extend-1:]

    def _set_bound_flux(self, flux_edge, w_rep, t):
        """
        Method that takes BC into account to set flux on edge
        Data is written to flux_edge, w_rep contains solution in the cell centers
        This is the flux of w over radius! So times 2 pi R for full flux w, or
        for C times 2 pi.
        """
        if self.bound_left == FLUX:
            if abs(self.grid_edge[0]) < 1e-10 and not (self.boundary_fib_left == 0):
                print 'ERROR: Flux boundary to left, but at center of radial symmetry'
                sys.exit(0)
            if self.boundary_fib_left == 0:
                flux_edge[0] = 0
            else:
                # a FVM cannot do pure Neumann condition, instead we set the
                # FVM flux D \partial_x C with the value we want for \partial_x C.
                # eg: \partial_x C = F, set D \partial_x C = D F.
                # for radial coordinates, this flux times r.
                flux_edge[0] =  -(self.porosity_domain[0] * self.diffusion_coeff[0] *
                            sp.exp(-self.diffusion_exp_fact[0] * w_rep[0]/self.grid_edge[0]) \
                        ) * self.boundary_fib_left * self.grid_edge[0]
        else:
            print 'ERROR: boundary type left not implemented'
            sys.exit(0)
        #calculate normal val with w_rep = C*r, instead of C:
        conc_out = None
        if self.use_extend:
            conc_out = w_rep[self.tot_edges_no_extend-1] / \
                self.grid[self.tot_edges_no_extend-1]
        flux_edge[self.tot_edges_no_extend-1] = self._bound_flux_uR(
                w_rep[self.tot_edges_no_extend-2] /
                self.grid_edge[self.tot_edges_no_extend-1], t, conc_out)
        #and correct, flux needed is r times the x coord flux
        flux_edge[self.tot_edges_no_extend-1] = flux_edge[self.tot_edges_no_extend-1]\
                            * self.grid_edge[self.tot_edges_no_extend-1]
        if self.use_extend:
            flux_edge[-1] = 0.

    def _bound_flux_uR(self, conc_r, t, conc_out=None):
        """
        Calculate the flux on the surface from concentration on surface with
        the BC
        """
        if self.bound_right == FLUX:
            # a FVM cannot do pure Neumann condition, instead we set the
            # FVM flux - D \partial_x C with the value we want for \partial_x C.
            # eg: \partial_x C = F, set - D \partial_x C = - D F
            return -(self.porosity_domain[self.tot_edges_no_extend-1] \
                * self.diffusion_coeff[self.tot_edges_no_extend-1] *
                sp.exp(-self.diffusion_exp_fact[self.tot_edges_no_extend-1] * conc_r) \
                        ) * self.boundary_fib_right
        elif self.bound_right == TRANSFER:
            # a transfer coeff to the right, which is a given flux of
            # h_tf * C
            return self.boundary_transf_right * conc_r \
                             * self.porosity_domain[self.tot_edges_no_extend-1]
        elif self.bound_right == EVAP:
            # flux S h_lg (C_sat(T) - C_free) H(C - C_bo)
            #print 'the value from the function: out_conc', self.out_conc(0., 0.)
            #self.out_conc = eval(self.cfg.get('boundary.out_conc'))
            if self.use_extend:
                eCy = conc_out
            else:
                eCy = self.out_conc(t, self.get_userdata())
            eCs = self.evap_satconc(self.temp)
            #determine effective surface for evaporation: S=n C_Boundary/rho_compound
            S = 1.
            if not self.only_compound:
                # compound in a matrix, we need to correct effective surface. 
                S = self.porosity_domain[self.tot_edges_no_extend-2] \
                    * conc_r / self.density_compound
            #return evaporative law
            return (S 
                    * self.evap_transfer * (eCs - eCy)
                    * Heaviside_oneside_smoothed(conc_r - self.evap_minbound, 
                                        eCs - eCy)
                   )

##    NO LONGER SUPPORTED AS EXTENDAREA NOT TAKEN INTO ACCOUNT!
##    def _set_bound_fluxu(self, flux_edge, conc_r, t):
##        """
##        Method that takes BC into account to set flux on edge
##        Data is written to flux_edge, conc_r contains solution in the cell centers
##        """
##        if self.bound_left == FLUX:
##            flux_edge[0] = -self.boundary_fib_left * self.diffusion_coeff[0]
##        else:
##            print 'ERROR: boundary type left not implemented'
##            sys.exit(0)
##        flux_edge[-1] = self._bound_flux_uR(conc_r[-1], t)

    def f_conc1_ode(self, t, w_rep):
        self.f_conc1_odes(t, w_rep, self.__tmp_diff_w_t)
        return self.__tmp_diff_w_t

    def f_conc1_odes(self, t, w_rep, diff_w_t):
        grid = self.grid
        n_cell = len(grid)
        #Initialize the left side of ODE equations
        #initialize the flux rate on the edge with replace 'w'
        flux_edge = self.__tmp_flux_edge
        self._set_bound_flux(flux_edge, w_rep, t)
        #Diffusion coefficient changes with the concentration changing
        #calculate flux rate in each edge of the domain
        flux_edge[1:self.tot_edges_no_extend-1] = (
        -(self.porosity_domain[:self.tot_edges_no_extend-2] * self.diffusion_coeff[:self.tot_edges_no_extend-2] * 
                            sp.exp(-self.diffusion_exp_fact[:self.tot_edges_no_extend-2] * w_rep[:self.tot_edges_no_extend-2]/self.grid[:self.tot_edges_no_extend-2]) \
                         + self.porosity_domain[1:self.tot_edges_no_extend-1] * self.diffusion_coeff[1:self.tot_edges_no_extend-1] * 
                            sp.exp(-self.diffusion_exp_fact[1:self.tot_edges_no_extend-1] * w_rep[1:self.tot_edges_no_extend-1]/self.grid[1:self.tot_edges_no_extend-1]))/2.\
                        * self.grid_edge[1:self.tot_edges_no_extend-1] \
                        * (w_rep[1:self.tot_edges_no_extend-1]/self.grid[1:self.tot_edges_no_extend-1] - w_rep[:self.tot_edges_no_extend-2]/self.grid[:self.tot_edges_no_extend-2])\
                        / ((self.delta_r[:self.tot_edges_no_extend-2] + self.delta_r[1:self.tot_edges_no_extend-1])/2.)
                        )
        if self.use_extend:
            #also flux in outside part
            flux_edge[self.tot_edges_no_extend:-1] = -(self.porosity_domain[self.tot_edges_no_extend-1:-1] * self.diffusion_coeff[self.tot_edges_no_extend-1:-1] \
                             + self.porosity_domain[self.tot_edges_no_extend:] * self.diffusion_coeff[self.tot_edges_no_extend:])/2.\
                            * self.grid_edge[self.tot_edges_no_extend:-1] \
                            * (w_rep[self.tot_edges_no_extend:]/self.grid[self.tot_edges_no_extend:] - w_rep[self.tot_edges_no_extend-1:-1]/self.grid[self.tot_edges_no_extend-1:-1])\
                            / ((self.delta_r[self.tot_edges_no_extend:] + self.delta_r[self.tot_edges_no_extend-1:-1])/2.)

        diff_w_t[:] = (flux_edge[:-1] - flux_edge[1:]) / self.delta_r[:] / self.porosity_domain[:]

##    def f_conc1_odeu(self, t, conc_r):
##        grid = self.grid
##        n_cell = len(grid)
##        #Initialize the left side of ODE equations
##        diff_u_t = self.__tmp_diff_w_t
##        #initialize the flux rate on the edge with replace 'w'
##        flux_edge = self.__tmp_flux_edge
##        self._set_bound_fluxu(flux_edge, conc_r, t)
##        #Diffusion coefficient changes with the concentration changing
##        #calculate flux rate in each edge of the domain
##        flux_edge[1:-1] = -(self.porosity_domain[:-1] * self.diffusion_coeff[:-1] * 
##                            sp.exp(-self.diffusion_exp_fact[:-1] * conc_r[:-1]) \
##                         + self.porosity_domain[1:] * self.diffusion_coeff[1:] * 
##                            sp.exp(-self.diffusion_exp_fact[1:] * conc_r[1:]))/2.\
##                        * self.grid_edge[1:-1] \
##                        * (conc_r[1:] - conc_r[:-1])\
##                        / ((self.delta_r[:-1] + self.delta_r[1:])/2.)
##        diff_u_t[:] = 2.*((flux_edge[:-1]-flux_edge[1:])
##                        /(self.grid_edge[1:]**2-self.grid_edge[:-1]**2)
##                        / self.porosity_domain[:])
##        return diff_u_t

    def solve_odes_init(self, clearmem = False):
        """
        Initialize the cvode solver
        """
        if not HAVE_ODES:
            raise Exception, 'Not possible to solve with given method, scikits.odes not available'
        self.initial_t = self.times[0]
        if clearmem:
            del self.times
        self.step_old_time = self.initial_t
        self.step_old_sol = self.initial_w1
        #data storage
        self.ret_y = sp.empty(len(self.initial_c1), float)

        n_cell = len(self.grid)
        self.__tmp_diff_w_t = sp.empty(n_cell, float)
        self.__tmp_flux_edge = sp.empty(n_cell+1, float)
        self.tstep = 0
        self.solve_odes_reinit()
        self.initialized = True

    def solve_odes_reinit(self):
        """
        Reinitialize the cvode solver to start again
        """
        #self.initial_t = self.times[0]
        if self. solver is None:
            self.solver = sc_ode('cvode', self.f_conc1_odes,
                            max_steps=50000, lband=1, uband=1)
        self.solver.init_step(self.step_old_time, self.step_old_sol)

    def solve_odes(self, run_per_step = None, viewend = True):
        self.initial_w1 = self.initial_c1 * self.grid
        #data storage, will give outofmem for long times!
        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
        tstep = 0
        self.conc1[tstep][:] = self.initial_c1
        for time in self.times[1:]:
            flag, realtime = self.solver.step(time, self.conc1[tstep+1])
            if flag != 0:
                print 'ERROR: unable to compute solution, flag', flag
                break
            if self.verbose:
                print 'INFO: fibermodel at t = ', realtime
            tstep += 1
            self.conc1[tstep][:] = self.conc1[tstep][:] / self.grid[:]
            conc_out = None
            if self.use_extend:
                conc_out = self.conc1[tstep][self.tot_edges_no_extend-1] 
            self.fiber_surface[tstep] = self.conc1[tstep][self.tot_edges_no_extend-2]
            self.flux_at_surface[tstep] = self._bound_flux_uR(
                        self.conc1[tstep][self.tot_edges_no_extend-2], time, conc_out)
            if run_per_step:
                run_per_step(self, time, self.conc1[tstep])

            #print 'mass = ', self.calc_mass(self.conc1[tstep])
        if viewend:
            self.view_sol(self.times, self.conc1)

    def do_step_odes(self, stoptime, needreinit=True):
        """
        Solve the fibermodel up to stoptime, continuing from the present
        state, return time and r * concentration after step
        It is needed that run_init and solve_init method have been called
        before calling this method.
        
        if needreinit = True, the solver is initialized first with the
            data present in step_old_time amd step_old_sol
            
        Return: concentration over the grid
        """
        assert stoptime > self.step_old_time, "%f > %f" % (stoptime, self.step_old_time)
        if not self.initialized:
            raise Exception, 'Solver ode not initialized'
        if needreinit:
            self.solve_odes_reinit()
        else:
            self.solver.set_options(tstop=stoptime)
        compute = True
        #even is step is large, we don't compute for a longer time than delta_t
        t = self.step_old_time
        while compute:
            t +=  self.delta_t
            if  t >= stoptime - self.delta_t/100.:
                t = stoptime
                compute = False
            flag, realtime = self.solver.step(t, self.ret_y)
            if flag < 0:
                raise Exception, 'could not find solution, flag %d' % flag
        
        self.step_old_time = realtime
        self.step_old_sol = self.ret_y
        self.surface_conc = self.step_old_sol[self.tot_edges_no_extend-2] / \
                                self.grid[self.tot_edges_no_extend-2]
        assert np.allclose(realtime, stoptime, atol=1e-6, rtol=1e-6)
        return realtime, self.ret_y / self.grid

##    def solve_ode_init(self):
##        """
##        Initialize the ode solver
##        """
##        self.initial_t = self.times[0]
##        self.step_old_time = self.initial_t
##        self.step_old_sol = self.initial_w1
##        #data storage
##        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
##        self.conc1[0][:] = self.initial_c1
##        n_cell = len(self.grid)
##        self.__tmp_diff_w_t = sp.empty(n_cell, float)
##        self.__tmp_flux_edge = sp.empty(n_cell+1, float)
##        
##        self.tstep = 0
##        self.solve_ode_reinit()
##        self.initialized = True
##
##    def solve_ode_reinit(self):
##        """
##        Reinitialize the ode solver to start again
##        """
##        self.initial_t = self.times[0]
##        self.solver = ode(self.f_conc1_ode).set_integrator('vode', 
##                            method = 'bdf',
##                            nsteps=50000)
##        self.solver.set_initial_value(self.step_old_sol, self.step_old_time)

##    def do_step_ode(self, step, needreinit=True):
##        """
##        Solve the fibermodel for one step, continuing from the present
##        state, return the time and r* concentration after step
##        It is needed that run_init and solve_init method have been called
##        before calling this method.
##        
##        """
##        if not self.initialized:
##            raise Exception, 'Solver ode not initialized'
##        if needreinit:
##            self.solve_ode_reinit()
##        else:
##            raise Exception, 'odew method requires reinit for stepwize operation'
##        curt = self.solver.t
##        #print 'curt value', curt
##        while self.solver.t < curt + step - self.delta_t /10. and \
##        self.solver.t < self.step_old_time + step - self.delta_t /10.:
##            self.solver.integrate(self.solver.t + self.delta_t)
##            if not self.solver.successful():
##                raise Exception, 'could not find solution'
##
##        #return the concentration after step
##        self.solver.integrate(curt + step)
##        if not self.solver.successful():
##            raise Exception, 'could not find solution'
##        self.step_old_time += step
##        self.step_old_sol = self.solver.y
##        assert abs(self.solver.t - self.step_old_time) < 1e-15, "%f %f %g" % (self.solver.t, self.step_old_time, self.solver.t- self.step_old_time)
##        return self.solver.t, self.solver.y
        
##    def solve_ode(self, run_per_step = None, viewend = True):
##        self.solve_ode_init()
##        endT = self.times[-1]
##        self.initial_w1 = self.initial_c1 * self.grid
##        tstep = 0
##        self.conc1[tstep][:] = self.initial_c1
##        while self.solver.successful() and self.solver.t < endT - self.delta_t /10.:
##            self.solver.integrate(self.solver.t + self.delta_t)
##            if self.verbose:
##                print 'INFO: fibermodel at t = ', self.solver.t
##            tstep += 1
##            self.conc1[tstep][:] = self.solver.y / self.grid
##            self.fiber_surface[tstep] = self.conc1[tstep][self.tot_edges_no_extend-2]
##            conc_out = None
##            if self.use_extend:
##                conc_out = self.conc1[tstep][self.tot_edges_no_extend-1] 
##            self.flux_at_surface[tstep] = self._bound_flux_uR(
##                        self.conc1[tstep][-1], self.solver.t + self.delta_t, conc_out)
##            #print 'mass = ', self.calc_mass(self.conc1[tstep])
##            if run_per_step:
##                run_per_step(self, self.solver.t, self.conc1[tstep])
##        if viewend:
##            self.view_sol(self.times, self.conc1)

    def solve_simple(self):
        """
        A simplified solution of the problem in which diffusion is neglected,
        so with M = C * V the mass: dM/dt = - flux_out, so
        M(t_e) = M(0) - int(flux_out, t=0..t_e)
        For the transfer condition we have
            dM/dt = - k M
              M(t) = M(0)*exp(-k*t/V)
        If also fluxout is removed as flux:
            M(t) = fluxout/k+(M0-fluxout/k)*exp(-k*t)
        fluxout is the amount over an edge, so it is 2 pi fluxout(r) * r_surf porosity_edge
        """
        V = self.volume
        M0 = self.simple_sol[0]
        k = 0.
        flux_out = 0.
        evap = False
        if self.bound_left == FLUX:
            flux_out += -2*sp.pi * self.grid_edge[-1]*self.porosity_domain[-1] \
                        * self.boundary_fib_left
        else:
            raise Exception, 'ERROR: boundary type left not implemented'
        if self.bound_right == FLUX:
            flux_out += 2*sp.pi * self.grid_edge[-1]*self.porosity_domain[-1] \
                        * self.boundary_fib_right 
        elif self.bound_right == TRANSFER:
            #flux in C is dC/dr=kC = k M/V, so flux in mass is this integrated
            # over edge, which gives 2 pi porosity_edge r_edge / V * k * M
            # which gives us the rate k in terms of M of:
            k = 2 * sp.pi * self.boundary_transf_right * self.grid_edge[-1] \
                    * self.porosity_domain[-1] / V
        elif self.bound_right == EVAP:
            coeffevap = 2 * sp.pi * self.grid_edge[-1] \
                        * self.porosity_domain[-1] \
                        * self.evap_transfer
            satevap = self.evap_satconc(self.temp)
            flux_outevap = lambda M, t:  coeffevap \
                                * (satevap - self.out_conc(t, self.get_userdata())) \
                                * Heaviside_oneside(M/V-self.evap_minbound, 
                                    satevap - self.out_conc(t, self.get_userdata()))
        else:
            raise Exception, 'ERROR: boundary type right not implemented'
        tstep = 0
        #a check to avoid errors by too large timestep!
        if self.bound_right == EVAP:
            first_step = self.delta_t * (flux_outevap(self.simple_sol[0], 
                                            self.times[0]) + flux_out)
            if abs(first_step/self.simple_sol[0]) > 0.05:
                print 'start mass:', self.simple_sol[0], 'change:', first_step
                raise Exception, 'ERROR, reduce time step of fiber model!' +\
                            ' Change is mass over first time step > 5%'
        for time in self.times[1:]:
            tstep += 1
            if self.bound_right == EVAP:
                #forward Euler on the dM/dt = -fluxout equation
                self.simple_sol[tstep] = self.simple_sol[tstep-1] - \
                    self.delta_t * (flux_outevap(self.simple_sol[tstep-1], 
                                            self.times[tstep-1])
                                    + flux_out)
            else:
                if k:
                    self.simple_sol[tstep] = flux_out/k+(M0-flux_out/k)*np.exp(-k*time)
                else:
                    self.simple_sol[tstep] = flux_out * time + M0
        
            #convert the mass to the average concentration valid over domain
            self.fiber_surface[tstep] = self.simple_sol[tstep] / V
            conc_out = None
            if self.use_extend:
                raise NotImplementedError, 'out conc needed, not known'
            self.flux_at_surface[tstep] = self._bound_flux_uR(
                                            self.fiber_surface[tstep], time, conc_out)
        if self.cfg.get('plot.plotmass'):
            self.view_time(self.times, self.simple_sol, 'Mass in cross section fiber')

##    def solve_odeu(self):
##        self.initial_t = self.times[0]
##        endT = self.times[-1]
##        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
##        r = ode(self.f_conc1_odeu).set_integrator('vode', method = 'bdf',nsteps = 10000)
##        r.set_initial_value(self.initial_c1, self.initial_t)#.set_f_params(2.0)
##        tstep = 0
##        self.conc1[tstep][:] = self.initial_c1
##        while r.successful() and r.t < endT - self.delta_t /10.:
##            r.integrate(r.t + self.delta_t)
##            tstep += 1
##            self.conc1[tstep][:] = r.y 
##            self.fiber_surface[tstep] = self.conc1[tstep][-1]
##            self.flux_at_surface[tstep] = self._bound_flux_uR(
##                                    self.conc1[tstep][-1], r.t + self.delta_t)
##            #print 'mass = ', self.calc_mass(self.conc1[tstep])
##        self.view_sol(self.times, self.conc1)

##    def solve_fipy(self):
##        #using fipy to solve 1D problem in fiber
##        self.solution_fiber = CellVariable(name = "fiber concentration", 
##                                mesh = self.mesh_fiber,
##                                value = self.initial_c1 * self.porosity_domain, hasOld = 1)
##        if self.plotevery:
##            self.viewer =  Viewer(vars = self.solution_fiber / self.porosity_domain, datamin=0., datamax= 1.1)
##        self.conc1 = np.empty((len(self.times), len(self.initial_c1)), float)
##
##        if self.bound_left == FLUX and self.bound_right == FLUX:
##            self.BCs_fiber = (FixedFlux(faces = self.mesh_fiber.getFacesRight(), 
##                                    value = self.boundary_fib_right),
##                              FixedFlux(faces = self.mesh_fiber.getFacesLeft(), 
##                                    value = -self.boundary_fib_left))
##        elif self.bound_left == FLUX and self.bound_right == TRANSFER:
##            self.BCs_fiber = (FixedFlux(faces = self.mesh_fiber.getFacesRight(), 
##                                       value = self.boundary_transf_right \
##                                            * self.solution_fiber.getFaceValue()),
##                              FixedFlux(faces = self.mesh_fiber.getFacesLeft(), 
##                                    value = -self.boundary_fib_left))
##        else:
##            print 'ERROR: boundary type left not implemented'
##            sys.exit(0)
##        ## TODO: in following diffusion is given as cellvariable, would a 
##        ##      FaceVariable already not be better?
##        self.eqX_fiber = TransientTerm() == DiffusionTerm(coeff = 
##                            self.diffusion_coeff * 
##                            sp.exp(-self.diffusion_exp_fact * self.solution_fiber / self.porosity_domain)*
##                            self.porosity_domain) 
##        tstep = 0
##        self.conc1[tstep][:] = self.initial_c1[:]
##        for time in self.times[1:]:
##            self.solve_fipy_sweep()
## ##           if self.viewer is not None:
## ##               self.viewer.plot()
##                #raw_input("please<return>.....")
##            tstep += 1
##            self.conc1[tstep][:] = self.solution_fiber.getValue()
##            self.fiber_surface[tstep] = self.conc1[tstep][-1]
##            self.flux_at_surface[tstep] = self._bound_flux_uR(
##                                                self.conc1[tstep][-1], time)
##            #if time == 200.0:
##            #    dump.write({'space_position': self.grid, 'conc1': self.conc1[tstep][:]},
##            #            filename = utils.OUTPUTDIR + os.sep + 'fipy_t1.gz', extension = '.gz')
##            #    print 'finish file'
##            #print 'mass = ', self.calc_mass(self.conc1[tstep])
##
##    def solve_fipy_sweep(self):
##        res = 1e+1
##        while res > 1e-8:
##            res = self.eqX_fiber.sweep(var = self.solution_fiber,
##                                        dt = self.delta_t,
##                                        boundaryConditions = self.BCs_fiber)
##        self.solution_fiber.updateOld()

    def solve_simple_init(self):
        """
        Initialize the simple solver
        """
        self.initial_t = self.times[0]
        self.step_old_time = self.initial_t
        self.step_old_sol = self.initial_w1
        self.simple_sol = np.empty(len(self.times), float)
        self.simple_sol[0] = self.calc_mass(self.initial_c1)
        self.tstep = 0
        self.initialized = True

    def solve(self):
        """
        Solve the diffusion process in the fiber. 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. The mesh in this 1-D condition is 
        uniform
        """
        
        def run_every_step(object, time, conc):
            if self.plotevery:
                if object.viewerplotcount == 0:
                    object.solution_view.setValue(conc[:self.tot_edges_no_extend-1])
                    object.viewer.axes.set_title('time %s' %str(time))
                    object.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'conc%s.png' % str(int(10*time)))
                object.viewerplotcount += 1
                object.viewerplotcount = self.viewerplotcount % self.plotevery

        if self.method == 'FVM':
            if self.submethod == 'fipy':
                raise NotImplementedError, 'this option %s is no longer supported' % self.submethod
                self.solve_fipy()
            elif self.submethod == 'cvode':
                self.solve_odes()
            elif self.submethod == 'cvode_step':                    
                #self.solve_odes_init()
                if self.plotevery:
                    self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber,
                            value = self.conc1[0][:self.tot_edges_no_extend-1])
                    self.viewer =  Matplotlib1DViewer(vars = self.solution_view, 
                                        datamin=0., 
                                        datamax=1.2 * self.conc1[0][:self.tot_edges_no_extend-1].max())
                    self.viewer.axes.set_title('time 0.0')
                    self.viewer.plot()
                    self.viewerplotcount = 1
                self.solve_odes(run_per_step =run_every_step, viewend=False)

            elif  self.submethod == 'odew':
                raise NotImplementedError, 'this option %s is no longer supported' % self.submethod
                self.solve_ode()
            elif  self.submethod == 'odew_step':
                raise NotImplementedError, 'this option %s is no longer supported' % self.submethod
                self.solve_ode_init()
                if self.plotevery:
                    self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber,
                            value = self.conc1[0][:])
                    self.viewer =  Matplotlib1DViewer(vars = self.solution_view, 
                                        datamin=0., 
                                        datamax=1.2 * self.conc1[0].max())
                    self.viewer.axes.set_title('time 0.0')
                    self.viewer.plot()
                    self.viewerplotcount = 1
                self.solve_ode(run_per_step =run_every_step, viewend=False)
            elif self.submethod == 'odeu':
                raise NotImplementedError, 'this option %s is no longer supported' % self.submethod
                self.solve_odeu()
            if self.verbose:
                print 'end mass = ', self.calc_mass(self.conc1[-1])
        elif self.method == 'SIMPLE':
            self.solve_simple()
        else:
            raise NotImplementedError, 'Method %s is not implemented' % self.method
        if self.verbose:
            print 'Finished the fiber calculation'
        if self.cfg.get('plot.plotflux'):
            self.view_time(self.times, self.flux_at_surface, 'Flux of DEET ($\mathrm{mg\cdot cm/s}$)')

    def do_step(self, stoptime, needreinit=True):
        """
        Solve the diffusion process in the fiber up to stoptime 
        &C/&t = 1/r * &(Dr&C/&r) / &r
        The diffusion coefficient is constant. The finite volume method is used to
        discretize the right side of equation. 
        The resulting time and r*concentration is returned
        """
        if  self.submethod in ['cvode', 'cvode_step']:
            return self.do_step_odes(stoptime, needreinit)
        #elif  self.submethod in ['odew', 'odew_step']:
        #    res = self.do_step_ode(step)
        else:
            raise Exception, 'Not supported option %s' % self.submethod

    def solve_init(self, clearmem = False):
        """
        Initialize the solvers so they can be solved stepwize
        """
        if self.method == 'FVM':
            if self.submethod == 'fipy':
                raise NotImplementedError, 'this option %s is no longer supported' % self.submethod
                self.solve_fipy()
            elif  self.submethod in ['cvode', 'cvode_step']:
                self.solve_odes_init(clearmem=clearmem)
            elif  self.submethod in ['odew', 'odew_step']:
                raise NotImplementedError, 'this option %s is no longer supported' % self.submethod
                self.solve_ode_init()
            elif self.submethod == 'odeu':
                raise NotImplementedError, 'this option %s is no longer supported' % self.submethod
                self.solve_ode_init()
            else:
                raise NotImplementedError
        elif self.method == 'SIMPLE':
            self.solve_simple_init()
        else:
            raise NotImplementedError, 'Method %s is not implemented' % self.method

    def view_sol(self, times, conc):
        """
        Show the solution in conc with times.
        conc[i][:] contains solution at time times[i]
        """
        self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber_plot,
                            value = conc[0][:self.tot_edges_no_extend-1])
        name = self.solution_view.name                    
        if self.plotevery:
            self.viewer =  Matplotlib1DViewer(vars = self.solution_view, datamin=0.,
                                datamax=conc.max()+0.20*conc.max())
        self.viewerplotcount = 0
        for time, con in zip(times[1:], conc[1:][:self.tot_edges_no_extend-1]):
            self.solution_view.setValue(con)
            if self.plotevery:
                self.viewer.axes.set_title('Fiber Conc vs radius at time %s' %str(time))
                self.viewer.axes.set_xlabel('Radius')
                self.viewer.axes.set_ylabel('Conc')
                if self.viewerplotcount == 0:
                   self.viewer.plot(filename=utils.OUTPUTDIR + os.sep + 'fiber%sconc%08.4f.png' % (name,time))
                #else:
                #    self.viewer.plot()
                    
                self.viewerplotcount += 1
                self.viewerplotcount = self.viewerplotcount % self.plotevery   

    def view_last_sol(self, title, time=None, conc=None):
        if conc is None:
            conc = self.step_old_sol / self.grid
        if time is None:
            time = self.step_old_time
        self.solution_view = CellVariable(name = "fiber concentration", 
                            mesh = self.mesh_fiber_plot,
                            value = conc[:self.tot_edges_no_extend-1])
        self.viewer =  Matplotlib1DViewer(vars = self.solution_view, datamin=0.,
                                datamax=conc.max()+0.20*conc.max())
        self.viewer.axes.set_title('Fiber Conc vs radius at time %s' %str(time) + 
                                    title)
        self.viewer.axes.set_xlabel('Radius')
        self.viewer.axes.set_ylabel('Conc')
        self.viewer.plot(filename=utils.OUTPUTDIR + os.sep +
                            'fiber_fiberconc_%08.4f_sec_%s.png' % 
                            (time, title.replace(',','').replace(' ','_')))

    def view_time(self, times, conc, title=None):
        draw_time = times/(3600.*24.*30.) # convert seconds to months
        draw_conc = conc *1.0e4
        plt.figure(num=None)
        plt.plot(draw_time, draw_conc, '-', color = 'red')
        plt.xlabel('Time (month)')
        plt.ylabel(title)
        plt.show()

    def dump_solution(self): 
        """write out the solution to disk for future use"""
        if self.method == 'FVM':
            fipywrite({'space_position': self.grid, 'conc': self.conc1},
                filename = utils.OUTPUTDIR + os.sep + 'sol_%s.gz' % self.submethod,
                extension = '.gz')
        fipywrite({'time_step':self.times, 'flux': self.flux_at_surface},
            filename = utils.OUTPUTDIR + os.sep + 'flux_boundary', 
            extension = '.gz')

    def run_init(self):
        self.create_mesh()
        self.initial_fiber()

    def run(self, wait=False, output=False):
        self.run_init()
        if not self.initialized:
            self.solve_init()
        mass_start = self.calc_mass(self.step_old_sol / self.grid)
        self.solve()
        mass_end = self.calc_mass(self.step_old_sol / self.grid)

        print 'mass start:', mass_start, ', mass end:', mass_end
        if output:
            self.dump_solution()
        if wait:
            raw_input("Finished fiber1d run")

    def __del__(self):
        if self.method == 'FVM' and self.submethod in ['cvode', 'cvode_step']:
            #remove the memory
            if self.solver:
                print 'del self.solver'
                del self.solver
            self.solver = None
