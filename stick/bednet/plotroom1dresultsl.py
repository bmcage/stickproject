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
    plotting for use in articles
    An array as saved by room1dmodel can be opened and replotted
"""
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
from __future__ import division, print_function
import os.path
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
from numpy import pi

FILETOLOAD = '/home/benny/stickproject/fabricmuslin_Deet.initest/bednetroom1d_sol.npz'

if not os.path.isfile(FILETOLOAD):
    print ('Unexisting file %s' % FILETOLOAD)
    sys.exit()
else:
    print ('Loading file %s'% FILETOLOAD)

with np.load(FILETOLOAD) as data:
    times = data['times']  #self.times
    sol = data['sol']  #self.sol
    tresh_sat = data['tresh_sat'] 
    saturation_conc = tresh_sat[0]
    treshold = tresh_sat[1]
    grid = data['grid_cellcenters'] #self.grid
    yarnmass = data['yarnmass']    #self.yarnmass
    totyarnmass = data['totyarnmass'] #self.totyarnmass
    totroommass = data['totroommass'] #self.totroommass

#determine at what distance in mm to plot concentration over time: 
x0 = [1]

plotdata = []
for xplot in x0:
    assert grid[0] < xplot < grid[-1], "%f < %f < %f "\
        "Not satisfied, observer out of domain" % (grid[0], 
                                            xplot, grid[-1])
    for ind, xcell in enumerate(grid):
        if xcell >= xplot:
            interpol_start = (xcell-xplot)/(grid[ind]-grid[ind-1])
            plotdata.append((ind-1, interpol_start))
            break

def view_sol(times, sol):
    ind = 0
    for ind, interpdat in enumerate(plotdata):
        xval = x0[ind]
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
        plt.plot(times, conc_in_point)
        #plt.ylim(0, maxv*1.1)
        plt.plot(times, np.ones(len(times)) * saturation_conc, 'k--')
        plt.plot(times, np.ones(len(times)) * treshold, 'b--')
        plt.draw()
        #plt.savefig('AIconc_%03.1f_mm' % xval + '.png')
    return ind

def view_sol_mass(ind):
    """
    Plot the evolution of the mass at current state of the solution
    """
    fignr = ind
    plt.ion()
    for ind, ymass in enumerate(yarnmass):
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
        print (len(times), len(ymass))
        plt.plot(times, ymass)
        fignr += 1

    plt.figure(fignr)
    plt.gca().set_xlabel('Time [s]')
    plt.gca().set_ylabel('Mass [$\mu$g]')
    plt.title('Mass AC in the bednet')
    #plt.plot([0,],[28.935,], 'r*')
    plt.plot(times, totyarnmass)
    fignr += 1
    plt.figure(fignr)
    plt.gca().set_xlabel('Time [s]')
    plt.gca().set_ylabel('Mass [$\mu$g]')
    plt.title('Mass AC in the room')
    plt.plot(times, totroommass)
    fignr += 1
    #plot to check mass conservation
    plt.figure(fignr)
    plt.gca().set_xlabel('Time [s]')
    plt.gca().set_ylabel('Mass [$\mu$g]')
    plt.title('Total Mass AC')
    plt.plot(times, totroommass+totyarnmass)
    return fignr

# we plot the data
plt.ion()
ind = view_sol(times, sol)
view_sol_mass(ind)

plt.show(block=True) #block=True)
