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

#BASEDIR = '/Users/Tine/stickproject/'
BASEDIR = '/home/benny/stickproject/'
PROBS = False  #set tot False if only one problem
#PROBTOLOAD = 'fabric.ini'
#PROBTOLOAD = 'fabricbednetY335_Deet.ini_50nmol8hour'
PROBTOLOAD = 'fabricbednetY335_replenishtest.ini'
#PROBTOLOAD = 'fabricbednetY335_releasetest.inirelease'
#ROBTOLOAD = 'fabric.inihigh'
#PROBTOLOAD = 'fabric.inilow'

#
GRAY = True
#all problems must be over the same grid !
SHOWLEGEND = True
if not PROBS:
    SHOWLEGEND = False
#PROBSTOLOAD = ['fabricbednetY335_Deet.ini_50nmol8hour', 
#    'fabricbednetY335_Deet.ini_100nmol8hour', 
#    'fabricbednetY335_Deet.ini_200nmol8hour', 
#    'fabricbednetY335_Deet.ini_400nmol8hour']
#LABELS = ['50 nmol', '100 nmol', '200 nmol', '400 nmol']
PROBSTOLOAD = ['fabric.inilow', 'fabric.inilowblock10',
               'fabric.inilowblockpoly', 'fabric.inilowblockpolyonlycomp']
LABELS = ['Case 1', 'Case 2', 'Case 3', 'Case 4']
#PROBSTOLOAD = ['fabricmuslin_Deet.ini_25nmol2min_b',
#    'fabricmuslin_Deet.ini_20nmol2min_b', 'fabricmuslin_Deet.ini_15nmol2min_b',
#    'fabricmuslin_Deet.ini_10nmol2min_b', 'fabricmuslin_Deet.ini_05nmol2min_b']
##LABELS = ['25 nmol', '20 nmol', '15 nmol', '10 nmol', '5 nmol']
ARG = '/bednetroom1d_solpart_%05d.npz'
INDEX = range(29) #range(4) # what dumped data to load
EVERY = 60 #1    # what time data to skip to reduce plotting time

#determine at what distance in mm to plot concentration over time: 
#x0 = [1, 5, 10, 500]
x0 = [1, 4]

double = True

if not PROBS:
    PROBSTOLOAD = [PROBTOLOAD]


def setAxLinesBW(ax):
    """
    Take each Line2D in the axes, ax, and convert the line style to be 
    suitable for black and white viewing.
    """
    MARKERSIZE = 3

    COLORMAP = {
        'b': {'marker': None, 'dash': (None,None)},
        'g': {'marker': None, 'dash': [3,3]},
        'r': {'marker': None, 'dash': [5,3,1,3]},
        'c': {'marker': None, 'dash': [1,3]},
        'm': {'marker': None, 'dash': [5,2,5,2,5,10]},
        'y': {'marker': None, 'dash': [5,3,1,2,1,10]},
       # 'k': {'marker': 'o', 'dash': (None,None)} #[1,2,1,10]}
        }

    extra = []
    if not ax.get_legend() is None:
        extra = ax.get_legend().get_lines()
    for line in ax.get_lines() + extra:
        origColor = line.get_color()
        line.set_color('black')
        if origColor in COLORMAP:
            line.set_dashes(COLORMAP[origColor]['dash'])
            line.set_marker(COLORMAP[origColor]['marker'])
            line.set_markersize(MARKERSIZE)

def setFigLinesBW(fig):
    """
    Take each axes in the figure, and for each line in the axes, make the
    line viewable in black and white.
    """
    for ax in fig.get_axes():
        if not ax is None:
            setAxLinesBW(ax)


prob = 0
ltimes = []
ltotyarnmass = []
ltotroommass = []
lyarnmass = []
lfconc_sta = []
lfconc_mid = []
lfconc_end = []
lsol = []
for PROBTOLOAD in PROBSTOLOAD:
  times = np.empty(0, float)
  sol = None
  yarnmass = None
  fiberconc_center = np.empty(0,float)
  fiberconc_middle = np.empty(0,float)
  fiberconc_surface = np.empty(0,float)
  fibertimes = np.empty(0, float)
  totyarnmass = np.empty(0, float)
  totroommass = np.empty(0, float)
  for index in INDEX:
    firstfile = BASEDIR + PROBTOLOAD + ARG % index
    if not os.path.isfile(firstfile):
        print ('Unexisting file %s' % firstfile)
        sys.exit()
    else:
        print ('Loading file %s'% firstfile)
    
    with np.load(firstfile) as data:
        times = np.append(times, data['times'][::EVERY])  #self.times
        if sol is None:
            sol = data['sol'][::EVERY]
        else:
            sol = np.append(sol, data['sol'][::EVERY], axis=0)  #self.sol
        tresh_sat = data['tresh_sat'] 
        saturation_conc = tresh_sat[0]
        treshold = tresh_sat[1]
        grid = data['grid_cellcenters'] #self.grid
            
        if yarnmass is None:
            yarnmass = [0] * len(data['yarnmass'])
            fconc_sta = [0] * len(data['yarnmass'])
            fconc_mid = [0] * len(data['yarnmass'])
            fconc_end = [0] * len(data['yarnmass'])
            for ind in range(len(data['yarnmass'])):
                yarnmass[ind] = np.empty(0, float)
                fconc_sta[ind] =  np.empty(0, float)
                fconc_mid[ind] =  np.empty(0, float)
                fconc_end[ind] =  np.empty(0, float)
        for ind, dt in enumerate(data['yarnmass']):
            yarnmass[ind] = np.append(yarnmass[ind], dt[::EVERY])  #self.yarnmass
            fconc_sta[ind] = np.append(fconc_sta[ind], data['fconc_sta'][ind][::EVERY])
            fconc_mid[ind] = np.append(fconc_mid[ind], data['fconc_mid'][ind][::EVERY])
            fconc_end[ind] = np.append(fconc_end[ind], data['fconc_end'][ind][::EVERY])
        totyarnmass = np.append(totyarnmass, data['totyarnmass'][::EVERY]) #self.totyarnmass
        totroommass = np.append(totroommass, data['totroommass'][::EVERY]) #self.totroommass
        
  #times are up to max time, not up to end of simulation, so trim
  ltimes += [times]
  lsol += [sol]
  ltotyarnmass += [totyarnmass]
  ltotroommass += [totroommass]
  lyarnmass += [yarnmass]
  lfconc_sta += [fconc_sta]
  lfconc_mid += [fconc_mid]
  lfconc_end += [fconc_end]
  prob += 1
#now we determine how to interpolate over the grid for x0
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
colors = 'bgrkmy'
lencolors = len(colors)
color = 0
savedir = {}

def view_sol(times, sol, label=None):
    global colors, lencolors, color
    ind = 0
    usecolor = colors[color%lencolors]
    color += 1
    for dataind, interpdat in enumerate(plotdata):
        xval = x0[dataind]
        cellstart, interpval = interpdat
        conc_in_point = (interpval * sol[:, dataind] 
                            + (1-interpval) * sol[:, dataind+1] )
        print ('conc in end point', conc_in_point[-1])
        plt.rc("font", family="serif")
        plt.rc("font", size=10)
        width = 4.5  #width in inches
        height = 1.4 #height in inches
        plt.rc("figure.subplot", left=(50/72.27)/width)
        plt.rc("figure.subplot", right=(width-10/72.27)/width)
        plt.rc("figure.subplot", bottom=(14/72.27)/height)
        plt.rc("figure.subplot", top=(height-7/72.27)/height)
        fig = plt.figure(ind)
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
        #plt.gca().yaxis.set_major_formatter(pylab.FormatStrFormatter('%e'))
        plt.title('Concentration at position %g mm' % xval)
        plt.plot(times, conc_in_point, label=label)
        #plt.ylim(0, maxv*1.1)
        plt.plot(times, np.ones(len(times)) * saturation_conc, 'k--')
        plt.plot(times, np.ones(len(times)) * treshold, 'k--')
        savedir['times'] = times
        savedir['concpos_%g'%xval] = conc_in_point
        if SHOWLEGEND:
            plt.legend()
            
        if GRAY:
            setFigLinesBW(fig)
        plt.draw()
        ind += 1
        #same plot in unit minutes !
        plt.rc("font", family="serif")
        plt.rc("font", size=10)
        width = 4.5  #width in inches
        height = 1.4 #height in inches
        plt.rc("figure.subplot", left=(50/72.27)/width)
        plt.rc("figure.subplot", right=(width-10/72.27)/width)
        plt.rc("figure.subplot", bottom=(14/72.27)/height)
        plt.rc("figure.subplot", top=(height-7/72.27)/height)
        fig = plt.figure(ind)
        plt.gca().set_xlabel('Time [min]')
        plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
        #plt.gca().yaxis.set_major_formatter(pylab.FormatStrFormatter('%e'))
        plt.title('Concentration at position %g mm' % xval)
        plt.plot(times/60, conc_in_point, usecolor, label=label)
            #if double:
            #plt.plot(times/60, 2*conc_in_point, usecolor+'--', label='2 x '+label)
        #plt.ylim(0,  treshold*1.1)
        ##plt.plot(times/60, np.ones(len(times)) * saturation_conc, 'k--')
        plt.plot(times/60, np.ones(len(times)) * treshold, 'k-')
        if SHOWLEGEND:
            plt.legend()
        if GRAY:
            setFigLinesBW(fig)
        plt.draw()
        ind +=1
        #same plot in unit hour !
        plt.rc("font", family="serif")
        plt.rc("font", size=10)
        width = 4.5  #width in inches
        height = 1.4 #height in inches
        plt.rc("figure.subplot", left=(50/72.27)/width)
        plt.rc("figure.subplot", right=(width-10/72.27)/width)
        plt.rc("figure.subplot", bottom=(14/72.27)/height)
        plt.rc("figure.subplot", top=(height-7/72.27)/height)
        fig = plt.figure(ind)
        plt.gca().set_xlabel('Time [hour]')
        plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
        #plt.gca().yaxis.set_major_formatter(pylab.FormatStrFormatter('%e'))
        plt.title('Concentration at position %g mm' % xval)
        plt.plot(times/60/60, conc_in_point, label=label)
        plt.plot(times/60/60, np.ones(len(times)) * saturation_conc, 'k--')
        plt.plot(times/60/60, np.ones(len(times)) * treshold, 'k--')
        if SHOWLEGEND:
            plt.legend()
        if GRAY:
            setFigLinesBW(fig)
        plt.draw()
        ind +=1
        #plt.savefig('AIconc_%03.1f_mm' % xval + '.png')
    return ind

def view_fiberconc(ind, fibertimes, yfiberconc_center, yfiberconc_middle,
                   yfiberconc_surf, label=None):
    fignr = ind
    plt.ion()
    fibnr = 0
    for fiberconc_center, fiberconc_middle, fiberconc_surf in zip(yfiberconc_center, yfiberconc_middle, yfiberconc_surf):
        fig = plt.figure(fignr)
        fignr += 1
        #center of fiber
        plt.rc("font", family="serif")
        plt.rc("font", size=10)
        #width = 4.5  #width in inches
        #height = 1.4 #height in inches
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Conc [$\mu$g/mm$^3$]')
        plt.title('Conc AI in fiber center')
        plt.plot(fibertimes, fiberconc_center, '-', label=label)
        if GRAY:
            setFigLinesBW(fig)
        if SHOWLEGEND:
            plt.legend()

        fig = plt.figure(fignr)
        fignr += 1
        #middle of fiber    
        plt.rc("font", family="serif")
        plt.rc("font", size=10)
        #width = 4.5  #width in inches
        #height = 1.4 #height in inches
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Conc [$\mu$g/mm$^3$]')
        plt.title('Conc AI in middle of fiber coating')
        plt.plot(fibertimes, fiberconc_middle, '-', label=label)
        if GRAY:
            setFigLinesBW(fig)
        if SHOWLEGEND:
            plt.legend()

        fig = plt.figure(fignr)
        fignr += 1
        #surface of fiber
        plt.rc("font", family="serif")
        plt.rc("font", size=10)
        #width = 4.5  #width in inches
        #height = 1.4 #height in inches
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Conc [$\mu$g/mm$^3$]')
        plt.title('Conc AI on fiber surface')
        plt.plot(fibertimes, fiberconc_surf, '-', label=label)
        if GRAY:
            setFigLinesBW(fig)
        if SHOWLEGEND:
            plt.legend()  
        savedir['fibertimes'] = fibertimes
        savedir['fiberconc_surf_%d' % fibnr] = fiberconc_surf
        fibnr += 1

def view_sol_mass(ind, times, yarnmass, totyarnmass, totroommass, label=''):
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
        fig = plt.figure(fignr)
        fignr += 1
        plt.gca().set_xlabel('Time [s]')
        plt.gca().set_ylabel('Mass [$\mu$g]')
        plt.title('Mass AC in yarn type %d' % ind)
        plt.plot(times, ymass, label=(label or '')+" yarn %d" % ind)
        if SHOWLEGEND:
            plt.legend()
        if GRAY:
            setFigLinesBW(fig)

    fig = plt.figure(fignr)
    fignr += 1
    plt.gca().set_xlabel('Time [h]')
    plt.gca().set_ylabel('Mass [$\mu$g]')
    plt.title('Mass AC in the textile')
    #plt.plot([0,],[28.935,], 'r*')
    plt.plot(times/60/60, totyarnmass, label=label)
    if GRAY:
        setFigLinesBW(fig)
    if SHOWLEGEND:
        plt.legend()
    fig = plt.figure(fignr)
    fignr += 1
    plt.gca().set_xlabel('Time [s]')
    plt.gca().set_ylabel('Mass [$\mu$g]')
    plt.title('Mass AC in the room')
    plt.plot(times, totroommass, label=label)
    if GRAY:
        setFigLinesBW(fig)
    if SHOWLEGEND:
        plt.legend()
    #plot to check mass conservation
    fig = plt.figure(fignr)
    fignr += 1
    plt.gca().set_xlabel('Time [s]')
    plt.gca().set_ylabel('Mass [$\mu$g]')
    plt.title('Total Mass AC')
    plt.plot(times, totroommass+totyarnmass, label=label)
    if GRAY:
        setFigLinesBW(fig)
    if SHOWLEGEND:
        plt.legend()
    savedir['totyarnmass'] = totyarnmass
    savedir['totroommass'] = totroommass
    return fignr

# we plot the data
plt.ion()
if len(lsol) == len(LABELS):
    pass
else:
    print (len(lsol), len(LABELS))
    LABELS = [None] * len(sol)

for times, sol, label in zip(ltimes, lsol, LABELS):
    ind = view_sol(times, sol, label)
    
for times, yarnmass, totyarnmass, totroommass, label in zip(ltimes, lyarnmass,
                                        ltotyarnmass, ltotroommass, LABELS):
    newind = view_sol_mass(ind, times, yarnmass, totyarnmass, totroommass, label)

for times, fconc_sta, fconc_mid, fconc_end, label in zip(ltimes, lfconc_sta, 
                                             lfconc_mid, lfconc_end, LABELS):
    view_fiberconc(newind, times, fconc_sta, fconc_mid, fconc_end, label)

#store an npz file with main data for multiple printing
if not PROBS:
    np.savez(BASEDIR + PROBTOLOAD + os.sep + 'figout.npz',
              **savedir)

plt.show(block=True) #block=True)
