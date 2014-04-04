#
# Copyright (C) 2010  B. Malengier
# Copyright (C) 2010  P.Li
# Copyright (C) 2014  T. Goessens
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
    plotting time vs log conc
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

BASEDIR = '/Users/Tine/stickproject/'
#PROBS = False  #set tot False if only one problem
PROBTOLOAD = 'fabric.ini'
#PROBTOLOAD = 'fabricbednetY335_Deet.ini_50nmol8hour'
#all problems must be over the same grid !
#PROBSTOLOAD = ['fabricbednetY335_Deet.ini_50nmol8hour',
#              'fabricbednetY335_Deet.ini_100nmol8hour',
#              'fabricbednetY335_Deet.ini_200nmol8hour',
#              'fabricbednetY335_Deet.ini_400nmol8hour']
#LABELS = ['50 nmol', '100 nmol', '200 nmol', '400 nmol']
#PROBSTOLOAD = ['fabricmuslin_Deet.ini_25nmol2min_b',
#    'fabricmuslin_Deet.ini_20nmol2min_b', 'fabricmuslin_Deet.ini_15nmol2min_b',
#    'fabricmuslin_Deet.ini_10nmol2min_b', 'fabricmuslin_Deet.ini_05nmol2min_b']
##LABELS = ['25 nmol', '20 nmol', '15 nmol', '10 nmol', '5 nmol']

ARGS = ['/fiberconccenter.txt',
        '/fiberconcmiddle.txt',
        '/fiberconcsurface.txt',
        '/yarnconccenter.txt',
        '/yarnconcsurface.txt',
        '/roomconcLEFT.txt',
        '/roomconcMIDDLE.txt',
        '/roomconcRIGHT.txt'
        ]
        
#INDEX = range(1) #range(4) # what dumped data to load
#EVERY = 60 #1    # what time data to skip to reduce plotting time
#if not PROBS:
#   PROBSTOLOAD = [PROBTOLOAD]

#prob = 0
j=0

for ARG in ARGS:
    n = len(ARGS)
    concentration = open(BASEDIR + PROBTOLOAD + ARG,'r')
    conc = concentration.readlines()
    globals()['times%s' % j] = np.empty(len(conc),float)
    globals()['sol%s' % j] = np.empty(len(conc),float)
    globals()['logconc%s' % j] = np.empty(len(conc),float)
    for i in range(len(conc)):
        conc[i] = conc[i].strip()
        con = conc[i].split()
        globals()['times%s' % j][i] = con[0]
        globals()['sol%s' % j][i] = con[1]
        if float(con[1])>0:
            globals()['logconc%s' % j][i] = math.log(float(con[1]))
        else:
            globals()['logconc%s' % j][i] = 0
    globals()['times%s' % j]= np.sort(globals()['times%s' % j])
    globals()['sol%s' % j]= np.sort(globals()['sol%s' % j])
    globals()['logconc%s' % j]= np.sort(globals()['logconc%s' % j])
    j+=1

plotdata = []
colors = 'bgrkmy'
lencolors = len(colors)
color = 0

times=[]
logconc = []
for j in range(len(ARGS)):
    times.append(globals()['times%s' %j])
    logconc.append(globals()['logconc%s' %j])

def view_sol(time, concentration):
    global colors, lencolors, color
    usecolor = colors[color%lencolors]
    color += 1
    plt.rc("font", family="serif")
    plt.rc("font", size=10)
    plt.figure(figsize=(8, 5))
    plt.gca().set_xlabel('Time [s]')
    plt.gca().set_ylabel('Concentration [$\mu$g/mm$^3$]')
    plt.title('Log(concentration) vs time')
    plt.plot(time[0], concentration[0],'r-',time[1], concentration[1],'g-',time[2], concentration[2],'b-', time[3], concentration[3], 'c-', time[4], concentration[4],'y-', time[5], concentration[5],'m-',time[6],concentration[6],'k-',time[7], concentration[7],'r--')
    plt.legend( ["fiberconccenter.txt",
                      "fiberconcmiddle.txt",
                      "fiberconcsurface.txt",
                      "yarnconccenter.txt",
                      "yarnconcsurface.txt",
                      "roomconcLEFT.txt",
                      "roomconcMIDDLE.txt",
                      "roomconcRIGHT.txt"
                 ],loc='upper center', bbox_to_anchor=(0.5, 0.85),
               ncol=3, fancybox=True, shadow=True
)
    plt.draw()

   # we plot the data
plt.ion()
view_sol(times, logconc)

#if len(lsol) == len(LABELS):
# pass
#else:
# print (len(lsol), len(LABELS))
#LABELS = [None] * len(sol)

#for times, sol, label in zip(ltimes, lsol, LABELS):
#ind = view_sol(times, sol, label)

#for times, yarnmass, totyarnmass, totroommass, label in zip(ltimes, lyarnmass,
#ltotyarnmass, ltotroommass, LABELS):
# view_sol_mass(ind, times, yarnmass, totyarnmass, totroommass, label)

plt.show(block=True)



