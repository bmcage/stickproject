#
# Copyright (C) 2011       Benny Malengier <bm@cage.ugent.be>
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
#

from __future__ import division
"""
This package contains experimental data, usefull for plotting.

Experiments 201104 from Department Organic Chemistry
"""

#---------------------------------------------------------------
#
# System imports
#
#---------------------------------------------------------------
import os
import numpy
import time

#---------------------------------------------------------------
#
# local imports
#
#---------------------------------------------------------------


#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------

LE = 'Liquid extraction'
PER = 'Permethrine'
DEET = 'DEET'
LE_1stPEAK = 3 #mL used for first peak
LE_2ndPEAK = 2 #mL used for first peak
LE_HEADERS = ('Compound', 'Sample', 'Weight g', 'Peak 1st LE', 'Peak 2nd LE')
LE_A5 = (   (DEET, 4, 0.045, 39304700, 2551721),
            (DEET, 5, 0.048, 35528900, 3361885),
            (DEET, 6, 0.048, 32224500, 2649236),
            (PER, 4, 0.045, 4496200, 619155),
            (PER, 5, 0.048, 6477600, 367440),
            (PER, 6, 0.048, 3203100, 606709)
        )
LE_C5 = (   (DEET, 1, 0.063, 1446700, 357011),
            (DEET, 2, 0.059, 3437700, 643069),
            (DEET, 3, 0.055, 3243400, 590819),
            (PER, 1, 0.063, 8854400, 1193029),
            (PER, 2, 0.059, 7332300, 974129),
            (PER, 3, 0.055, 7305000, 823552)
        )
LE_C3 = (   (DEET, 7, 0.061, 24436400,1857128),
            (DEET, 8, 0.061, 24059000,2568056),
            (DEET, 9, 0.062, 36600600,21992888)
        )
LE_A3 = (   (DEET, 10, 0.051, 36613200,2896294),
            (DEET, 11, 0.049, 26140500, 2143200),
            (DEET, 12, 0.068, 62181000, 2896294),
        )
LE_C2 = (   (DEET, 16, 0.053, 9957200,791148),
            (DEET, 17, 0.064, 20680900,1806349),
            (DEET, 18, 0.072, 13272600,1018263),
            (PER, 16, 0.053, 9037400,1234738),
            (PER, 17, 0.064, 14296600,1663628),
            (PER, 18, 0.072, 14850600,1928021)
        )
LE_A2 = (   (DEET, 19, 0.042, 8676100,502327),
            (DEET, 20, 0.049, 5261600,260656),
            (DEET, 21, 0.057, 2133600,284562),
            (PER, 19, 0.042, 5273700,436837),
            (PER, 20, 0.049, 4885400,479561),
            (PER, 21, 0.057, 3656000,756018)
        )
LE_C6 = (   (DEET, 22, 0.049, 109706800,3305107),
            (DEET, 23, 0.058, 193078400,11983694),
            (DEET, 24, 0.05, 155041100,6057975),
            (PER, 22, 0.049, 21092900,839111),
            (PER, 23, 0.058, 6047500,1058782),
            (PER, 24, 0.05, 12701100,573436)
        )
SAMPLES = {
    'A5': (LE_A5, 'A5 - DEET+Per, Foulard'),
    'A3': (LE_A3, 'A3 - DEET,     Foulard'),
    'A2': (LE_A2, 'A2 - Per,      Foulard'),
    'C5': (LE_C5, 'C5 - DEET+Per, PBA'),
    'C3': (LE_C3, 'C3 - DEET,     PBA'),
    'C2': (LE_C2, 'C2 - Per,      PBA'),
    'C6': (LE_C6, 'C6 - Per+DEET, PBA'),
    }
    

#---------------------------------------------------------------
#
# Local functions
#
#---------------------------------------------------------------

def calib_DEET(x):
    """ Conversion of peak area to microg/mL"""
    return x / 1416.6 / 1000.

def calib_PER(x):
    """ Conversion of peak area to microg/mL"""
    return x / 221.06 / 1000.

def comp_weight(xlst, type=None):
    """ compute ng/mL of the peaks
        xlst should be the different peak areas measured"""
    if type == DEET:
        weight = calib_DEET(numpy.array(xlst))
    elif type == PER:
        weight = calib_PER(numpy.array(xlst))
    else:
        raise Exception, 'Unknown type'
    return weight

def comp_content_per_gramtextile(ledata):
    """ compute microg of compount per gram of material, based on 
        given Liquid Extraction data ledata
        Returns compound in microg/g and efficiency first peak
    """
    peaks = numpy.array([ledata[3], ledata[4]], float)
    weight = comp_weight(peaks, ledata[0])
    weight[0] *= LE_1stPEAK
    weight[1] *= LE_2ndPEAK
    total_amount = weight[0] + weight[1]
    eff = weight[0]/total_amount
    return total_amount/ledata[2], eff

#---------------------------------------------------------------
#
# Classes
#
#---------------------------------------------------------------

#---------------------------------------------------------------

def main(detail=False):
    samples = SAMPLES.keys()
    samples.sort()
    for key in samples:
        data = SAMPLES[key]
        print 'Values for sample: %s' % data[1]
        avgDEET = 0.;
        nrDEET = 0
        avgPER = 0.;
        nrPER = 0
        for ledata in data[0]:
            result = comp_content_per_gramtextile(ledata)
            if detail:
                print '  Sample %02d: eff 1stLE %f; content %s =  %f microg/g' % (
                        ledata[1], result[1], ledata[0], result[0])
            if ledata[0] == DEET:
                avgDEET += result[0]
                nrDEET += 1
            elif ledata[0] == PER:
                avgPER += result[0]
                nrPER += 1
        if nrDEET:
            print ' Avg DEET= %f microg/g' % (avgDEET/nrDEET)
        if nrPER:
            print ' Avg PERM= %f microg/g' % (avgPER/nrPER)

if __name__ == '__main__': 
    main(detail=False)