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
Experiments 201106  "
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
LE_2ndPEAK = 2 #mL used for second peak
LE_1stPEAKextra = 20 #mL used for first peak for extra exp
LE_2ndPEAKextra = 20 #mL used for second peak for extra exp
EXTRA = ['35','36']  #samples that use LE extra
PERIOD = ['march2011', 'midjune2011', 'endjune2011']
PERIOD2DATE = { #begin to end date (weeknumber, year)
    'march2011': [(9,2011), (11,2011)],
    'midjune2011': [(25,2011), (25,2011)], #15 weeks later
    'endjune2011': [(26,2011), (26,2011)],
    }
LE_HEADERS = ('Compound', 'Sample', 'Weight g', 'Peak 1st LE', 'Peak 2nd LE', 'period')
LE_1  = (   (DEET, 9,  0.0436, 0, 0, 'midjune2011'),
            (DEET, 10, 0.0564, 0, 0, 'midjune2011'),
            (DEET, 18, 0.0515, 0, 0, 'midjune2011'),
        )
LE_A5 = (   (DEET, 4, 0.045, 39304700, 2551721, 'march2011'),
            (DEET, 5, 0.048, 35528900, 3361885, 'march2011'),
            (DEET, 6, 0.048, 32224500, 2649236, 'march2011'),
            (PER, 4, 0.045, 4496200, 619155, 'march2011'),
            (PER, 5, 0.048, 6477600, 367440, 'march2011'),
            (PER, 6, 0.048, 3203100, 606709, 'march2011'),
            (DEET, 19, 0.0476, 290499, 871726, 'endjune2011'),
            (DEET, 20, 0.0529, 297978, 692612, 'endjune2011'),
            (DEET, 21, 0.0515, 272392, 840634, 'endjune2011'),
            (DEET, 22, 0.0554, 279433, 856935, 'endjune2011'),
            (PER, 19, 0.0476, 11528690, 362629, 'endjune2011'),
            (PER, 20, 0.0529, 12301786, 277688, 'endjune2011'),
            (PER, 21, 0.0515, 11689620, 365262, 'endjune2011'),
            (PER, 22, 0.0554, 13294682, 326642, 'endjune2011'),
            (DEET, 35, 0.0447, 32982, 0, 'endjune2011'),
            (PER, 35, 0.0447, 1267861, 3550, 'endjune2011'),
        )
LE_C5 = (   (DEET, 1, 0.063, 1446700, 357011, 'march2011'),
            (DEET, 2, 0.059, 3437700, 643069, 'march2011'),
            (DEET, 3, 0.055, 3243400, 590819, 'march2011'),
            (PER, 1, 0.063, 8854400, 1193029, 'march2011'),
            (PER, 2, 0.059, 7332300, 974129, 'march2011'),
            (PER, 3, 0.055, 7305000, 823552, 'march2011'),
            (DEET, 23, 0.0528, 1466525, 32925, 'endjune2011'),
            (DEET, 24, 0.061, 1438981, 41870, 'endjune2011'),
            (DEET, 25, 0.0541, 646714, 12923, 'endjune2011'),
            (DEET, 26, 0.047, 578745, 17730, 'endjune2011'),
            (PER, 23, 0.0528, 12329718, 169823, 'endjune2011'),
            (PER, 24, 0.061, 12376012, 210634, 'endjune2011'),
            (PER, 25, 0.0541, 12541304, 121612, 'endjune2011'),
            (PER, 26, 0.047, 12656640, 170574, 'endjune2011'),
        )
LE_C3 = (   (DEET, 7, 0.061, 24436400,1857128, 'march2011'),
            (DEET, 8, 0.061, 24059000,2568056, 'march2011'),
            (DEET, 9, 0.062, 36600600,21992888, 'march2011'),
            (DEET, 15, 0.0457, 828062,62087, 'midjune2011'),
            (DEET, 16, 0.0496, 792785,52594, 'midjune2011'),
            (DEET, 17, 0.0456, 387493,1793852, 'endjune2011'),
            (PER, 15, 0.0457, 131076,169592, 'midjune2011'),
            (PER, 16, 0.0496, 132592,0, 'midjune2011'),
            (PER, 17, 0.0456, 109862,0, 'endjune2011'),
        )
LE_A3 = (   (DEET, 10, 0.051, 36613200,2896294, 'march2011'),
            (DEET, 11, 0.049, 26140500, 2143200, 'march2011'),
            (DEET, 12, 0.068, 62181000, 2896294, 'march2011'),
            (DEET, 11, 0.045, 796834, 53678, 'midjune2011'),
            (DEET, 12, 0.0494, 810956, 51353, 'midjune2011'),
            (DEET, 13, 0.0439, 603561, 45985, 'midjune2011'),
            (DEET, 14, 0.0476, 527676, 45453, 'midjune2011'),
            (PER, 11, 0.045, 465270, 11367, 'midjune2011'),
            (PER, 11, 0.0494, 909316, 10470, 'midjune2011'),
            (PER, 11, 0.0439, 303554, 4118, 'midjune2011'),
            (PER, 11, 0.0476, 317038, 4352, 'midjune2011'),
        )
LE_C2 = (   (DEET, 16, 0.053, 9957200,791148, 'march2011'),
            (DEET, 17, 0.064, 20680900,1806349, 'march2011'),
            (DEET, 18, 0.072, 13272600,1018263, 'march2011'),
            (PER, 16, 0.053, 9037400,1234738, 'march2011'),
            (PER, 17, 0.064, 14296600,1663628, 'march2011'),
            (PER, 18, 0.072, 14850600,1928021, 'march2011'),
            (PER, 5, 0.0498, 19794606,540734, 'midjune2011'),
            (PER, 6, 0.0544, 15702738,543919, 'endjune2011'),
            (PER, 7, 0.0544, 15077898,334366, 'endjune2011'),
            (PER, 8, 0.0447, 12612794,259176, 'endjune2011'),
        )
LE_A2 = (   (DEET, 19, 0.042, 8676100,502327, 'march2011'),
            (DEET, 20, 0.049, 5261600,260656, 'march2011'),
            (DEET, 21, 0.057, 2133600,284562, 'march2011'),
            (PER, 19, 0.042, 5273700,436837, 'march2011'),
            (PER, 20, 0.049, 4885400,479561, 'march2011'),
            (PER, 21, 0.057, 3656000,756018, 'march2011'),
            (PER, 1, 0.0566, 18905485, 469107, 'midjune2011'),
            (PER, 2, 0.0564, 17929555, 569509, 'endjune2011'),
            (PER, 3, 0.0515, 21990415, 651812, 'endjune2011'),
            (PER, 4, 0.0455, 20142143, 615055, 'endjune2011'),
        )
LE_A6 = (   (DEET, 27, 0.0569, 1670529,69932, 'endjune2011'),
            (DEET, 28, 0.0484, 1456575,69502, 'endjune2011'),
            (DEET, 29, 0.0488, 1777403,67690, 'endjune2011'),
            (DEET, 30, 0.0564, 1818404,51930, 'endjune2011'),
            (PER, 27, 0.0569, 20377690,587595, 'endjune2011'),
            (PER, 28, 0.0484, 17634452,462370, 'endjune2011'),
            (PER, 29, 0.0488, 26968453,742955, 'endjune2011'),
            (PER, 30, 0.0564, 25484190,479484, 'endjune2011'),
        )
LE_C6 = (   (DEET, 22, 0.049, 109706800,3305107, 'march2011'),
            (DEET, 23, 0.058, 193078400,11983694, 'march2011'),
            (DEET, 24, 0.05, 155041100,6057975, 'march2011'),
            (PER, 22, 0.049, 21092900,839111, 'march2011'),
            (PER, 23, 0.058, 6047500,1058782, 'march2011'),
            (PER, 24, 0.05, 12701100,573436, 'march2011'),
            (DEET, 31, 0.0472, 2791424,110614, 'endjune2011'),
            (DEET, 32, 0.0517, 2839494,151767, 'endjune2011'),
            (DEET, 33, 0.0487, 2800502,99872, 'endjune2011'),
            (DEET, 34, 0.0478, 2620122,119003, 'endjune2011'),
            (PER, 31, 0.0472, 10213386,258442, 'endjune2011'),
            (PER, 32, 0.0517, 10011091,411195, 'endjune2011'),
            (PER, 33, 0.0487, 8962571,149877, 'endjune2011'),
            (PER, 34, 0.0478, 8089477,213776, 'endjune2011'),
            (DEET, 36, 0.0506, 103010,0, 'endjune2011'),
            (PER, 36, 0.0506, 790418,0, 'endjune2011'),
        )
SAMPLES = {
    '1' : (LE_1,  'Untreated sample'),
    'A5': (LE_A5, 'A5 - DEET+Per, Foulard'),
    'A3': (LE_A3, 'A3 - DEET,     Foulard'),
    'A2': (LE_A2, 'A2 - Per,      Foulard'),
    'A6': (LE_A6, 'A6 - Per+DEET, Foulard'),
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

def calib_DEET_may2011(x):
    """ Conversion of peak area to microg/mL"""
    return x / 1416.6 / 1000.

def calib_DEET_midjune2011(x):
    """ Conversion of peak area to microg/mL"""
    return x / 1139.8 / 1000.

def calib_DEET_endjune2011(x):
    """ Conversion of peak area to microg/mL"""
    return x / 811.17 / 1000.

def calib_PER_may2011(x):
    """ Conversion of peak area to microg/mL"""
    return x / 221.06 / 1000.

def calib_PER_midjune2011(x):
    """ Conversion of peak area to microg/mL"""
    return x / 533.7 / 1000.

def calib_PER_endjune2011(x):
    """ Conversion of peak area to microg/mL"""
    return x / 404.25 / 1000.

def period2calib(period, compound):
    if period == 'march2011':
        if compound == DEET:
            return calib_DEET_may2011
        elif compound == PER:
            return calib_PER_may2011
    elif period == 'midjune2011':
        if compound == DEET:
            return calib_DEET_midjune2011
        elif compound == PER:
            return calib_PER_midjune2011
    elif period == 'endjune2011':
        if compound == DEET:
            return calib_DEET_endjune2011
        elif compound == PER:
            return calib_PER_endjune2011

def comp_weight(xlst, type=None, time=None):
    """ compute ng/mL of the peaks
        xlst should be the different peak areas measured"""
    weight = period2calib(time, type)(numpy.array(xlst))
    return weight

def comp_content_per_gramtextile(ledata):
    """ compute microg of compount per gram of material, based on 
        given Liquid Extraction data ledata
        Returns compound in microg/g and efficiency first peak
    """
    peaks = numpy.array([ledata[3], ledata[4]], float)
    weight = comp_weight(peaks, ledata[0], ledata[5])
    if str(ledata[1]) in EXTRA:
        weight[0] *= LE_1stPEAKextra
        weight[1] *= LE_2ndPEAKextra
    else:        
        weight[0] *= LE_1stPEAK
        weight[1] *= LE_2ndPEAK
    total_amount = weight[0] + weight[1]
    if total_amount == 0.:
        return 0., 0.
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
        avgDEET = {};
        nrDEET = {}
        avgPER = {};
        nrPER = {}
        for ledata in data[0]:
            result = comp_content_per_gramtextile(ledata)
            if detail:
                print '  Sample %02d: eff 1stLE %f; content %s =  %f microg/g'\
                      ' Period %s' % (
                        ledata[1], result[1], ledata[0], result[0], ledata[5])
            if ledata[0] == DEET:
                if ledata[5] not in avgDEET:
                    avgDEET[ledata[5]] = 0.
                    nrDEET[ledata[5]] = 0
                avgDEET[ledata[5]] += result[0]
                nrDEET[ledata[5]]  += 1
            elif ledata[0] == PER:
                if ledata[5] not in avgPER:
                    avgPER[ledata[5]] = 0.
                    nrPER[ledata[5]] = 0
                avgPER[ledata[5]]  += result[0]
                nrPER[ledata[5]]  += 1
        for key in nrDEET:
            print ' Avg DEET= %f microg/g in period %s' % \
                    (avgDEET[key]/nrDEET[key], key)
        for key in nrPER:
            print ' Avg PERM= %f microg/g in period %s' % \
                    (avgPER[key]/nrPER[key], key)
        if ('midjune2011' in nrDEET) and ('endjune2011' in nrDEET):
            print ' Avg DEET= %f microg/g in June' % \
                    ((avgDEET['midjune2011']+avgDEET['endjune2011'])/ \
                     (nrDEET['midjune2011']+nrDEET['endjune2011']))
        if ('midjune2011' in nrPER) and ('endjune2011' in nrPER):
            print ' Avg PER= %f microg/g in June' % \
                    ((avgPER['midjune2011']+avgPER['endjune2011'])/ \
                     (nrPER['midjune2011']+nrPER['endjune2011']))

if __name__ == '__main__': 
    main(detail=False)