#! /usr/bin/env python

# Copyright (C) 2010  P. Li, B. Malengier
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
    Main program that reads in ini file for yarn 2D simulation and decides 
    how to handle it
"""
#-------------------------------------------------------------------------
#
# Python modules
#
#-------------------------------------------------------------------------
from __future__ import division
import sys, os, shutil
import getopt
import time
import numpy as N
import scipy as S

#-------------------------------------------------------------------------
#
# local modules
#
#-------------------------------------------------------------------------
from lib.utils.utils import set_outputdir
import const
import yarn2d.config as conf

#-------------------------------------------------------------------------
#
# Initialization
#
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
#
# the main prog
#
#-------------------------------------------------------------------------

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        options, leftargs = getopt.getopt(argv[1:],
                                          conf.SHORTOPTS, conf.LONGOPTS)
    except getopt.GetoptError, msg:
        print msg
        # return without filling anything if we could not parse the args
        print "Error parsing the arguments: %s " % argv[1:]
        sys.exit(0)
    if leftargs:
        print 'yarn2d.py does not understand argument %s' % leftargs
        sys.exit(0)

    inifile = conf.INIFILE_DEFAULT
    outputdir = const.DATA_DIR
    onlymesh = False
    writeini = False
    for opt_ix in range(len(options)):
        option, value = options[opt_ix]
        if option in ( '-i', '--inifile'):
            inifile = value
        elif option in ('-o', '--outputdir'):
            outputdir = value
        elif option in ('--mesh'):
            onlymesh = True
        elif option in ('--write-ini'):
            writeini = True
    
    #Parse ini file to obtain parameters.
    cfg = conf.Yarn2dConfigManager.get_instance(inifile)

    #create outputdir if not existing
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    #create outputdir for this run, remove if existing
    outputdir = outputdir + os.sep + os.path.basename(inifile)    
    #determine whether we make a new file for yarn2d or read the old file; the
    #default value is 'False', which means the new file will be generated when
    #the programme runs
    read_old_file = cfg.get('general.read')
    if not read_old_file:
        if not os.path.isdir(outputdir):
            os.mkdir(outputdir)
        else:
            shutil.rmtree(outputdir)
            os.mkdir(outputdir)
        set_outputdir(outputdir)
        #store the ini file in the outputdir so the experiment can be repeated
        shutil.copy(inifile, outputdir)
    else:
        set_outputdir(outputdir)
        print "ready to read the old file"
        
    if writeini:
        print "Writing out ini file cleaned.ini to outputdir %s" % outputdir
        cfg.save(outputdir + os.sep + 'cleaned.ini')
        sys.exit()
    if not hasattr(cfg, 'onlymesh'):
        cfg.onlymesh = onlymesh

    #create the correct model, and run it
    from yarn2d.yarn2dmodel import Yarn2DModel
    model = Yarn2DModel(cfg)

    #pass further execution to the mode
    model.run()

if __name__ == "__main__":
    main()