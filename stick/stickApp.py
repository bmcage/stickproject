#!/usr/bin/env python

# Copyright (C) 2009  B. Malengier
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
    Main program that decides what part of STICK to run
"""
#-------------------------------------------------------------------------
#
# Python modules
#
#-------------------------------------------------------------------------
from __future__ import division
import sys

#-------------------------------------------------------------------------
#
# local modules
#
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#
# Initialization
#
#-------------------------------------------------------------------------

ARGS = sys.argv
# module with the program, main() will be executed
PROGS = {
    'fiber1d': 'stick.fiber1d.fiberprog',
    'yarn1d': 'stick.yarn1d.yarn1dprog',
    'yarn2d': 'stick.yarn2d.yarn2dprog',
    'bednet': 'stick.bednet.bednetprog',
    'pcm': 'stick.pcm.pcmprog',
    'room': 'stick.room.roomprog',
    'fiberfabric': 'stick.fiberfabric.fiberfabricprog',
    }

#-------------------------------------------------------------------------
#
# the main prog
#
#-------------------------------------------------------------------------

def main(prog=None, arg=None):
    if prog is None:
        if len(ARGS) < 2 or not ARGS[1] in PROGS:
            print "ERROR: first argument must be the program to run."
            print " possibilities: " + ' '.join(PROGS)
            sys.exit()
        remove = None
        for ind, val in enumerate(ARGS):
            if val == '--backend':
                backend = ARGS[ind+1]
                if backend in ['ps', 'PS', 'eps']:
                    import matplotlib
                    print 'pyplot backend', matplotlib.get_backend()
                    matplotlib.use('PS') # 'cairo.pdf', 'cairo.png'
                    import const
                    const.FIGFILEEXT = '.eps'
                    print 'changed to backend', matplotlib.get_backend()
                remove = ind
                break
        if remove:
            del ARGS[remove+1]
            del ARGS[remove]
        modname = PROGS[ARGS[1]]
    else:
        try:
            modname = PROGS[prog]
        except:
            print "ERROR: prog must be the program to run."
            print " possibilities: " + ' '.join(PROGS)
            sys.exit()
    _topmodule = __import__(modname)
    _realmodule = sys.modules[modname]
    if prog is None:
        #run the program, shift ARGS by two (program and progtype
        _realmodule.main(ARGS[2:])
    else:
        _realmodule.main(arg)

if __name__ == "__main__":
    main()
