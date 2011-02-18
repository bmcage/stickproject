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
    'fiberode': 'fiberode.fiberprog',
    'fiberfipy': 'fiberfipy.fiberprog',
    'yarn1d': 'yarn1d.yarn1dprog',
    }

#-------------------------------------------------------------------------
#
# the main prog
#
#-------------------------------------------------------------------------

def main():
    if not ARGS[1] in PROGS:
        print "ERROR: first argument must be the program to run."
        print " possibilities: " + ' '.join(PROGS)
        sys.exit()
    
    modname = PROGS[ARGS[1]]
    _topmodule = __import__(modname)
    #run the program, shift ARGS by one
    _realmodule = sys.modules[modname]
    _realmodule.main(ARGS[1:])

if __name__ == "__main__":
    main()

