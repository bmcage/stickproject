#
# Copyright (C) 2009  Benny Malengier
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

"""
Provides constants for other modules

To access:
import const
const.VALUE
"""

#-------------------------------------------------------------------------
#
# Standard python modules
#
#-------------------------------------------------------------------------
import os
import platform

#-------------------------------------------------------------------------
#
# Website
#
#-------------------------------------------------------------------------
URL_HOMEPAGE    = "http://gitorious.org/stickproject"

#-------------------------------------------------------------------------
#
# paths
#
#-------------------------------------------------------------------------
USER_HOME = os.path.expanduser('~') 
DATA_DIR  = USER_HOME + os.sep + 'stickproject'
#edit following to point to root of paths in ini files
ROOT_DIR = os.path.expanduser('~') + os.sep + 'git' + os.sep + 'diffusionit'\
                                   + os.sep + 'data'

# dirs that need to be created
USER_DIRLIST = (DATA_DIR,)

#-------------------------------------------------------------------------
#
# About box information
#
#-------------------------------------------------------------------------
PROGRAM_NAME   = "STICK"
VERSION        = "0.0.1a"
COPYRIGHT_MSG  = u"\u00A9 2010 Tinne Goessens\n" \
                 u"\u00A9 2010 Pei Li\n" \
                 u"\u00A9 2006-2010 Benny Malengier\n"
COMMENTS       = "Sophisticated Textile Information Computing Kit"
PYTHONVERSION  = platform.python_version_tuple()

#-------------------------------------------------------------------------
#
# Constants
#
#-------------------------------------------------------------------------
