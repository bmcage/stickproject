# Copyright (C) 2009  Pavol Kison
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
    Utilities
"""

#-------------------------------------------------------------------------
#
# python modules
#
#-------------------------------------------------------------------------
from __future__ import division
import os
from subprocess import Popen, PIPE
import numpy as np
AXES3D = False
try:
    from mpl_toolkits.mplot3d import Axes3D
    AXES3D = True
except:
    pass
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------
#
# Stick modules
#
#-------------------------------------------------------------------------
import stick.const as const

#-------------------------------------------------------------------------
#
# Utility functions
#
#-------------------------------------------------------------------------

OUTPUTDIR = ''

def set_outputdir(odir):
    global OUTPUTDIR
    OUTPUTDIR = odir

def write_info(elapsed_time, cpu_time):
    filename=os.path.join(OUTPUTDIR, 'info')
    outfile=open(filename,'w')
    outfile.write('\nTime to calculate solution\n')
    outfile.write('    elapsed time : %g\n' % elapsed_time )
    outfile.write('    cpu time     : %g\n' % cpu_time )
    outfile.close()
    
def write_sol(scipy_array, time = ''):
    print "writing solution at time %s to disk"  % str(time)
    filename=os.path.join(OUTPUTDIR, 'sol%s.npy' % str(time))
    #text save,  ## we do quick binary instead
    #from scipy.io import write_array
    #write_array(filename, positions)
    np.save(filename, scipy_array)

def plot_posfig(positions, time='', show=False):
    """
    Plot a figure as given by position, a 3xN array. If show True, an image is
    shown, if False, the picture is written out at sol<time>.png
    """
    if not AXES3D:
        print 'utils.py: Cannot plot 3D images, toolkit not installed'
        return
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_aspect('equal')
    ax.set_autoscale_on(False)

    ax.plot(positions[0, :], positions[1, :], positions[2, :], 'o-',
            label='parametric curve')
    if show:
        plt.show()
    else:
        plt.savefig(os.path.join(OUTPUTDIR, 'sol%s.png' % str(time)))

def read_sol(filepath):
    """
    reads a npy file in, assuming a format to know the time
    Returns: time, array
    """
    filename = os.path.split(filepath)[1]
    try:
        time = float(filename[3:-4])
        data = np.load(filepath)
    except NotImplementedError:
        return None, None
    return time, data

def open_file_with_default_application( file_path ):
    """
    Launch a program to open an arbitrary file. The file will be opened using 
    whatever program is configured on the host as the default program for that 
    type of file.
    """
    
    norm_path = os.path.normpath( file_path )
    
    if not os.path.exists(norm_path):
        print "%s does not exist" % file_path
        return
        
    if os.sys.platform == 'win32':
        try:
            os.startfile(norm_path)
        except WindowsError, msg:
            print "Error Opening File. " + str(msg)
    else:
        search = os.environ['PATH'].split(':')
        for path in search:
            prog = os.path.join(path, 'xdg-open')
            if os.path.isfile(prog):
                os.spawnvpe(os.P_NOWAIT, prog, [prog, norm_path], os.environ)
                return

def add_root_to_file(file):
    """
    If file is not an absolute path, add the root dir
    """
    if os.path.isabs(file):
        return file
    else:
        return const.ROOT_DIR + os.sep + file