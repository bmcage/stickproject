#
# Copyright (C) 2012  B. Malengier
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
Starting Mayavi instance for the solution
"""
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------
from __future__ import division
import os
import os.path
import shutil
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time
import subprocess
import tempfile
import time

#-------------------------------------------------------------------------
#
# Local Imports
#
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#
#Fipy Imports
#-------------------------------------------------------------------------
from fipy import *
from fipy.viewers.viewer import AbstractViewer

#-------------------------------------------------------------------------
#
# show solution
#
#-------------------------------------------------------------------------
ARGS = sys.argv

class MayaviClient(AbstractViewer):
    """
    The `MayaviClient` uses the Mayavi_ python plotting package.

    .. _Mayavi: http://code.enthought.com/projects/mayavi

    """
    
    def __init__(self, file, title=None, daemon_file=None, fps=1.0, **kwlimits):
        """
        Create a `MayaviClient`.
        
        :Parameters:
          file
            a vtk file of a `CellVariable` 
          title
            displayed at the top of the `Viewer` window
          xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
          daemon_file
            the path to the script to run the separate MayaVi viewer process.
            Defaults to "fipy/viewers/mayaviViewer/mayaviDaemon.py"
          fps
            frames per second to attempt to display
        """
        self.fps = fps
        
        self.vtkdir = tempfile.mkdtemp()
        self.vtkcellfname = file
        self.vtkfacefname = None
        self.vtklockfname = os.path.join(self.vtkdir, "lock")

        from fipy.viewers.vtkViewer import VTKCellViewer, VTKFaceViewer

        try:
            self.vtkCellViewer = VTKCellViewer(vars=vars)
            cell_vars = self.vtkCellViewer.vars
        except TypeError:
            self.vtkCellViewer = None
            cell_vars = []

##        try:
##            self.vtkFaceViewer = VTKFaceViewer(vars=vars)
##            face_vars = self.vtkFaceViewer.vars
##        except TypeError:
##            self.vtkFaceViewer = None
##            face_vars = []

##        AbstractViewer.__init__(self, vars=cell_vars + face_vars, title=title, **kwlimits)
        
##        self.plot()

        from pkg_resources import Requirement, resource_filename
        daemon_file = (daemon_file 
                       or resource_filename(Requirement.parse("FiPy"), 
                                "fipy/viewers/mayaviViewer/mayaviDaemon.py"))
        
        cmd = ["python", 
               daemon_file,
               "--lock",
               self.vtklockfname,
               "--fps",
               str(self.fps)]

        if self.vtkcellfname:
            cmd += ["--cell", self.vtkcellfname]
            
##        if self.vtkFaceViewer is not None:
##            cmd += ["--face", self.vtkfacefname]
            
                
        if 'xmin' in kwlimits:
            cmd += ["--xmin" , str(xmin)]
##        cmd += self._getLimit('xmin')
##        cmd += self._getLimit('xmax')
##        cmd += self._getLimit('ymin')
##        cmd += self._getLimit('ymax')
##        cmd += self._getLimit('zmin')
##        cmd += self._getLimit('zmax')
##        cmd += self._getLimit('datamin')
##        cmd += self._getLimit('datamax')

        self.daemon = subprocess.Popen(cmd)
        
    def __del__(self):
        for fname in [self.vtklockfname]:
            if fname and os.path.isfile(fname):
                os.unlink(fname)
        os.rmdir(self.vtkdir)

def main():
    if len(ARGS) < 2 or not os.path.isfile(ARGS[1]):
        print "ERROR: first argument must be file to show."
        sys.exit()
    file = ARGS[1]
    mc = MayaviClient(file)
    raw_input('Press Enter to close')

from fipy.viewers.viewer import AbstractViewer

if __name__ == "__main__":
    main()
