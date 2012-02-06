#!/usr/bin/env python
import sys
import os
#
from fipy import DiffusionTerm
#------------------------------------------------------------------------------
#
#Update the DiffusionTermNoCorrection from fipy
#
#------------------------------------------------------------------------------

class MyDiffusionTermNoCorrection(DiffusionTerm):
    def _getNormals(self, mesh):
        return mesh._getFaceNormals()