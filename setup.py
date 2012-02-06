#! /usr/bin/env python
"""
Stickproject is a toolbox to model textile related physics
"""

import os
import sys

import setuptools

from common import *

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path,
        license = LICENSE,
        download_url = DOWNLOAD_URL,
        long_description = LONG_DESCRIPTION)
    config.add_subpackage(DISTNAME)
    return config

def setup_package():

    from numpy.distutils.core import setup
    setup(name=DISTNAME,
        version = VERSION,
        maintainer = MAINTAINER,
        maintainer_email = MAINTAINER_EMAIL,
        description = DESCRIPTION,
        url = URL,
        license = LICENSE,
        configuration = configuration,
        install_requires = 'scipy',
        zip_safe = False
        )
    return

if __name__ == '__main__':
    setup_package()
