#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('stick', parent_package, top_path)
    config.add_subpackage('bednet')
    config.add_subpackage('fabric1d')
    config.add_subpackage('fiber')
    config.add_subpackage('fiber1d')
    config.add_subpackage('lib')
    config.add_subpackage('model')
    config.add_subpackage('utils')
    config.add_subpackage('yarn')
    config.add_subpackage('yarn1d')
    config.add_subpackage('yarn2d')
    return config

if __name__ == '__main__':
    print('This is the wrong setup.py file to run')
