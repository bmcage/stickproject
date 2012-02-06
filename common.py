descr   = """
stick is a toolbox for textile modelling.

LICENSE: the license is GPLv2 or later.
"""

DISTNAME            = 'stick'
DESCRIPTION         = 'A python module for textile modelling'
LONG_DESCRIPTION    = descr
MAINTAINER          = 'maintainer is B. Malengier'
MAINTAINER_EMAIL    = 'benny.malengier@gmail.org'
URL                 = 'http://cage.ugent.be/~bm/progs.html'
LICENSE             = 'GPL v2 or later'

DOWNLOAD_URL        = URL

MAJOR = 0
MINOR = 1
MICRO = 0
DEV = True

CLASSIFIERS = [
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPLv2 License',
        'Topic :: Scientific/Engineering']

def build_verstring():
    return '%d.%d.%d' % (MAJOR, MINOR, MICRO)

def build_fverstring():
    if DEV:
        return build_verstring() + 'dev'
    else:
        return build_verstring()

def write_version(fname):
    f = open(fname, "w")
    f.writelines("version = '%s'\n" % build_verstring())
    f.writelines("dev =%s\n" % DEV)
    f.writelines("full_version = '%s'\n" % build_fverstring())
    f.close()

VERSION = build_fverstring()
INSTALL_REQUIRE = 'scipy'
