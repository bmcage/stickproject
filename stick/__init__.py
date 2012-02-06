__all__ = ['test'] + [s for s in dir() if not s.startswith('_')]

#add unittesting
try:
    from numpy.testing import Tester
    test = Tester().test
except:
    #testing could not be loaded, old numpy version
    pass
