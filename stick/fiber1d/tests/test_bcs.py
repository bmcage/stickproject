from __future__ import division
from numpy.testing import TestCase, run_module_suite
from numpy import allclose
from fiber.config import FiberConfigManager
from fiber1d.fibermodel import FiberModel

INI_FIBER = u"""
[general]
method = %(method)s
submethod = %(submethod)s
read = False
verbose = False

[fiber]
radius_pure_fiber = 0.052
form = 'circle'
nrlayers = 2
mean_deviation = 0.0074
internal_diffusion = %(intdif)s
diffusion_coef = 2e-9

[fiberlayer_0]
thickness = 0.00085
n_edge = 41
porosity_layer = 1.0
diffusion_coef = 5.2e-9
diffusion_polymer_exp_factor = 0.
init_conc = 'lambda x: 0.0'

[fiberlayer_1]
thickness = 0.00085
n_edge = 41
porosity_layer = 1.0
diffusion_coef = 7e-9
diffusion_polymer_exp_factor = 0.
init_conc = 'lambda x: 0.7'

[boundary]
type_right = %(boundtype)s
boundary_fib_right = %(boundary_fib_right)s
transfer_right = %(tf)s
evap_satconc = %(evapsatc)s
evap_transfer = %(evaptf)s
evap_minbound = %(evapmin)s
out_conc = %(outc)s

[time]
time_period = 100.
dt = 5.0

[plot]
plotevery = 0
plotmass = False
plotflux = False

"""

class TestFiberBcs(TestCase):
    """
    Check fiber boundary conditions
    """
    def test_zero_flux(self):
        """ we set zero flux cond and check that mass is conserved"""
        for met, smet, inner in [('FVM', 'cvode', 'False'), 
                          ('FVM', 'cvode_step', 'False'), 
                          ('SIMPLE', 'standard', 'False'), 
                          ('FVM', 'odew', 'False'),
                          ('FVM', 'odew_step', 'False'),
                          ('FVM', 'cvode', 'True'), 
                          ('FVM', 'cvode_step', 'True'), 
                          ('SIMPLE', 'standard', 'True'), 
                          ('FVM', 'odew', 'True'),
                          ('FVM', 'odew_step', 'True')
                         ]:
            inifile = INI_FIBER % {
                'intdif': inner,
                'method': "'%s'" % met,
                'submethod': "'%s'" % smet,
                'boundtype': "'flux'",
                'boundary_fib_right': '0.',
                'tf': '0.',
                'evapsatc': "'lambda T: 1.'",
                'evaptf': '0.',
                'evapmin': '0.',
                'outc' : "'lambda t: 0.'",
                }
            cfg = FiberConfigManager.get_instance('fiber.ini', realdatastr=inifile)
            model = FiberModel(cfg)
            #pass further execution to the model
            model.run_init()
            model.solve_init()
            if met == 'SIMPLE':
                startmass = model.simple_sol[0]
            else:
                startmass = model.calc_mass(model.initial_c1)
            model.solve()
            if met == 'SIMPLE':
                endmass = model.simple_sol[-1]
            else:
                endmass = model.calc_mass(model.conc1[-1])
            ok = allclose([startmass], [endmass], atol=1e-7, rtol=1e-4)
            #print 'Info: start mass %f, end mass %f' % (startmass, endmass)
            assert ok, 'Info: start mass %f, end mass %f for %s' % (startmass, endmass, smet)
            FiberConfigManager.delete('fiber.ini')
            del model
            del cfg

    def test_neumann_bc(self):
        """ we set neumann right cond and check that this is indeed the case"""
        for met, smet, inner in [('FVM', 'cvode', 'False'), 
                          ('FVM', 'cvode_step', 'False'), 
                        #  ('SIMPLE', 'standard', 'False'), 
                          ('FVM', 'odew', 'False'),
                          ('FVM', 'odew_step', 'False'),
                         ]:
            neumannleft = 1.
            neumannright = 1.5
            inifile = INI_FIBER % {
                'intdif': inner,
                'method': "'%s'" % met,
                'submethod': "'%s'" % smet,
                'boundtype': "'flux'",
                'boundary_fib_right': '%f' % neumannright,
                'tf': '0.',
                'evapsatc': "'lambda T: 1.'",
                'evaptf': '0.',
                'evapmin': '0.',
                'outc' : "'lambda t: 0.'",
                }
            cfg = FiberConfigManager.get_instance('fiber.ini', realdatastr=inifile)
            # as we start with an init condition that does not satisfy the 
            # bc, we need a longer time before we satisfy bc somewhat
            cfg.set('time.time_period', 400.)
            # we also check left boundary
            cfg.set('boundary.boundary_fib_left', neumannleft)
            model = FiberModel(cfg)
            #pass further execution to the model
            model.run_init()
            model.solve_init()
            if met == 'SIMPLE':
                startmass = model.simple_sol[0]
            else:
                startmass = model.calc_mass(model.initial_c1)
            model.solve()
            if met == 'SIMPLE':
                endmass = model.simple_sol[-1]
                raise NotImplementedError
            else:
                endmass = model.calc_mass(model.conc1[-1])
                
                leftderiv = (model.conc1[-1][1] - model.conc1[-1][0])/((model.delta_r[0]+model.delta_r[1])/2)
                rightderiv = (model.conc1[-1][-1] - model.conc1[-1][-2])/((model.delta_r[-1]+model.delta_r[-2])/2)
            #print model.conc1
            #print model.conc1[-1]
            #print rightderiv, neumann
            ok = allclose([neumannright, neumannleft], [rightderiv, leftderiv], atol=1e-1, rtol=1e-1)
            #print 'Info: start mass %f, end mass %f' % (startmass, endmass)
            assert ok, 'Info: neuman cond %s, deriv %s for %s' % \
                        (str([neumannright, neumannleft]), 
                         str([rightderiv, leftderiv]), smet)
            print 'Info: neuman cond %s, deriv %s for %s' % \
                        (str([neumannright, neumannleft]), 
                         str([rightderiv, leftderiv]), smet)
            FiberConfigManager.delete('fiber.ini')
            del model
            del cfg

if __name__ == "__main__":
    try:
        run_module_suite()
    except NameError:
        test = TestFiberBcs()
        test.test_zero_flux()
        test.test_neumann_bc()
