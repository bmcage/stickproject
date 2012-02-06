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
transfer_right = %(tf)s
evap_satconc = %(evapsatc)s
evap_transfer = %(evaptf)s
evap_minbound = %(evapmin)s
out_conc = %(outc)s

[time]
time_period = 50.
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
        for met, smet in [('FVM', 'cvode'), 
                          ('FVM', 'cvode_step'), 
                          ('SIMPLE', 'standard'), 
                          ('FVM', 'odew'),
                          ('FVM', 'odew_step')
                         ]:
            inifile = INI_FIBER % {
                'intdif': 'False',
                'method': "'%s'" % met,
                'submethod': "'%s'" % smet,
                'boundtype': "'flux'",
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
            assert ok, 'Info: start mass %f, end mass %f' % (startmass, endmass)
            FiberConfigManager.delete('fiber.ini')
            del model
            del cfg

if __name__ == "__main__":
    try:
        run_module_suite()
    except NameError:
        test = TestFiberBcs()
        test.test_zero_flux()
