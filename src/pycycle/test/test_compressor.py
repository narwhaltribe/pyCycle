import unittest

from openmdao.core.problem import Problem

from pycycle.api import Compressor
from test_util import assert_rel_error
import pycycle.flowstation

class CompressorTestCase(unittest.TestCase):

    def setUp(self): 
        self.comp = Compressor()

    def tearDown(self):
        comp = None

    def test_compressor(self):
        comp = self.comp

        comp.params['PR_des'] = 12.47
        comp.params['MNexit_des'] = 0.4
        comp.params['eff_des'] = 0.8
        self.comp.params['Fl_I:W'] = 1.08
        flowstation.setTotalTP(self.comp.params, self.comp.Fl_I_data, 630.74523, 0.0271945)
        self.comp.params['Fl_I:Mach'] = 0.6
        comp.params['design'] = True

        p = Problem(root=comp)
        p.setup()
        p.run()

        TOL = 0.001
        assert_rel_error(self,comp.Fl_O.W, 1.08, TOL)
        assert_rel_error(self,comp.Fl_O.Pt, 0.33899, TOL)
        assert_rel_error(self,comp.Fl_O.Tt, 1424.01, TOL)
        assert_rel_error(self,comp.Fl_O.rhos, 0.000594, TOL)
        assert_rel_error(self,comp.Fl_O.Mach, 0.4 ,TOL)
        assert_rel_error(self,comp.Fl_O.area, 364.7, TOL)
        assert_rel_error(self,comp.pwr, 303.2, TOL)
        assert_rel_error(self,comp.eff_poly, 0.8545, TOL)

        # run off design
        p = Problem(root=comp)
        p.setup()
        p.run()

        # values should remain unchanged in off-design at design condition
        assert_rel_error(self,comp.Fl_O.W, 1.08,TOL)
        assert_rel_error(self,comp.Fl_O.Pt, .33899, TOL)
        assert_rel_error(self,comp.Fl_O.Tt, 1424.01, TOL)
        assert_rel_error(self,comp.Fl_O.rhos, .000594, TOL)
        assert_rel_error(self,comp.Fl_O.Mach, .4 ,TOL)
        assert_rel_error(self,comp.Fl_O.area, 364.7, TOL)
        assert_rel_error(self,comp.pwr, 303.2, TOL)
        assert_rel_error(self,comp.eff_poly, .8545, TOL)

        # try changing something
        comp.params['Fl_I:W'] *= 1.1
        p = Problem(root=comp)
        p.setup()
        p.run()
        assert_rel_error(self, comp.PR, 13.52995, TOL)

if __name__ == "__main__":
    unittest.main()
