import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group

from pycycle.api import Compressor
from test_util import assert_rel_error
from pycycle.flowstation import FlowStation

class CompressorTestCase(unittest.TestCase):

    def setUp(self): 
        self.comp = Compressor()

    def tearDown(self):
        comp = None

    def test_compressor(self):
        comp = self.comp
        g = Group()
        g.add('comp', comp)
        p = Problem(root=g)
        p.setup()

        comp.params['PR_des'] = 12.47
        comp.params['MNexit_des'] = 0.4
        comp.params['eff_des'] = 0.8
        comp.params['flow_in:in:W'] = 1.08
        comp.params['flow_in:in:Tt'] = 630.74523
        comp.params['flow_in:in:Pt'] = 0.0271945
#        comp.solve_totals_TP('flow_in')
#        flowstation.set_total_TP(comp.params, comp.Fl_I_data, 630.74523, 0.0271945, flowstation.SET_BY_NONE)
        comp.params['flow_in:in:Mach'] = 0.6
#        comp.flow_in.solve_statics()
#        flowstation.set_static(comp.params, comp.Fl_I_data, flowstation.SET_BY_Mach)
        comp.params['design'] = True

        p.run()

        TOL = 0.001
        assert_rel_error(self, comp.unknowns['flow_out:W'], 1.08, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:Pt'], 0.33899, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:Tt'], 1424.01, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:rhos'], 0.000594, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:Mach'], 0.4 ,TOL)
        assert_rel_error(self, comp.unknowns['flow_out:area:out'], 364.7, TOL)
        assert_rel_error(self, comp.unknowns['pwr'], 303.2, TOL)
        assert_rel_error(self, comp.unknowns['eff_poly'], 0.8545, TOL)

        # run off design
        comp.params['design'] = False
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
