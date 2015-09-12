import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group

from pycycle.components import Compressor
from test_util import assert_rel_error
from pycycle import flowstation

class CompressorTestCase(unittest.TestCase):
    def test_compressor(self):
        comp = Compressor()
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
        comp.params['flow_in:in:Mach'] = 0.6
        comp.params['design'] = True

        p.run()

        TOL = 0.001
        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 1.08, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.33899, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 1424.01, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.000594, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.4 ,TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 364.7, TOL)
        assert_rel_error(self, comp.unknowns['pwr'], 303.2, TOL)
        assert_rel_error(self, comp.unknowns['eff_poly'], 0.8545, TOL)

        # run off design
        comp.params['design'] = False
        p.run()

        # values should remain unchanged in off-design at design condition
        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 1.08, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.33899, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 1424.01, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.000594, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.4 ,TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 364.7, TOL)
        assert_rel_error(self, comp.unknowns['pwr'], 303.2, TOL)
        assert_rel_error(self, comp.unknowns['eff_poly'], 0.8545, TOL)

        # try changing something
        comp.params['flow_in:in:W'] *= 1.1
        p.run()
        assert_rel_error(self, comp.unknowns['PR'], 13.52995, TOL)

if __name__ == "__main__":
    unittest.main()
