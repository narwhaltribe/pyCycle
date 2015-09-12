import unittest

from openmdao.core.group import Group
from openmdao.core.problem import Problem

from test_util import assert_rel_error
from pycycle.components import SplitterBPR, SplitterW

class SplitterTestCase(unittest.TestCase):
    def check(self, comp):
        TOL = .001
        assert_rel_error(self, comp.unknowns['flow_out_1:out:W'], 1.08, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_1:out:Pt'], 0.0271945, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_1:out:Tt'], 630.75, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_1:out:Mach'], 1.0, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_1:out:area'], 1877.2, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_1:out:rhos'], 0.0000737216, TOL)

        assert_rel_error(self, comp.unknowns['flow_out_2:out:W'], 2.407, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_2:out:Pt'], 0.0271945, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_2:out:Tt'], 630.75, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_2:out:Mach'], 1.0, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_2:out:area'], 4183.4, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_2:out:rhos'], 0.0000737216, TOL)

    def test_splitterBPR(self):
        g = Group()
        p = Problem(root=g)
        comp = g.add('comp', SplitterBPR())
        p.setup(check=False)
        
        comp.params['flow_in:in:W'] = 3.48771299
        comp.params['flow_in:in:Tt'] = 630.74523
        comp.params['flow_in:in:Pt'] = 0.0271945
        comp.params['flow_in:in:Mach'] = 1.0

        comp.params['BPR'] = 2.2285
        comp.params['MNexit1_des'] = 1.0
        comp.params['MNexit2_des'] = 1.0
        comp.params['design'] = True
        p.run()
        self.check(comp)

        comp.params['design'] = False
        comp.params['flow_out_1:in:is_super'] = True
        comp.params['flow_out_2:in:is_super'] = True
        p.run()
        self.check(comp)

        comp.params['flow_in:in:W'] *= 0.95
        comp.params['flow_out_1:in:is_super'] = False
        comp.params['flow_out_2:in:is_super'] = False
        p.run()
        TOL = 0.001
        assert_rel_error(self, comp.unknowns['flow_out_1:out:Mach'], 0.76922, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_2:out:Mach'], 0.76922, TOL)

    def test_splitterW(self):
        g = Group()
        p = Problem(root=g)
        comp = g.add('comp', SplitterW())
        p.setup(check=False)
        
        comp.params['W1_des'] = 1.08
        comp.params['MNexit1_des'] = 1.0
        comp.params['MNexit2_des'] = 1.0
        comp.params['design'] = True

        comp.params['flow_in:in:W'] = 3.48771299
        comp.params['flow_in:in:Tt'] = 630.74523
        comp.params['flow_in:in:Pt'] = 0.0271945
        comp.params['flow_in:in:Mach'] = 1.0
        p.run()
        self.check(comp)

        comp.params['design'] = False
        comp.params['flow_out_1:in:is_super'] = True
        comp.params['flow_out_2:in:is_super'] = True
        p.run()
        self.check(comp)

        comp.params['flow_in:in:W'] *= 0.95
        comp.params['flow_out_1:in:is_super'] = False
        comp.params['flow_out_2:in:is_super'] = False
        p.run()
        TOL = 0.001
        assert_rel_error(self, comp.unknowns['flow_out_1:out:Mach'], 0.76922, TOL)
        assert_rel_error(self, comp.unknowns['flow_out_2:out:Mach'], 0.76922, TOL)

if __name__ == '__main__':
    unittest.main()
