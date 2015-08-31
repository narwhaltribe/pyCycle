import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group

from test_util import assert_rel_error
from pycycle.inlet import Inlet
from pycycle import flowstation

class InletTestCase(unittest.TestCase):
    def test_inlet(self):
        comp = Inlet()
        g = Group()
        g.add('comp', comp)
        p = Problem(root=g)
        p.setup()

        comp.params['ram_recovery'] = 1.0
        comp.params['MNexit_des'] = 0.6
        comp.params['flow_in:in:W'] = 1.08
        comp.params['flow_in:in:Tt'] = 630.75
        comp.params['flow_in:in:Pt'] = 0.0272
        comp.params['flow_in:in:Mach'] = 1.0
        comp.params['design'] = True
        p.run()

        TOL = 0.005
        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 1.080, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.0272, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 630.75, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.000098, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.6, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 2230.8, TOL)

        #check off design
        comp.params['design'] = False
        p.run()

        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 1.080, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.0272, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 630.75, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.000098, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.6, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 2230.8, TOL)

        #vary something
        comp.params['flow_in:in:W'] = 0.9
        p.run()

        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 0.9, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.0272, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 630.75, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.45955, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 2230.8, TOL)
        
if __name__ == "__main__":
    unittest.main()
