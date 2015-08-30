import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group

from test_util import assert_rel_error
from pycycle.duct import Duct
from pycycle import flowstation

class DuctTestCase(unittest.TestCase):
    def test_duct(self): 
        comp = Duct()
        g = Group()
        g.add('comp', comp)
        p = Problem(root=g)
        p.setup()

        comp.params['dPqP'] = 0.0
        comp.params['Q_dot'] = -237.8
        comp.params['MNexit_des'] = 0.4
        comp.params['flow_in:in:W'] = 1.080
        comp.params['flow_in:in:Tt'] = 1424.01
        comp.params['flow_in:in:Pt'] = 0.34
        comp.params['flow_in:in:Mach'] = 0.4
        comp.params['design'] = True

        p.run()
        
        TOL = 0.005
        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 1.080, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.34, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 540.0, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.001566, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.4, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 221.4, TOL)

        # check off design 
        comp.params['design'] = False
        p.run()
        
        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 1.080, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.34, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 540.0, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.001566, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.4, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 221.4, TOL)

        # vary something
        comp.params['dPqP'] = 0.1
        p.run()

        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 1.080, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.306, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 540.0, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.0013783, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 0.4572, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 221.4, TOL)
        
if __name__ == "__main__":
    unittest.main()
    
