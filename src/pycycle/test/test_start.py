import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group

from pycycle.start import FlowStart
from test_util import assert_rel_error

class StartTestCase(unittest.TestCase):
    def test_start(self): 
        comp = FlowStart()
        g = Group()
        g.add('comp', comp)
        p = Problem(root=g)
        p.setup()

        comp.params['W'] = 3.488
        comp.params['Pt'] = 0.0272
        comp.params['Tt'] = 630.75
        comp.params['Mach'] = 1.0

        p.run()

        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 3.488,.005)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], .0272, .005)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 630.75, .005)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], .000074, .005)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 1.00,.005)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 6060.6, .005)
 
if __name__ == "__main__":
    unittest.main()
    
