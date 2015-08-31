import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group
from openmdao.components.constraint import ConstraintComp
from openmdao.components import ParamComp
from openmdao.drivers.scipy_optimizer import ScipyOptimizer

from test_util import assert_rel_error
from pycycle.heat_exchanger import HeatExchanger
from pycycle import flowstation

class HXTestCase(unittest.TestCase):
    def test_heat_exchanger(self): 
        comp = HeatExchanger()
        g = Group()
        g.add('comp', comp)
        p = Problem(root=g)
        p.setup()

        comp.params['flow_in:in:Tt'] = 1423.8
        comp.params['flow_in:in:Pt'] = 0.302712118187
        comp.params['flow_in:in:W'] = 1.0
        comp.params['dPqP'] = 0.0
        comp.params['design'] = True
        p.run()
        
        TOL = 0.005
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 539.94, TOL)
        assert_rel_error(self, comp.unknowns['Qreleased'], 327.22, TOL)
        assert_rel_error(self, comp.unknowns['Qabsorbed'], 327.22, TOL)
        assert_rel_error(self, comp.unknowns['Qmax'], 335.1, TOL)
        assert_rel_error(self, comp.unknowns['T_cold_out'], 749.96, TOL)
        
        #check off design 
        comp.params['design'] = False
        
        p.run()

        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 539.94, TOL)
        assert_rel_error(self, comp.unknowns['Qreleased'], 327.22, TOL)
        assert_rel_error(self, comp.unknowns['Qabsorbed'], 327.22, TOL)
        assert_rel_error(self, comp.unknowns['Qmax'], 335.1, TOL)
        assert_rel_error(self, comp.unknowns['T_cold_out'], 749.96, TOL)
        
if __name__ == "__main__":
    unittest.main()
