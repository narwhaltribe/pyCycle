import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group
from openmdao.drivers.scipy_optimizer import ScipyOptimizer
from openmdao.components.exec_comp import ExecComp
from openmdao.components.param_comp import ParamComp
from openmdao.components.constraint import ConstraintComp

from pycycle.nozzle import Nozzle
from pycycle.start import FlowStart
from pycycle.cycle_component import CycleComponent
from pycycle import flowstation
from test_util import assert_rel_error

TOL = 0.001

class NozzleTestCaseResonable(unittest.TestCase):
    def setUp(self):
        self.comp = comp = Nozzle()
        g = Group()
        g.add('comp', comp)
        self.p = p = Problem(root=g)
        p.setup(check=False)
        
        comp.params['flow_in:in:W'] = 100.0
        comp.params['flow_in:in:Tt'] = 700.0
        comp.params['flow_in:in:Pt'] = 50.0
        comp.params['flow_in:in:Mach'] = 0.40
        comp.params['back_Ps'] = 15.0
        comp.params['design'] = True

        p.run()

    def tearDown(self):
        self.comp = None
        self.p = None

    def test_nozzle_off_design(self):
        assert_rel_error(self, self.comp.unknowns['flow_out:out:W'], 100.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Pt'], 50.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Tt'], 700.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Mach'], 1.432, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:area'], 112.88, TOL)
        assert_rel_error(self, self.comp.unknowns['Athroat_des'], 99.59, TOL)
        assert_rel_error(self, self.comp.unknowns['Aexit_des'], 112.88, TOL)

        #off design calcs
        self.comp.params['design'] = False

        self.p.run()

        assert_rel_error(self, self.comp.unknowns['flow_out:out:W'], 100.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Pt'], 50.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Tt'], 700.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Mach'], 1.432, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:area'], 112.88, TOL)
        assert_rel_error(self, self.comp.unknowns['Athroat_des'], 99.59, TOL)
        assert_rel_error(self, self.comp.unknowns['Aexit_des'], 112.88, TOL)

    def test_nozzle_under(self):
        self.comp.params['back_Ps'] = 14.0
        self.comp.params['design'] = False

        self.p.run()
        
        assert_rel_error(self, self.comp.unknowns['flow_out:out:W'], 100.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Pt'], 50.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Tt'], 700.0, TOL)
        assert_rel_error(self, self.comp.unknowns['Athroat_dmd'], 99.59, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:area'], 112.88, TOL)
        self.assertEqual(self.comp.unknowns['switchRegime'], 'UNDEREXPANDED')

    def test_nozzle_over(self):
        self.comp.params['back_Ps'] = 16.0
        self.comp.params['design'] = False
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'], 'OVEREXPANDED')
        assert_rel_error(self, self.comp.unknowns['flow_out:out:W'], 100.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Pt'], 50.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Tt'], 700.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Ps'], 15.0, TOL)
        
        self.comp.params['back_Ps'] = 19.0
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'], 'OVEREXPANDED')
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Ps'], 15.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Mach'], 1.432, TOL)
        
        self.comp.params['back_Ps'] = 24.0
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'], 'OVEREXPANDED')
        
        self.comp.params['back_Ps'] = Pt=27.0
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'],'OVEREXPANDED')

        self.comp.params['back_Ps'] = 32.0
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'],'OVEREXPANDED')

    def test_nozzle_normal_shock(self):
        self.comp.params['back_Ps'] = 35.0
        self.comp.params['design'] = False
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'], 'NORMAL_SHOCK')
        assert_rel_error(self, self.comp.unknowns['flow_out:out:W'], 100.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Pt'], 50.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Tt'], 700., TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Ps'], 35.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Mach'], 0.697, TOL)
        
        self.comp.params['back_Ps'] = 37.0
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'], 'NORMAL_SHOCK')
        assert_rel_error(self, self.comp.unknowns['flow_out:out:W'], 100.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Pt'], 50.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Tt'], 700.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Ps'], 37.0, TOL)
        assert_rel_error(self, self.comp.unknowns['flow_out:out:Mach'], 0.662, TOL)

    def test_nozzle_perfect_expand(self):
        self.comp.params['design'] = False
        self.comp.params['back_Ps'] = 14.99
        self.p.run()
        self.assertEqual(self.comp.unknowns['switchRegime'], 'PERFECTLY_EXPANDED')

#class NozzleTestCaseMassFlowIter(unittest.TestCase):
#    def test_mass_flow_iter(self):
#        g = Group()
#        start = g.add('start', FlowStart())
#        ref = g.add('ref', FlowStart())
#        nozzle = g.add('nozzle', Nozzle())
#
#        CycleComponent.connect_flows(g, 'start.flow_out', 'nozzle.flow_in')
#        CycleComponent.connect_flows(g, 'ref.flow_out', 'nozzle.flow_ref')
#
#        g.add('W_start', ParamComp('W', 100.0))
#        g.connect('W_start.W', 'start.W')
#        g.add('WqAexit_resid', ExecComp('r = WqAexit - WqAexit_dmd', WqAexit=0.0, WqAexit_dmd=0.0))
#        g.connect('nozzle.WqAexit', 'WqAexit_resid.WqAexit')
#        g.connect('nozzle.WqAexit_dmd', 'WqAexit_resid.WqAexit_dmd')
#
#        p = Problem(root=g)
#        p.setup(check=False)
#
#        start.params['W'] = 100.0
#        start.params['Tt'] = 700.0
#        start.params['Pt'] = 50.0
#        start.params['Mach'] = 0.4
#
#        ref.params['W'] = 1.0
#        ref.params['Tt'] = 518.67
#        ref.params['Pt'] = 15.0
#        ref.params['Mach'] = 0.0
#        comp.params['back_Ps'] = flowstation.solve(W=1.0, Tt=518.67, Pt=15.0).Ps
#
#        start.params['design'] = ref.params['design'] = nozzle.params['design'] = True
#
#        p.run()
#
#        p.driver = ScipyOptimizer()
#        p.driver.options['optimizer'] = 'COBYLA'
#        p.driver.add_param('W_start.W', low=1e-15, high=1000)
#        p.driver.add_objective('WqAexit_resid.r')
#        
#        p.setup(check=False)
#
#        start.params['Tt'] = 700.0
#        start.params['Pt'] = 50.0
#        start.params['Mach'] = 0.4
#
#        ref.params['W'] = 1.0
#        ref.params['Tt'] = 518.67
#        ref.params['Mach'] = 0.0
#        comp.params['back_Ps'] = flowstation.solve(W=1.0, Tt=518.67, Pt=15.0).Ps
#
#        start.params['design'] = ref.params['design'] = nozzle.params['design'] = False
#        ref.params['Pt'] = 39.0
#
#        p.run()
#
#        TOL = 0.001
#        self.assertEqual(nozzle.unknowns['switchRegime'], 'UNCHOKED')
#        assert_rel_error(self, nozzle.unknowns['flow_out:out:W'], 96.03, TOL)
#        assert_rel_error(self, nozzle.unknowns['flow_out:out:Mach'], 0.607, TOL)
#        
#        # set W = 80.80, MN throat = 0.562, MN exit = 0.470
#        comp.params['back_Ps'] = flowstation.solve(W=1.0, Tt=518.67, Pt=15.0).Ps
#        ref.params['Pt'] = 43.0
#
#        p.run()
#        
#        self.assertEqual(nozzle.unknowns['switchRegime'], 'UNCHOKED')
#        assert_rel_error(self, nozzle.unknowns['flow_out:out:W'], 80.80, TOL)
#        assert_rel_error(self, nozzle.unknowns['flow_out:out:Mach'], 0.47, TOL)

class NozzleTestCaseVeryLowTemp(unittest.TestCase):
    def test_nozzle_very_low_temperatures(self): 
        comp = Nozzle()
        g = Group()
        g.add('comp', comp)
        p = Problem(root=g)
        p.setup(check=False)

        comp.params['flow_in:in:W'] = 0.639
        comp.params['flow_in:in:Tt'] = 540.0
        comp.params['flow_in:in:Pt'] = 0.34
        comp.params['flow_in:in:Mach'] = 0.4
        comp.params['back_Ps'] = 0.0272
        comp.params['design'] = True
        
        p.run()

        TOL = 0.01 #this test needs larger tollerance due to exteremely low temperatures
        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 0.639, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.34, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 540.0, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 2.7092, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 264.204, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], .000177443, TOL)
        assert_rel_error(self, comp.unknowns['Fg'], 38.98, TOL)

        #off design calcs
        comp.params['design'] = False
        p.run()

        self.assertEqual(comp.unknowns['switchRegime'], 'PERFECTLY_EXPANDED')
        assert_rel_error(self, comp.unknowns['flow_out:out:W'], 0.639, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Pt'], 0.34, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Tt'], 540.0, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:Mach'], 2.7092, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:area'], 264.204, TOL)
        assert_rel_error(self, comp.unknowns['flow_out:out:rhos'], 0.000177443, TOL)

        comp.params['back_Ps'] = 0.03
        p.run()

        self.assertEqual(comp.unknowns['switchRegime'], 'OVEREXPANDED')

        comp.params['back_Ps'] = 0.026
        p.run()

        self.assertEqual(comp.unknowns['switchRegime'], 'UNDEREXPANDED')
 
if __name__ == "__main__":
    unittest.main()
