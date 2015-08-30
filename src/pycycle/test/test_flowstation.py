import unittest

from openmdao.core.problem import Problem
from openmdao.core.group import Group

from test_util import assert_rel_error
from pycycle import flowstation

from pycycle.cycle_component import CycleComponent
from pycycle.flowstation import GAS_CONSTANT

class DummyComp(CycleComponent):
    def __init__(self):
        super(DummyComp, self).__init__()
        self._add_flowstation('flow')

    def solve_nonlinear(self, params, unknowns, resids):
        self._solve_flow_vars('flow', params, unknowns)

class CycleComponentTestCase(unittest.TestCase):
    def setUp(self): 
        '''Initialization function called before every test function''' 
        self.g = Group()
        self.comp1 = DummyComp()
        self.comp2 = DummyComp()
        self.g.add('comp1', self.comp1)
        self.g.add('comp2', self.comp2)
        self.p = Problem(root=self.g)
        self.p.setup()

    def tearDown(self):
        '''Function called after every test function'''
        self.g = None
        self.p = None

    def test_copy_flow(self): 
        self.comp1.params['flow:in:W'] = 100.0
        self.comp1.params['flow:in:Tt'] = 518.0
        self.comp1.params['flow:in:Pt'] = 15.0
        CycleComponent.copy_from(self.comp1, 'flow', self.comp2, 'flow')

        self.p.run()

        TOL = 0.0001
        assert_rel_error(self, self.comp2.unknowns['flow:out:Tt'], 518.0, TOL)
        assert_rel_error(self, self.comp2.unknowns['flow:out:Pt'], 15.0, TOL)

class FlowStationTestCase(unittest.TestCase):
#        self.fs = FlowStation()
#
#        self.fs.W = 100
#        self.fs.setDryAir()
#        self.fs.setTotalTP(518, 15)

    def test_solve_total_TP(self):
        flow = flowstation.solve(Pt=15.0, Tt=518.0, W=100.0)

        self.assertAlmostEqual(flow.Pt, 15.0, places=2)
        self.assertAlmostEqual(flow.Tt, 518, places=2)
        self.assertAlmostEqual(flow.ht, -6.32357, places=4) #Tom says the ht values will be different
        self.assertAlmostEqual(flow.rhot, .07812, places=4)

    def test_solve_total_hP(self):
        flow = flowstation.solve(ht=-6.32357, Pt=15.0, W=100.0)

        self.assertAlmostEqual(flow.Pt, 15.0, places=2)
        self.assertAlmostEqual(flow.Tt, 518, places=2)
        self.assertAlmostEqual(flow.ht, -6.32357, places=4) #Tom says the ht values will be different
        self.assertAlmostEqual(flow.rhot, .07812, places=4)

    def test_solve_total_sP(self):
        s = flowstation.solve(Pt=15.0, Tt=518.0, W=100.0).s
        flow = flowstation.solve(s=s, Pt=15.0, W=100.0)

        self.assertAlmostEqual(flow.Pt, 15.0, places=2)
        self.assertAlmostEqual(flow.Tt, 518, places=2)
        self.assertAlmostEqual(flow.ht, -6.32357, places=4) #Tom says the ht values will be different
        self.assertAlmostEqual(flow.rhot, .07812, places=4)
    
    def test_delh(self):
        ht = flowstation.solve(Pt=15.0, Tt=518.0, W=100.0).ht
        flow = flowstation.solve(Tt=1000.0, Pt=40.0, W=100.0)
        diffh = flow.ht - ht

        assert_rel_error(self, diffh, 117.4544, .0001)
        
    def test_dels(self):
        s = flowstation.solve(Pt=15.0, Tt=518.0, W=100.0).s
        flow = flowstation.solve(Tt=1000.0, Pt=40.0, W=100.0)
        diffs = flow.s - s
        assert_rel_error(self, diffs, .092609, .0001)
        
# TODO implement the following:
#    def test_set_WAR(self):
#        self.fs.setWAR( 0.02 )
#        self.fs.setTotalTP(1000, 15); 
#        assert_rel_error(self,self.fs.Pt, 15., .0001)
#        assert_rel_error(self,self.fs.Tt, 1000, .0001)
#        assert_rel_error(self,self.fs.WAR, 0.02, .0001)
#        assert_rel_error(self,self.fs.FAR, 0, .0001)
#        assert_rel_error(self,self.fs.ht, -.11513, .0001)
#
#    def test_setDryAir(self):
#        self.fs.setDryAir()
#        self.fs.setTotalTP(1000, 15); 
#        assert_rel_error(self,self.fs.Pt, 15.,.0001)
#        assert_rel_error(self,self.fs.Tt, 1000, .0001)
#        assert_rel_error(self,self.fs.WAR, 0.0, .0001)
#        assert_rel_error(self,self.fs.FAR, 0, .0001)
#        assert_rel_error(self,self.fs.ht, 111.129, .0001)
#        assert_rel_error(self,self.fs.WAR, 0, .0001)
#        assert_rel_error(self,self.fs.FAR, 0, .0001)



#class TestBurn(unittest.TestCase): 
#    def setUp(self):
#        self.fs = FlowStation()
#        self.fs.setDryAir()
#        self.fs.setTotalTP(1100, 400)
#        self.fs.W = 100
#        self.fs.burn(4,2.5, -642)  
#        
#    #all test cases use the same checks here, so just re-use
#    def _assert(self): 
#        #print (self.fs._flow)  assert_rel_error(self,self.fs.W, 102.5, places=2)
#        assert_rel_error(self,self.fs.FAR, .025, .0001)   
#        assert_rel_error(self,self.fs.Pt, 400, .0001)
#        assert_rel_error(self,self.fs.Tt, 2669.69, .0001)
#        assert_rel_error(self,self.fs.ht, 117.171, .0001) 
#        assert_rel_error(self,self.fs.rhot, .404265, .0001)
#        assert_rel_error(self,self.fs.W, 102.5, .0001)
#        assert_rel_error(self,self.fs.gamt, 1.293336, .0001)
#
#    def test_burn(self): 
#
#        self._assert()
#
#        
#    def test_add( self ):
#        self.fs1 = FlowStation()
#
#        self.fs1.setDryAir()
#        self.fs1.setTotalTP(1100, 15)
#        self.fs1.W = 100.
#        self.fs1.setWAR( .02 )
#        self.fs1.setTotalTP(1100, 400)
#        self.fs.add( self.fs1 )
#        assert_rel_error(self,self.fs.Tt, 1932.471, .0001)
#        assert_rel_error(self,self.fs.W, 202.5, .0001)
#        assert_rel_error(self,self.fs.FAR, .012623, .001)   
#        assert_rel_error(self,self.fs.Pt, 400, .001)
#        assert_rel_error(self,self.fs.ht, 71.83056, .0001)
#        assert_rel_error(self,self.fs.gamt, 1.3200, .0001)
#        assert_rel_error(self,self.fs.WAR, .00990099, .0001)  
                
class TestStatics(unittest.TestCase): 
    def setUp(self):
        self.W = 100.
#        self.fs.setDryAir()
        self.Tt = 1100.0
        self.Pt = 400.0

    # all test cases use the same checks here, so just re-use
    def _assert(self, flow):
        assert_rel_error(self, flow.area, 32.006, .0001)
        assert_rel_error(self, flow.Mach, .3, .0001)
        assert_rel_error(self, flow.Ps, 376.194, .0001)
        assert_rel_error(self, flow.Ts, 1081.781, .0001)
        assert_rel_error(self, flow.Vflow, 479.519, .0001)
        assert_rel_error(self, flow.rhos, .93826, .0001)
        assert_rel_error(self, flow.gams, 1.37596, .0001)

    def test_solve_Mach(self):
        flow = flowstation.solve(W=self.W, Tt=self.Tt, Pt=self.Pt, Mach=0.3)
        self._assert(flow)

    def test_solve_area(self):
        flow = flowstation.solve(W=self.W, Tt=self.Tt, Pt=self.Pt, area=32.006)
        self.assertLess(flow.Mach, 1.0)
        self._assert(flow)

    def test_solve_Ps(self):
        flow = flowstation.solve(W=self.W, Tt=self.Tt, Pt=self.Pt, Ps=376.194)
        self._assert(flow)
        
# TODO implement the following:
#    def test_setStaticTsPsMN(self):
#        self.fs.setStaticTsPsMN(1081.802, 376.219, .3)
#        self._assert()

    def test_solve_sub(self):
        flow = flowstation.solve(W=self.W, Tt=self.Tt, Pt=self.Pt, area=32.006, is_super=False)
        self.assertLess(flow.Mach, 1.0)

    def test_solve_super(self):
        flow = flowstation.solve(W=self.W, Tt=self.Tt, Pt=self.Pt, area=32.006, is_super=True)
        self.assertGreater(flow.Mach, 1.0)
        
# TODO implement the following
#class HotH2(unittest.TestCase): 
#    def testMN( self ):
#        self.fs = FlowStation()
#        self.fs.setReactant(6)
#        self.fs.W = 100
#        self.fs.setTotalTP( 5000, 540 );
#        self.fs.Mach = 3.   
#        assert_rel_error(self,self.fs.Mach, 3., .0001)   
#        assert_rel_error(self,self.fs.Ts, 2052.78, .0001)   
#        assert_rel_error(self,self.fs.Ps, 13.032, .0001)
#        assert_rel_error(self,self.fs.Vflow, 24991.2, .0001)
#        assert_rel_error(self,self.fs.rhos, .001192, .0001)        
#        assert_rel_error(self,self.fs.gams, 1.370663, .0001)         

if __name__ == "__main__":
    unittest.main()
    
