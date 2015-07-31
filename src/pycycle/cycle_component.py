import math 

from openmdao.core.component import Component

class CycleComponent(Component): 

    def __init__(self): 
        super(CycleComponent, self).__init__()
        self.add_param('design', False, desc='flag to indicate that the calculations are design conditions')
        self.run_design = False

    def _design_fired(self): 
        self.run_design = True

    def solve_nonlinear(self,*args,**kwargs): 
        super(CycleComponent, self).run(*args,**kwargs)
        self.run_design = False
