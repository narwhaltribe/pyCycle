import math 

from openmdao.core.component import Component

class CycleComponent(Component): 

    def __init__(self): 
        super(CycleComponent, self).__init__()
        self.add_param('design', False, desc='flag to indicate that the calculations are design conditions')
