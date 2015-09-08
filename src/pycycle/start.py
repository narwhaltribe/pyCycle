from pycycle import flowstation
from pycycle.cycle_component import CycleComponent

class FlowStart(CycleComponent):
    '''Flow initialization'''
    def __init__(self):
        super(FlowStart, self).__init__()
        self.add_param('W', 1.0, desc='mass flow rate', units='lbm/s')
        self.add_param('Pt', 14.7, desc='total pressure', units='psi')
        self.add_param('Tt', 518.0, desc='total temperature', units='degR') # TODO standardize units across files
        self.add_param('Mach', 0.1, desc='Mach number')
        self.add_output('area_des', 0.0, desc='flow area at design conditions', units='inch**2')
        self._add_flowstation('flow_out') # outgoing flow at specified conditions

    def solve_nonlinear(self, params, unknowns, resids):
        self._clear_unknowns('flow_out', unknowns)
        unknowns['flow_out:out:Tt'] = params['Tt']
        unknowns['flow_out:out:Pt'] = params['Pt']
        unknowns['flow_out:out:W'] = params['W']
        unknowns['flow_out:out:Mach'] = params['Mach']
        print 'Mach', params['Mach']
        self._solve_flow_vars('flow_out', params, unknowns)
        if params['design']: 
            unknowns['area_des'] = unknowns['flow_out:out:area']

class FlowStartStatic(CycleComponent):
    def __init__(self):
        super(FlowStartStatic, self).__init__()
        self.add_param('W', 1.0, desc='mass flow rate', units='lbm/s')
        self.add_param('Ps', 14.7, desc='static pressure', units='psi')
        self.add_param('Ts', 518.0, desc='static temperature', units='degR')
        self.add_param('Mach', 0.1, desc='Mach number')
        self._add_flowstation('flow_out') # outgoing flow at specified conditions

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['flow_out:out:Ts'] = params['Ts']
        unknowns['flow_out:out:Ps'] = params['Ps']
        unknowns['flow_out:out:Mach'] = params['Mach']
        unknowns['flow_out:out:W'] = params['W']
        self._solve_flow_vars('flow_out', params, unknowns) # TODO solve by TsPsMN
