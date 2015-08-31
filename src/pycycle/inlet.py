from pycycle.cycle_component import CycleComponent

class Inlet(CycleComponent):
    '''The inlet takes in air at a given flow rate and mach number, and diffuses it down to a slower mach number and larger area'''
    def __init__(self):
        super(Inlet, self).__init__()
        self.add_param('ram_recovery', 1.0, desc='fraction of the total pressure retained')
        self.add_param('MNexit_des', 0.6, desc='Mach number at exit')

        self.add_output('A_capture', 0.0, desc='area at entrance plane', units='inch**2')
        self.add_output('F_ram', 0.0, desc='ram drag', units='lbf')

        self._add_flowstation('flow_in')
        self._add_flowstation('flow_out')

    def solve_nonlinear(self, params, unknowns, resids):
        self._clear_unknowns('flow_in', unknowns)
        self._clear_unknowns('flow_out', unknowns)
        self._solve_flow_vars('flow_in', params, unknowns)
        Pt_out = unknowns['flow_in:out:Pt'] * params['ram_recovery']
        unknowns['flow_out:out:W'] = unknowns['flow_in:out:W']
        unknowns['flow_out:out:Tt'] = unknowns['flow_in:out:Tt']
        unknowns['flow_out:out:Pt'] = Pt_out
        unknowns['F_ram'] = unknowns['flow_in:out:W'] * unknowns['flow_in:out:Vflow'] / 32.174 #lbf
        if params['design']:
            unknowns['flow_out:out:Mach'] = params['MNexit_des']
            self._solve_flow_vars('flow_out', params, unknowns)
            self._exit_area_des = unknowns['flow_out:out:area']
            unknowns['A_capture'] = unknowns['flow_in:out:area']
        else: 
            unknowns['flow_out:out:area'] = self._exit_area_des
            self._solve_flow_vars('flow_out', params, unknowns)
