from pycycle.cycle_component import CycleComponent

class Duct(CycleComponent):
    '''The inlet takes in air at a given flow rate and mach number, and diffuses it down to a slower mach number and larger area'''
    def __init__(self):
        super(Duct, self).__init__()
        self.add_param('dPqP', 0.0, desc='pressure differential as a fraction of incoming pressure')
        self.add_param('Q_dot', 0.0, desc='heat flow rate into (positive) or out of (negative) the air', units='Btu/s')
        self.add_param('MNexit_des', 0.6, desc='Mach number at the exit of the inlet')
        self._add_flowstation('flow_in')
        self._add_flowstation('flow_out')

    def solve_nonlinear(self, params, unknowns, resids):
        self._clear_unknowns('flow_in', unknowns)
        self._clear_unknowns('flow_out', unknowns)
        self._solve_flow_vars('flow_in', params, unknowns)
        Pt_out = unknowns['flow_in:out:Pt'] * (1.0 - params['dPqP'])
        q = params['Q_dot'] / params['flow_in:in:W']
        unknowns['flow_out:out:ht'] = unknowns['flow_in:out:ht'] + q
        unknowns['flow_out:out:Pt'] = Pt_out
        unknowns['flow_out:out:W'] = params['flow_in:in:W']
        if params['design']: 
            unknowns['flow_out:out:Mach'] = self.params['MNexit_des']
            self._solve_flow_vars('flow_out', params, unknowns)
            self._exit_area_des = unknowns['flow_out:out:area']
        else: 
            unknowns['flow_out:out:area'] = self._exit_area_des
            self._solve_flow_vars('flow_out', params, unknowns)
