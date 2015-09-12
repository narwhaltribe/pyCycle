from pycycle.cycle_component import CycleComponent

class SplitterBPR(CycleComponent):
    '''Takes a single incoming air stream and splits it into two separate ones based on a given bypass ratio'''
    def __init__(self):
        super(SplitterBPR, self).__init__()
        self._add_flowstation('flow_in')
        self._add_flowstation('flow_out_1')
        self._add_flowstation('flow_out_2')

        self.add_param('BPR', 2.0, desc='ratio of mass flow in flow_out_1 to flow_out_2')
        self.add_param('MNexit1_des', 0.4, desc='Mach number at design conditions for flow_out_1')
        self.add_param('MNexit2_des', 0.4, desc='Mach number at design conditions for flow_out_2')
        
        self.add_output('BPR_des', 0.0, desc='bypass ratio of splitter at design conditions')

    def solve_nonlinear(self, params, unknowns, resids):
        self._clear_unknowns('flow_in', unknowns)
        self._clear_unknowns('flow_out_1', unknowns)
        self._clear_unknowns('flow_out_2', unknowns)
        self._solve_flow_vars('flow_in', params, unknowns)
        unknowns['flow_out_1:out:W'] = unknowns['flow_in:out:W'] / (params['BPR'] + 1.0)
        unknowns['flow_out_2:out:W'] = unknowns['flow_out_1:out:W'] * params['BPR']
        unknowns['flow_out_1:out:Tt'] = unknowns['flow_out_2:out:Tt'] = unknowns['flow_in:out:Tt']
        unknowns['flow_out_1:out:Pt'] = unknowns['flow_out_2:out:Pt'] = unknowns['flow_in:out:Pt']
        if params['design']:
            unknowns['flow_out_1:out:Mach'] = params['MNexit1_des']
            unknowns['flow_out_2:out:Mach'] = params['MNexit2_des']
            self._solve_flow_vars('flow_out_1', params, unknowns)
            self._solve_flow_vars('flow_out_2', params, unknowns)
            self._exit_area_1_des = unknowns['flow_out_1:out:area']
            self._exit_area_2_des = unknowns['flow_out_2:out:area']
            unknowns['BPR_des'] = params['BPR']
        else:
            unknowns['flow_out_1:out:area'] = self._exit_area_1_des
            unknowns['flow_out_2:out:area'] = self._exit_area_2_des
            self._solve_flow_vars('flow_out_1', params, unknowns)
            self._solve_flow_vars('flow_out_2', params, unknowns)

class SplitterW(CycleComponent):
    '''Takes a single incoming air stream and splits it into two separate ones based on a given mass flow for the Fl_O1'''
    def __init__(self):
        super(SplitterW, self).__init__()
        self._add_flowstation('flow_in')
        self._add_flowstation('flow_out_1')
        self._add_flowstation('flow_out_2')

        self.add_param('W1_des', 0.44, desc='design mass flow in flow_out_1', units='lbm/s')
        self.add_param('MNexit1_des', 0.4, desc='Mach number at design conditions for flow_out_1')
        self.add_param('MNexit2_des', 0.4, desc='Mach number at design conditions for flow_out_2')

    def solve_nonlinear(self, params, unknowns, resids):
        self._clear_unknowns('flow_in', unknowns)
        self._clear_unknowns('flow_out_1', unknowns)
        self._clear_unknowns('flow_out_2', unknowns)
        self._solve_flow_vars('flow_in', params, unknowns)
        unknowns['flow_out_1:out:Tt'] = unknowns['flow_out_2:out:Tt'] = unknowns['flow_in:out:Tt']
        unknowns['flow_out_1:out:Pt'] = unknowns['flow_out_2:out:Pt'] = unknowns['flow_in:out:Pt']
        if params['design']:
            unknowns['flow_out_1:out:W'] = params['W1_des']
            unknowns['flow_out_2:out:W'] = unknowns['flow_in:out:W'] - params['W1_des']
            self._BPR_des = unknowns['flow_out_2:out:W'] / unknowns['flow_out_1:out:W']
            unknowns['flow_out_1:out:Mach'] = params['MNexit1_des']
            unknowns['flow_out_2:out:Mach'] = params['MNexit2_des']
            self._solve_flow_vars('flow_out_1', params, unknowns)
            self._solve_flow_vars('flow_out_2', params, unknowns)
            self._exit_area_1_des = unknowns['flow_out_1:out:area']
            self._exit_area_2_des = unknowns['flow_out_2:out:area']
        else:
            unknowns['flow_out_1:out:W'] = unknowns['flow_in:out:W'] / (self._BPR_des + 1.0)
            unknowns['flow_out_2:out:W'] = unknowns['flow_out_1:out:W'] * self._BPR_des
            unknowns['flow_out_1:out:area'] = self._exit_area_1_des
            unknowns['flow_out_2:out:area'] = self._exit_area_2_des
            self._solve_flow_vars('flow_out_1', params, unknowns)
            self._solve_flow_vars('flow_out_2', params, unknowns)
