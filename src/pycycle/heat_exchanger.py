'''
preHeatEx.py -  (Run this before heatExchanger2.py)
    Performs inital energy balance for a basic heat exchanger design

    Originally built by Scott Jones in NPSS, ported and augmented by Jeff Chin

NTU (effectiveness) Method
    Determine the heat transfer rate and outlet temperatures when the type and size of the heat exchanger is specified.

    NTU Limitations
    1) Effectiveness of the chosen heat exchanger must be known (empirical)

    Compatible with OpenMDAO v1.0.5 Alpha
'''

import math

from pycycle.cycle_component import CycleComponent

class HeatExchanger(CycleComponent):
    '''Calculates output temperatures for water and air, and heat transfer, for a given water flow rate for a water-to-air heat exchanger'''
    def __init__(self):
        super(HeatExchanger, self).__init__()
        self.add_param('W_cold', 0.992, desc='mass flow rate of cold water', units='lbm/s')
        self.add_param('Cp_cold', 0.9993, desc='specific heat of cold water', units='Btu/(lbm*R)')
        self.add_param('T_cold_in', 518.58, desc='temperature of water into heat exchanger', units='R')
        self.add_param('effectiveness', 0.9765, desc='heat exchange effectiveness')
        self.add_param('MNexit_des', 0.6, desc='Mach number at exit')
        self.add_param('dPqP', 0.1, desc='pressure differential as a fraction of incoming pressure')
        self.add_param('T_hot_out', 1400.0, desc='temperature of air out of heat exchanger', units='R')
        self.add_param('T_cold_out', 518.0, desc='temperature of water out of heat exchanger' units='R')
        self._add_flowstation('flow_in')

        self.add_output('Qreleased', 0.0, desc='energy released', units='hp')
        self.add_output('Qabsorbed', 0.0, desc='energy absorbed', units='hp')
        self.add_output('LMTD', 0.0, desc='logarithmic mean temperature difference')
        self.add_output('Qmax', 0.0, desc='theoretical maximum possible heat transfer', units='hp')
        self.add_output('resid_Qmax', 0.0, desc='residual of max*effectiveness', units='hp')
        self.add_output('resid_e_balance', 0.0, desc='residual of energy balance', units='hp')
        self._add_flowstation('flow_out')

    def solve_nonlinear(self, params, unknowns, resids):
        self._solve_flow_vars('flow_in', params, unknowns)
        W_cold_Cp_Min = min(params['W_cold'] * params['Cp_cold'], unknowns['flow_in:out:W'] * unknowns['flow_in:out:Cp'])
        unknowns['Qmax'] = W_cold_Cp_Min * (unknowns['flow_in:out:Tt'] - params['T_cold_in']) * 1.4148532 #BTU/s to hp
        unknowns['Qreleased'] = unknowns['flow_in:out:W'] * unknowns['flow_in:out:Cp'] * (unknowns['flow_in:out:Tt'] - params['T_hot_out']) * 1.4148532
        unknowns['Qabsorbed'] = params['W_cold'] * params['Cp_cold'] * (params['T_cold_out'] - params['T_cold_in']) * 1.4148532
        try: 
            unknowns['LMTD'] = (params['T_hot_out'] - unknowns['flow_in:out:Tt'] + params['T_cold_out'] - params['T_cold_in']) / math.log((params['T_hot_out'] - params['T_cold_in']) / (unknowns['flow_in:out:Tt'] - params['T_cold_out']))
        except ZeroDivisionError: 
            unknowns['LMTD'] = 0.0
        unknowns['resid_Qmax'] = unknowns['Qreleased'] - params['effectiveness'] * unknowns['Qmax']
        unknowns['resid_e_balance'] = unknowns['Qreleased'] - unknowns['Qabsorbed']

        unknowns['flow_out:out:Tt'] = params['T_hot_out']
        unknowns['flow_out:out:Pt'] = unknowns['flow_in:out:Pt'] * (1.0 - unknowns['dPqP'])
        unknowns['flow_out:out:W'] = unknowns['flow_in:out:W']
        if params['design']:
            unknowns['flow_out:out:Mach'] = params['MNexit_des']
            self._solve_flow_vars('flow_out', params, unknowns)
            self._exit_area_des = unknowns['flow_out:out:area']
        else: 
            unknowns['flow_out:out:area'] = self._exit_area_des
            self._solve_flow_vars('flow_out', params, unknowns)
