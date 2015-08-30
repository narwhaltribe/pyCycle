import math 

from pycycle import flowstation
from pycycle.flowstation import GAS_CONSTANT
from pycycle.cycle_component import CycleComponent

class Compressor(CycleComponent): 
    '''Basis for axial compressor performance calculations (without flowstations)''' 

    def __init__(self):
        super(Compressor, self).__init__()
        self.add_param('PR_des', 12.47, desc='Pressure ratio at design conditions')
        self.add_param('MNexit_des', 0.4, desc='mach number at the compressor exit at design conditions')
        self.add_param('eff_des', 0.95, desc='adiabatic efficiency at the design condition')
        self.add_param('hub_to_tip', 0.4, desc='ratio of hub radius to tip radius')
        self.add_param('op_slope', 0.85, desc='slope of operating line (pressure/efficiency)')
        self.add_output('PR', 0.0, desc='pressure ratio at operating conditions')
        self.add_output('eff', 0.0, desc='adiabatic efficiency at the operating condition')
        self.add_output('eff_poly', 0.0, desc='polytropic efficiency at the operating condition')
        self.add_output('pwr', 0.0, units='hp', desc='power required to run the compressor at the operating condition')
        self.add_output('tip_radius', 0.0, units='inch', desc='radius at the tip of the compressor')
        self.add_output('hub_radius', 0.0, units='inch', desc='radius at the tip of the compressor')
        self._add_flowstation('flow_in')
        self._add_flowstation('flow_out')

    def _op_line(self, params, Wc):
        '''Relationship between compressor pressure ratio and mass flow''' 
        b = 1.0 - params['op_slope'] # scaled PR and Wc at design are both 1
        # assume a linear op line, with given slope
        norm_PR = params['op_slope'] * (Wc / self._Wc_des) + b 
        return norm_PR * params['PR_des']

    def solve_nonlinear(self, params, unknowns, resids):
        self._clear_unknowns('flow_in', unknowns)
        self._clear_unknowns('flow_out', unknowns)
        self._solve_flow_vars('flow_in', params, unknowns)
        unknowns['flow_out:out:W'] = params['flow_in:in:W']
        if params['design']:
            # Design Calculations
            Pt_out = unknowns['flow_in:out:Pt'] * params['PR_des']
            unknowns['PR'] = params['PR_des']
            ideal_ht = flowstation.solve(s=unknowns['flow_in:out:s'], Pt=Pt_out).ht
            ht_out = (ideal_ht - unknowns['flow_in:out:ht']) / params['eff_des'] + unknowns['flow_in:out:ht']
            unknowns['flow_out:out:ht'] = ht_out
            unknowns['flow_out:out:Pt'] = Pt_out
            unknowns['flow_out:out:Mach'] = params['MNexit_des']
            self._solve_flow_vars('flow_out', params, unknowns)
            self._exit_area_des = unknowns['flow_out:out:area']
            self._Wc_des = unknowns['flow_in:out:Wc']
        else:
            # Assumed Op Line Calculation
            unknowns['PR'] = self._op_line(params, unknowns['flow_in:out:Wc'])
            unknowns['eff'] = params['eff_des'] # TODO: add in eff variation with W
            # Operational Conditions
            Pt_out = unknowns['flow_in:out:Pt'] * unknowns['PR']
            ideal_ht = flowstation.solve(s=unknowns['flow_in:out:s'], Pt=Pt_out).ht
            ht_out = (ideal_ht - unknowns['flow_in:out:ht']) / unknowns['eff'] + unknowns['flow_in:out:ht']
            unknowns['flow_out:out:ht'] = ht_out
            unknowns['flow_out:out:Pt'] = Pt_out
            unknowns['flow_out:out:area'] = self._exit_area_des # causes Mach to be calculated based on fixed area
            self._solve_flow_vars('flow_out', params, unknowns)
            self._solve_flow_vars('flow_out', params, unknowns)
        C = flowstation.GAS_CONSTANT * math.log(unknowns['PR'])
        delta_s = unknowns['flow_out:out:s'] - unknowns['flow_in:out:s']
        unknowns['eff_poly'] = C / (C + delta_s)
        unknowns['pwr'] = params['flow_in:in:W'] * (unknowns['flow_out:out:ht'] - unknowns['flow_in:out:ht']) * 1.4148532 # btu/s to hp
        unknowns['tip_radius'] = (unknowns['flow_out:out:area'] / math.pi / (1 - params['hub_to_tip'] ** 2)) ** 0.5
        unknowns['hub_radius'] = params['hub_to_tip'] * unknowns['tip_radius']
