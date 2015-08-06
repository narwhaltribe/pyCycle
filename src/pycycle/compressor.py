import math 

from openmdao.core.component import Component

from pycycle import flowstation
from pycycle.cycle_component import CycleComponent

class Compressor(CycleComponent): 
    '''Axial Compressor performance calculations''' 

    def __init__(self):
        super(Compressor, self).__init__()
        self.add_param('PR_des', 12.47, desc='Pressure ratio at design conditions')
        self.add_param('MNexit_des', 0.4, desc='mach number at the compressor exit at design conditions')
        self.add_param('eff_des', 0.95, desc='adiabatic efficiency at the design condition')
        self.add_param('hub_to_tip', 0.4, desc='ratio of hub radius to tip radius')
        self.add_param('op_slope', 0.85, desc='slope of operating line (pressure/efficiency)') # TODO add description
        self.Fl_I_data = flowstation.init_fs_tree(self.add_param, 'Fl_I') # incoming air stream to compressor
        self.add_output('PR', 0.0, desc='pressure ratio at operating conditions')
        self.add_output('eff', 0.0, desc='adiabatic efficiency at the operating condition')
        self.add_output('eff_poly', 0.0, desc='polytropic efficiency at the operating condition')
        self.add_output('pwr', 0.0, units='hp', desc='power required to run the compressor at the operating condition')
        self.Fl_O_data = flowstation.init_fs_tree(self.add_output, 'Fl_O') # outgoing air stream from compressor
        self.add_output('tip_radius', 0.0, units='inch', desc='radius at the tip of the compressor')
        self.add_output('hub_radius', 0.0, units='inch', desc='radius at the tip of the compressor')

    def _op_line(self, params, Wc): 
        '''Relationship between compressor pressure ratio and mass flow''' 
        b = 1 - params['op_slope'] # scaled PR and Wc at design are both 1
        # assume a linear op line, with given slope
        norm_PR = params['op_slope'] * (Wc / self._Wc_des) + b 
        return norm_PR*self.PR_des

    def solve_nonlinear(self, params, unknowns, resids):
        fs_ideal = flowstation.init_fs_standalone()
        unknowns['Fl_O:W'] = params['Fl_I:W']
        flowstation.set_static(unknowns, self.Fl_O_data, flowstation.SET_BY_NONE)
        if self.params['design']:
            # Design Calculations
            Pt_out = params['Fl_I:Pt'] * params['PR_des']
            unknowns['PR'] = params['PR_des']
            flowstation.set_total_sP(fs_ideal[0], fs_ideal[1], params['Fl_I:s'], Pt_out, flowstation.SET_BY_NONE)
            ht_out = (fs_ideal[0][':ht'] - params['Fl_I:ht']) / params['eff_des'] + params['Fl_I:ht']
            flowstation.set_total_hP(unknowns, self.Fl_O_data, ht_out, Pt_out, flowstation.SET_BY_NONE)
            unknowns['Fl_O:Mach'] = params['MNexit_des']
            flowstation.set_static(unknowns, self.Fl_O_data, flowstation.SET_BY_Mach)
            self._exit_area_des = unknowns['Fl_O:area']
            self._Wc_des = params['Fl_I:Wc']
        else:
            # Assumed Op Line Calculation
            unknowns['PR'] = self._op_line(params, params['Fl_I:Wc'])
            unknowns['eff'] = parameters['eff_des'] # TODO: add in eff variation with W
            # Operational Conditions
            Pt_out = params['Fl_I:Pt'] * unknowns['PR']
            flowstation.set_total_sP(fs_ideal[0], fs_ideal[1], params['Fl_I:s'], Pt_out, flowstation.SET_BY_NONE)
            ht_out = (fs_ideal[0][':ht'] - params['Fl_I:ht']) / unknowns['eff'] + params['Fl_I:ht']
            flowstation.set_total_hP(unknowns, self.Fl_O_data, ht_out, Pt_out, flowstation.SET_BY_NONE)
            unknowns['Fl_O:area'] = self._exit_area_des # causes Mach to be calculated based on fixed area
        C = flowstation.GAS_CONSTANT * math.log(unknowns['PR'])
        delta_s = unknowns['Fl_O:s'] - params['Fl_I:s']
        unknowns['eff_poly'] = C / (C + delta_s)
        unknowns['pwr'] = params['Fl_I:W'] * (unknowns['Fl_O:ht'] - params['Fl_I:ht']) * 1.4148532 # btu/s to hp
        unknowns['tip_radius'] = (unknowns['Fl_O:area'] / math.pi / (1 - params['hub_to_tip'] ** 2)) ** 0.5
        unknowns['hub_radius'] = params['hub_to_tip'] * unknowns['tip_radius']

if __name__ == '__main__': 
    from openmdao.core.problem import Problem
    from openmdao.core.group import Group
    
    g = Group()
    g.add('comp', Compressor())
    p = Problem(root=g)
    p.setup()
    p.run()
