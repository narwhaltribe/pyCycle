from pycycle import flowstation
from pycycle.cycle_component import CycleComponent

class Nozzle(CycleComponent): 
    '''Calculates the gross thrust for a convergent-divergent nozzle, assuming an ideally expanded exit condition'''
    def __init__(self):
        super(Nozzle, self).__init__()
        self._add_flowstation('flow_ref')
        self._add_flowstation('flow_in')
        self._add_flowstation('flow_out')

        self.add_param('dPqP', 0.0, desc='ratio of change in total pressure to incoming total pressure')

        self.add_output('Athroat_dmd', 0.0, desc='demand throat area at operating conditions')
        self.add_output('Athroat_des', 0.0, desc='throat area at design conditions')
        self.add_output('Aexit_des', 0.0, desc='exit area at design conditions', units='inch**2')
        self.add_output('Fg', 0.0, desc='gross thrust', units='lbf')
        self.add_output('PR', 0.0, desc='ratio between total and static pressures at exit')
        self.add_output('AR', 0.0, desc='ratio of exit area to throat area')
        self.add_output('WqAexit', 0.0, desc='mass flow per unit area at operating conditions', units='lbm/(s*inch**2)')
        self.add_output('WqAexit_dmd', 0.0, desc='demand mass flow per unit area at operating conditions', units='lbm/(s*inch**2)')
        self.add_output('switchRegime', '', desc='operating regime') # 'UNCHOKED', 'NORMAL_SHOCK', 'UNDEREXPANDED', 'PERFECTLY_EXPANDED', or 'OVEREXPANDED'

    @staticmethod
    def shockPR(Mach, gam):
        '''Calculates stagnation pressure ratio across a normal shock wave'''
        rhoR = (gam + 1) / 2 * Mach ** 2 / (1 + (gam - 1) / 2 * Mach ** 2) # density ratio
        PsR_recip = (gam + 1) / (2 * gam * Mach ** 2 - (gam - 1)) # reciprocal of static pressure ratio
        return rhoR ** (gam / (gam - 1)) * PsR_recip ** (1 / (gam - 1))

    def solve_nonlinear(self, unknowns, params, resids): 
        self._clear_unknowns('flow_ref', unknowns)
        self._clear_unknowns('flow_in', unknowns)
        self._clear_unknowns('flow_out', unknowns)
        self._solve_flow_vars('flow_ref', params, unknowns)
        self._solve_flow_vars('flow_in', params, unknowns)
        Pt_out = (1.0 - params['dPqP']) * unknowns['flow_in:out:Pt']
        flow_throat = flowstation.solve(Tt=unknowns['flow_in:out:Tt'], Pt=Pt_out, Mach=1.0, W=unknowns['flow_in:out:W'])
        unknowns['Athroat_dmd'] = flow_throat.area
        unknowns['flow_out:out:W'] = unknowns['flow_in:out:W']
        if params['design']:
            flow_exit_ideal = flowstation.solve(W=params['flow_out:in:W'], Tt=unknowns['flow_in:out:W'], Pt=Pt_out, Ps=unknowns['flow_ref:out:Ps'])
            unknowns['flow_out:out:Tt'] = unknowns['flow_in:out:Tt']
            unknowns['flow_out:out:Pt'] = Pt_out
            unknowns['flow_out:out:Mach'] = flow_exit_ideal.Mach
            self._solve_flow_vars('flow_out', params, unknowns)
            # Design Calculations at throat
            unknowns['Athroat_des'] = flow_throat.area
            # Design calculations at exit
            unknowns['Aexit_des'] = flow_exit_ideal.area
            unknowns['switchRegime'] = 'PERFECTLY_EXPANDED'
        else:
            # Find subsonic solution, curve 4
            flow_out_subsonic = flowstation.solve(W=unknowns['flow_out:out:W'], Tt=unknowns['flow_in:out:Tt'], Pt=unknowns['flow_in:out:Pt'], is_super=False, area=unknowns['Aexit_des'])
            if flow_out_subsonic.Mach > 1.0:
                raise Exception('invalid nozzle subsonic solution')
            # Find supersonic solution, curve 5
            flow_out_supersonic = flowstation.solve(W=unknowns['flow_out:out:W'], Tt=unknowns['flow_in:out:Tt'], Pt=unknowns['flow_in:out:Pt'], is_super=True, area=unknowns['Aexit_des'])
            # normal shock at nozzle exit, curve c
            Pt_exit = Nozzle.shockPR(flow_out_supersonic.Mach, flow_throat.gams) * flow_throat.Pt
            flow_out_shock = flowstation.solve(W=unknowns['flow_out:out:W'], Tt=flow_throat.Tt, Pt=Pt_exit, is_super=False, area=unknowns['Aexit_des'])
            # find correct operating regime
            # curves 1 to 4
            if unknowns['flow_ref:out:Ps'] >= flow_out_subsonic.Ps:
                unknowns['switchRegime'] = 'UNCHOKED'
                flow_throat = flowstation.solve(Tt=unknowns['flow_in:out:Tt'], Pt=Pt_out, W=unknowns['flow_in:out:W'], area=unknowns['Athroat_des'])
                unknowns['flow_out:out:Tt'] = flow_throat.Tt
                unknowns['flow_out:out:Pt'] = flow_throat.Pt
                unknowns['flow_out:out:area'] = unknowns['Aexit_des']
            # between curves 4 and c
            elif unknowns['flow_ref:out:Ps'] >= flow_out_shock.Ps:
                unknowns['switchRegime'] = 'NORMAL_SHOCK'
                unknowns['flow_out:out:Ps'] = unknowns['flow_ref:out:Ps']
            # between curves c and 5
            elif unknowns['flow_ref:out:Ps'] > flow_out_supersonic.Ps:
                unknowns['switchRegime'] = 'OVEREXPANDED'
                unknowns['flow_out:out:is_super'] = True
                unknowns['flow_out:out:Tt'] = flow_throat.Tt
                unknowns['flow_out:out:Pt'] = flow_throat.Pt
                unknowns['flow_out:out:area'] = unknowns['Aexit_des']
            # between curves 5 and e
            else:
                unknowns['switchRegime'] = 'UNDEREXPANDED'
                unknowns['flow_out:out:is_super'] = True
                unknowns['flow_out:out:Tt'] = flow_throat.Tt
                unknowns['flow_out:out:Pt'] = flow_throat.Pt
                unknowns['flow_out:out:area'] = unknowns['Aexit_des']
            self._solve_flow_vars('flow_out', params, unknowns)
            if abs(unknowns['flow_ref:out:Ps'] - flow_out_supersonic.Ps) / unknowns['flow_ref:out:Ps'] < 0.001:
                unknowns['switchRegime'] = 'PERFECTLY_EXPANDED'
        unknowns['Fg'] = unknowns['flow_out:out:W'] * unknowns['flow_out:out:Vflow'] / 32.174 + unknowns['flow_out:out:area'] * (unknowns['flow_out:out:Ps'] - unknowns['flow_ref:out:Ps'])
        unknowns['PR'] = flow_throat.Pt / unknowns['flow_out:out:Ps']
        unknowns['AR'] = unknowns['flow_out:out:area'] / flow_throat.area
        if unknowns['switchRegime'] == 'UNCHOKED':
            unknowns['WqAexit'] = unknowns['flow_in:out:W'] / unknowns['flow_ref:out:Ps']
            unknowns['WqAexit_dmd'] = unknowns['flow_in:out:W'] / unknowns['flow_out:out:Ps']
        else:
            unknowns['WqAexit'] = unknowns['flow_in:out:W'] / unknowns['Athroat_des']
            unknowns['WqAexit_dmd'] = unknowns['flow_in:out:W'] / unknowns['Athroat_dmd']
