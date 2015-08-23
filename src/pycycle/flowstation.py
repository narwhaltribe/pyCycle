from os.path import dirname, join
from Cantera import *
from scipy.optimize import newton
from openmdao.core.component import Component

import pycycle

GAS_CONSTANT = 0.0685592 # Btu/lbm-R

class FlowStation(Component):
    '''Thermodynamics analysis of fluid flow'''
    reactant_names = []
    reactant_splits =[]
    # add a reactant that can be mixed in
    @staticmethod
    def add_reactant(reactants, splits):
        assert len(reactant_names) == len(reactant_splits)
        assert len(reactants) == len(splits)
        reactant_names.append(reactants)
        reactant_splits.append(splits)
    # For right now, all FlowStations are Air/Fuel FlowStations
    add_reactant(['N2', 'O2', 'AR', 'CO2'], [0.755184, 0.231416, 0.012916, 0.000485])
    add_reactant(['H2O'], [1.0])
    add_reactant(['CH2', 'CH'], [0.922189, 0.07781])
    add_reactant(['C', 'H'], [0.86144, 0.13856])
    add_reactant(['Jet-A(g)'], [1.0])
    add_reactant(['H2'], [1.0])

    _directory = dirname(pycycle.__file__)
    _prop_file = join(directory, 'gri1000.cti')

    def __init__(self):
        def add_duo(name, desc, units):
            self.add_param(name, 0.0, desc=desc, units=units)
            self.add_output(name, 0.0, desc=desc, units=units)
        add_duo('ht', 'total enthalpy', 'Btu/lbm')
        add_duo('Tt', 'total temperature', 'degR')
        add_duo('Pt', 'total pressure', 'lbf/inch**2')
        add_duo('rhot', 'total density', 'lbm/ft**3')
        add_duo('s', 'entropy', 'Btu/(lbm*R)')
        add_duo('hs', 'static enthalpy', 'Btu/lbm')
        add_duo('Ts', 'static temperature', 'degR')
        add_duo('Ps', 'static pressure', 'lbf/inch**2')
        add_duo('rhos', 'static density', 'lbm/ft**3')
        add_duo('Vflow', 'velocity', 'ft/s')
        add_duo('Vsonic', 'speed of sound', 'ft/s')
        add_duo('Mach', 'Mach number', None)
        add_duo('area', 'flow area', 'inch**2')
        self.add_param('W', 0.0, desc='weight flow', units='lbm/s')
        self.add_param('FAR', 0.0, desc='fuel-to-air ratio')
        self.add_param('WAR', 0.0, desc='water-to-air ratio')
        self.add_param('is_super', False, desc='selects preference for supersonic versus subsonic solution when setting area')
        self.add_output('gams', 0.0, 'static gamma')
        self.add_output('gamt', 0.0, 'total gamma')
        self.add_output('Cp', 0.0, desc='specific heat at constant pressure', units='Btu/lbm*degR')
        self.add_output('Cv', 0.0, desc='specific heat at constant volume', units='Btu/lbm*degR')
        self.add_output('Wc', 0.0, desc='corrected weight flow', units='lbm/s')
        #properties file path
        self._flow = importPhase(_prop_file)
        self._species = [1.0, 0, 0, 0, 0, 0, 0, 0]
        self._set_comp()

    def _set_comp(self):
        comp_names = []
        fracts = []
        assert len(reactant_names) == len(reactant_splits)
        for i in range(len(reactant_names)):
            assert len(reactant_names[i]) == len(reactant_splits[i])
            for j in range(len(reactant_names[i])):
                name = reactant_names[i][j]
                fract = reactant_splits[i][j] * _species[i]
                if fract > 1e-5:
                    if name in comp_names:
                        fracts[comp_names.index(name)] += fract
                    else:
                        comp_names.append(name)
                        fracts.append(fract)
        temp_comp = ''
        for i in range(len(comp_names)):
            temp_comp += '%s:%s ' % (comp_names[i], fracts[i])
        self._flow.setMassFractions(temp_comp)

    def set_dry_air(self, params):
        '''Set the composition to dry air'''
        self._species[0] = 1.0
        params['WAR'] = 0.0
        params['FAR'] = 0.0
        self._set_comp()
        
    def set_reactant(self, i):
        '''Set the composition to pure mixture of one of the reactants'''
        self._species = [0 for x in range(len(reactant_names))]
        self._species[i - 1] = 1
              
    def set_WAR(self, params, unknowns):
        '''Set the compositon to air with water'''
        params['FAR'] = 0.0
        data._species[0] = 1.0 / (1.0 + params['WAR'])
        data._species[1] = params['WAR'] / (1.0 + params['WAR'])
        self._set_comp()
        self.solve_statics(params, unknowns)
        
    def _total_calcs(self, unknowns):
        unknowns['ht']     = self._flow.enthalpy_mass() * 0.0004302099943161011
        unknowns['s']      = self._flow.entropy_mass() * 0.000238845896627
        unknowns['rhot']   = self._flow.density() * 0.0624
        unknowns['Tt']     = self._flow.temperature() * 9.0 / 5.0
        unknowns['Cp']     = self._flow.cp_mass() * 2.388459e-4
        unknowns['Cv']     = self._flow.cv_mass() * 2.388459e-4
        unknowns['gamt']   = unknowns['Cp'] / unknowns['Cv']
        self.solve_statics(params, unknowns)
        unknowns['Vsonic'] = math.sqrt(unknowns['gams'] * GasConstant * self._flow.temperature() / self._flow.meanMolecularWeight()) * 3.28084

    def solve_total_TP(self, params, unknowns):
        '''Set total conditions based on T an P'''
        self._flow.set(T=params['Tt'] * 5.0 / 9.0, P=params['Pt'] * 6894.75729)
        self._flow.equilibrate('TP')
        self._total_calcs(unknowns)

    def solve_total_hP(self, params, unknowns):
        '''Set total conditions based on h and P'''
        self._flow.set(H=params['ht'] / 0.0004302099943161011, P=params['Pt'] * 6894.75729)
        self._flow.equilibrate('HP')
        self._total_calcs(unknowns)

    def solve_total_sP(self, params, unknowns):
        '''Set total condition based on S and P'''
        self._flow.set(S=params['s'] / 0.000238845896627, P=params['Pt'] * 6894.75729)
        self._flow.equilibrate('SP')
        self._total_calcs(unknowns)

# TODO implement self.add() and self.copy_from(), if appropriate

    def burn(self, params, unknowns, fuel, Wfuel, hfuel):
        '''Burn a fuel with this station'''
        W_1 = unknowns['W'] if unknowns['W'] else params['W']
        unknowns['W'] = W_1 + Wfuel 
        for i in range(0, len(self._species)):
            if (fuel - 1) == i:
                self._species[i] = (W_1 * self._species[i] + Wfuel) / unknowns['W']
            else:
                self._species[i] = (W_1 * self._species[i]) / unknowns['W']
        unknowns['ht'] = (W_1 * (unknowns['ht'] if unknowns['ht'] else params['ht']) + Wfuel * hfuel) / unknowns['W']
        FAR = unknowns['FAR'] if unknowns['FAR'] else params['FAR']
        air1 = W_1 * (1.0 / (1.0 + FAR + params['WAR']))
        unknowns['FAR'] = (air1 * FAR + Wfuel) / air1
        self._set_comp(unknowns)
        self._flow.set(T=2660 * 5 / 9, P=params['Pt'] * 6894.75729)
        self._flow.equilibrate('TP')
        self._flow.set(H=unknowns['ht'] / 0.0004302099943161011, P=params['Pt'] * 6894.75729)
        data._flow.equilibrate('HP')
        unknowns['Tt']   = self._flow.temperature() * 9.0 / 5.0
        unknowns['s']    = self._flow.entropy_mass() * 0.000238845896627
        unknowns['rhot'] = self._flow.density() * 0.0624
        unknowns['gamt'] = self._flow.cp_mass() / self._flow.cv_mass()

    def solve_statics_Mach(self, params, unknowns, Mach=None):
        '''Set the statics based on Mach'''
        Mach = Mach if Mach else params['Mach']
        def f(Ps):
            unknowns['Ps'] = Ps
            self.solve_statics_Ps(self, params, unknowns, Ps=Ps)
            return unknowns['Mach'] - Mach
        Ps_guess = params['Pt'] * (1.0 + (unknowns['gamt'] - 1.0) / 2.0 * mach_target ** 2) ** (unknowns['gamt'] / (1.0 - unknowns['gamt'])) * 0.9
        newton(f, Ps_guess)

    def solve_statics_Ps(self, params, unknowns, Ps=None):
        '''Set the statics based on pressure'''
        Ps = Ps if Ps else params['Ps']
        s = unknowns['s'] if unknowns['s'] else params['s']
        Ts = unknowns['Ts'] if unknowns['Ts'] else params['Ts']
        ht = unknowns['ht'] if unknowns['ht'] else params['ht']
        rhos = unknowns['rhos'] if unknowns['rhos'] else params['rhos']
        def f(Ts_tmp):
            self._flow.set(T=Ts_tmp * 5.0 / 9.0, P=Ps * 6894.75729) # 6894.75729 Pa/psi
            self._flow.equilibrate('TP')
            return s - self._flow.entropy_mass() * 0.000238845896627 # 0.0002... kCal/N-m
        newton(f, Ts)
        unknowns['Ts']     = self._flow.temperature() * 9.0 / 5.0
        unknowns['rhos']   = self._flow.density() * 0.0624
        unknowns['gams']   = self._flow.cp_mass() / self._flow.cv_mass()
        unknowns['hs']     = self._flow.enthalpy_mass() * 0.0004302099943161011
        unknowns['Vflow']  = math.sqrt((778.169 * 32.1740 * 2 * (ht - unknowns['hs']))) # 778.169 lb-f / J; 32.1740 ft/s^2 = g
        unknowns['Vsonic'] = math.sqrt(unknowns['gams'] * GasConstant * self._flow.temperature() / self._flow.meanMolecularWeight()) * 3.28084
        unknowns['Mach']   = unknowns['Vflow'] / unknowns['Vsonic']
        unknowns['area']   = params['W'] / (rhos * unknowns['Vflow']) * 144.0

    def solve_statics_area(self, params, unknowns):
        '''Set the statics based on area'''
        Pt = unknowns['Pt'] if unknowns['Pt'] else params['Pt']
        Ps_guess = Pt * (1.0 + (unknowns['gamt'] - 1.0) / 2.0) ** (unknowns['gamt'] / (1.0 - unknowns['gamt'])) # at mach 1
        self.solve_statics_Mach(params, unknowns, Mach=1.0)
        Ps_M1 = unknowns['Ps']
        # find the subsonic solution first
        guess = (Pt + Ps_M1) / 2.0
        def f(Ps):
            unknowns['Ps'] = Ps
            self.solve_statics_Ps(params, unknowns, Ps=Ps)
            return params['area'] - params['W'] / (unknowns['rhos'] * unknowns['Vflow']) * 144.0
        newton(f, guess)
        # if you want the supersonic one, just keep going with a little lower initial guess    
        if params['is_super']:
            # jsg: wild guess of 1/M_subsonic
            mach_guess = 1.0 / unknowns['Mach']
            Ps_guess = Pt * (1.0 + (unknowns['gamt'] - 1.0) / 2.0 * mach_guess ** 2) ** (unknowns['gamt'] / (1.0 - unknowns['gamt']))
            newton(f, Ps_guess)

    def solve_statics(self, params, unknowns):
        '''Determine which static calc to use'''
        Tt = unknowns['Tt'] if unknowns['Tt'] else params['Tt']
        if Tt and params['Pt']: # if non zero
            W = unknowns['W'] if unknowns['W'] else params['W']
            unknowns['Wc'] = math.sqrt(W * (Tt / 518.67)) / (params['Pt'] / 14.696)
        if params['Mach']:
            self.solve_statics_Mach(params, unknowns)
        elif params['area']:
            self.solve_statics_area(params, unknowns)
        elif params['Ps']:
            self.solve_statics_Ps(params, unknowns)
        else:
            unknowns['Ps']    = params['Pt']
            unknowns['Ts']    = unknowns['Tt'] if unknowns['Tt'] else params['Tt']
            unknowns['rhos']  = unknowns['rhot'] if unknowns['rhot'] else params['rhot']
            unknowns['gams']  = unknowns['gamt'] if unknowns['gamt'] else params['gamt']
            unknowns['hs']    = unknowns['ht'] if unknowns['ht'] else params['ht']
            unknowns['Vflow'] = 0.0
            unknowns['Mach']  = 0.0

    def solve_statics_Ts_Ps_MN(self, params, unknowns, Ts, Ps, MN):
        '''Set the statics based on Ts, Ps, and MN'''
        # UPDGRADED TO USE LOOPS
        # do this twice beacause gamt changes
        for n in range(2):
            Ts_tmp = unknowns['Ts'] if unknowns['Ts'] else params['Ts']
            gamt = unknowns['gamt'] if unknowns['gamt'] else params['gamt']
            unknowns['Tt'] = Ts_tmp * (1.0 + (gamt - 1.0) / 2.0 * MN ** 2)
            unknowns['Pt'] = Ps * (1.0 + (gamt - 1.0) / 2.0 * MN ** 2) ** (gamt / gamt - 1.0))
            set_total_TP(variables, data, params['Tt'], params['Pt'])
        solve_statics_Mach(params, unknowns)
        unknowns['area'] = params['W'] / (unknowns['rhos'] * unknowns['Vflow']) * 144.0
