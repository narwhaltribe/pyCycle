from os.path import dirname, join
from collections import namedtuple
from Cantera import *
from scipy.optimize import newton
from openmdao.core.component import Component

import pycycle

GAS_CONSTANT = 0.0685592 # Btu/lbm-R

Output = namedtuple('Output', ['ht', 'Tt', 'Pt', 's', 'hs', 'Ts', 'Ps', 'Mach', 'area', 'Vsonic', 'Vflow', 'rhos', 'rhot', 'gams', 'gamt', 'Cp', 'Cv', 'Wc'])

REACTANT_NAMES = []
REACTANT_SPLITS =[]
def add_reactant(reactants, splits):
    '''Add a reactant that can be mixed in'''
    assert len(REACTANT_NAMES) == len(REACTANT_SPLITS)
    assert len(reactants) == len(splits)
    REACTANT_NAMES.append(reactants)
    REACTANT_SPLITS.append(splits)
# For right now, all are Air/Fuel 
add_reactant(['N2', 'O2', 'AR', 'CO2'], [0.755184, 0.231416, 0.012916, 0.000485])
add_reactant(['H2O'], [1.0])
add_reactant(['CH2', 'CH'], [0.922189, 0.07781])
add_reactant(['C', 'H'], [0.86144, 0.13856])
add_reactant(['Jet-A(g)'], [1.0])
add_reactant(['H2'], [1.0])

_DIRECTORY = dirname(pycycle.__file__)
_PROP_FILE = join(_DIRECTORY, 'gri1000.cti')
def _init_flow(species=[1.0, 0, 0, 0, 0, 0, 0, 0]):
    flow = importPhase(_PROP_FILE)
    _set_comp(flow, species)
    return flow

def _set_comp(flow, species):
    comp_names = []
    fracts = []
    assert len(REACTANT_NAMES) == len(REACTANT_SPLITS)
    for i in range(len(REACTANT_NAMES)):
        assert len(REACTANT_NAMES[i]) == len(REACTANT_SPLITS[i])
        for j in range(len(REACTANT_NAMES[i])):
            name = REACTANT_NAMES[i][j]
            fract = REACTANT_SPLITS[i][j] * species[i]
            if fract > 1e-5:
                if name in comp_names:
                    fracts[comp_names.index(name)] += fract
                else:
                    comp_names.append(name)
                    fracts.append(fract)
    temp_comp = ''
    for i in range(len(comp_names)):
        temp_comp += '%s:%s ' % (comp_names[i], fracts[i])
    flow.setMassFractions(temp_comp)

#    def set_dry_air(self, params=None):
#        '''Set the composition to dry air'''
#        params = params if params is not None else self.params
#        self._species[0] = 1.0
#        params['WAR'] = 0.0
#        params['FAR'] = 0.0
#        self._set_comp()
#        
#    def set_reactant(self, i):
#        '''Set the composition to pure mixture of one of the reactants'''
#        self._species = [0 for x in range(len(REACTANT_NAMES))]
#        self._species[i - 1] = 1
#              
#    def set_WAR(self, params=None, unknowns=None):
#        '''Set the compositon to air with water'''
#        params = params if params is not None else self.params
#        unknowns = unknowns if unknowns is not None else self.unknowns
#        params['FAR'] = 0.0
#        data._species[0] = 1.0 / (1.0 + params['WAR'])
#        data._species[1] = params['WAR'] / (1.0 + params['WAR'])
#        self._set_comp()
#        self.solve_statics(params, unknowns)
        
def _total_calcs(flow, W, hs, Ts, Ps, Mach, area, is_super):
    ht = flow.enthalpy_mass() * 0.0004302099943161011
    s = flow.entropy_mass() * 0.000238845896627
    rhot = flow.density() * 0.0624
    Tt = flow.temperature() * 9.0 / 5.0
    Cp = flow.cp_mass() * 2.388459e-4
    Cv = flow.cv_mass() * 2.388459e-4
    Pt = flow.P / 6894.75729
    gamt = Cp / Cv
    out = solve_statics(flow, W=W, hs=hs, Ts=Ts, Ps=Ps, Mach=Mach, area=area, is_super=is_super, ht=ht, Pt=Pt, s=s, rhot=rhot, Tt=Tt, Cp=Cp, Cv=Cv, gamt=gamt)
    Vsonic = math.sqrt(out.gams * GasConstant * flow.temperature() / flow.meanMolecularWeight()) * 3.28084
    return Output(ht=ht, Tt=Tt, Pt=Pt, s=s, hs=out.hs, Ts=out.Ts, Ps=out.Ps, Mach=out.Mach, area=out.area, Vsonic=Vsonic, Vflow=out.Vflow, rhos=out.rhos, rhot=rhot, gams=out.gams, gamt=gamt, Cp=Cp, Cv=Cv, Wc=out.Wc) 

def solve_totals(Pt, Tt=-1.0, ht=-1.0, s=-1.0, W=-1.0, hs=-1.0, Ts=-1.0, Ps=-1.0, Mach=-1.0, area=-1.0, is_super=False):
    '''Calculate total conditions based on T, h, or s, and P'''
    assert Tt != -1 or ht != -1 or s != -1
    flow = _init_flow()
    if Tt != -1:
        flow.set(T=Tt * 5.0 / 9.0, P=Pt * 6894.75729)
        flow.equilibrate('TP')
    elif ht != -1:
        flow.set(H=ht / 0.0004302099943161011, P=Pt * 6894.75729)
        flow.equilibrate('HP')
    else:
        flow.set(S=s / 0.000238845896627, P=Pt * 6894.75729)
        flow.equilibrate('SP')
    return _total_calcs(flow, W=W, hs=hs, Ts=Ts, Ps=Ps, Mach=Mach, area=area, is_super=is_super)

# TODO implement self.burn()
#    def burn(self, params, unknowns, fuel, Wfuel, hfuel):
#        '''Burn a fuel with this station'''
#        W_1 = params['W']
#        params['W'] = W_1 + Wfuel 
#        for i in range(0, len(self._species)):
#            if (fuel - 1) == i:
#                self._species[i] = (W_1 * self._species[i] + Wfuel) / params['W']
#            else:
#                self._species[i] = (W_1 * self._species[i]) / params['W']
#        unknowns['ht'] = (W_1 * (unknowns['ht'] if unknowns['ht'] else params['ht']) + Wfuel * hfuel) / params['W']
#        air1 = W_1 * (1.0 / (1.0 + params['FAR'] + params['WAR']))
#        params['FAR'] = (air1 * FAR + Wfuel) / air1
#        self._set_comp(unknowns)
#        self._flow.set(T=2660 * 5 / 9, P=params['Pt'] * 6894.75729)
#        self._flow.equilibrate('TP')
#        self._flow.set(H=unknowns['ht'] / 0.0004302099943161011, P=params['Pt'] * 6894.75729)
#        data._flow.equilibrate('HP')
#        unknowns['Tt']   = self._flow.temperature() * 9.0 / 5.0
#        unknowns['s']    = self._flow.entropy_mass() * 0.000238845896627
#        unknowns['rhot'] = self._flow.density() * 0.0624
#        unknowns['gamt'] = self._flow.cp_mass() / self._flow.cv_mass()

def solve_statics_Mach(Mach, Pt, gamt, Ts, ht, s, Tt, W):
    '''Calculate the statics based on Mach'''
    out = [None] # Makes out[0] a reference
    def f(Ps):
        out[0] = solve_statics_Ps(Ps=Ps, s=s, Ts=Ts, Tt=Tt, ht=ht, W=W)
        return out[0].Mach - Mach
    Ps_guess = Pt * (1.0 + (gamt - 1.0) / 2.0 * Mach ** 2) ** (gamt / (1.0 - gamt])) * 0.9
    newton(f, Ps_guess)
    return out[0]

def solve_statics_Ps(Ps, s, Ts, Tt, ht, W):
    '''Calculate the statics based on pressure'''
    flow = _init_flow()
    def f(Ts):
        flow.set(T=Ts * 5.0 / 9.0, P=Ps * 6894.75729) # 6894.75729 Pa/psi
        flow.equilibrate('TP')
        return s - flow.entropy_mass() * 0.000238845896627 # 0.0002... kCal/N-m
    Ts_guess = Ts if Ts != -1 else Tt
    newton(f, Ts_guess)
    Ts = flow.temperature() * 9.0 / 5.0
    rhos = flow.density() * 0.0624
    gams = flow.cp_mass() / flow.cv_mass()
    hs = flow.enthalpy_mass() * 0.0004302099943161011
    Vflow = math.sqrt((778.169 * 32.1740 * 2 * (ht - hs))) # 778.169 lb-f / J; 32.1740 ft/s^2 = g
    Vsonic = math.sqrt(gams * GasConstant * flow.temperature() / flow.meanMolecularWeight()) * 3.28084
    Mach = Vflow / Vsonic
    area = W / (rhos * Vflow) * 144.0
    return Output(ht=ht, Tt=Tt, Ps=Ps, s=s, Ts=Ts, W=W, rhos=rhos, gams=gams, hs=hs, Vflow=Vflow, Vsonic=Vsonic, Mach=Mach, area=area)

def solve_statics_area(area, Pt, gamt, Ts, ht, s, Tt, W, is_super):
    '''Calculate the statics based on area'''
    flow = _init_flow()
    Ps_guess = Pt * (1.0 + (gamt - 1.0) / 2.0) ** (gamt / (1.0 - gamt)) # at mach 1
    statics_M1 = solve_statics_Mach(Mach=1.0, Pt=Pt, gamt=gamt, Ts=Ts, ht=ht, s=s, Tt=Tt, W=W)
    # find the subsonic solution first
    guess = (Pt + statics_M1.Ps) / 2.0
    out = [None] # Makes out[0] a reference
    def f(Ps):
        out[0] = solve_statics_Ps(Ps=Ps, s=s, Ts=Ts, Tt=Tt, ht=ht, W=W)
        return area - W / (out[0].rhos * out[0].Vflow) * 144.0
    newton(f, guess)
    # if you want the supersonic one, just keep going with a little lower initial guess    
    if is_super:
        # jsg: wild guess of 1/M_subsonic
        Mach_guess = 1.0 / out[0].Mach
        Ps_guess = Pt * (1.0 + (gamt - 1.0) / 2.0 * Mach_guess ** 2) ** (gamt / (1.0 - gamt))
        newton(f, Ps_guess)
    return out[0]

def solve_statics(Tt=-1.0, Pt=-1.0, W=-1.0, Mach=-1.0, area=-1.0, Ps=-1.0, gamt=-1.0, rhot=-1.0, Ts=-1.0, ht=-1.0, s=-1.0, W=-1.0, is_super=False):
    '''Determine which static calc to use'''
    if Tt > 0 and Pt > 0: # if non zero
        Wc = math.sqrt(W * (Tt / 518.67)) / (Pt / 14.696)
    else:
        Wc = -1.0
    if Mach != -1:
        out = solve_statics_Mach(Mach, Pt=Pt, gamt=gamt, Ts=Ts, Tt=Tt, ht=ht, s=s, W=W)
    elif area != -1:
        out = solve_statics_area(area, Pt=Pt, gamt=gamt, Ts=Ts, ht=ht, s=s, Tt=Tt, W=W, is_super=is_super)
    elif Ps != -1:
        out = solve_statics_Ps(Ps, s=s, Ts=Ts, Tt=Tt, ht=ht, W=W)
    else:
        return Output(Ps=Pt, Ts=Tt, rhos=rhot, gams=gamt, hs=ht, Vflow=0.0, Mach=0.0, Pt=Pt, Tt=Tt, W=W, area=area, gamt=gamt, rhot=rhot, ht=ht, s=s, Wc=Wc)
    return out._replace(Wc=Wc)

# TODO implement self.solve_statics_Ts_Ps_MN()
#    def solve_statics_Ts_Ps_MN(self, params, unknowns):
#        '''Set the statics based on Ts, Ps, and MN'''
#        # UPDGRADED TO USE LOOPS
#        # do this twice beacause gamt changes
#        for n in range(2):
#            unknowns['Tt:out'] = params['Ts:in'] * (1.0 + (unknowns['gamt'] - 1.0) / 2.0 * params['Mach:in'] ** 2)
#            unknowns['Pt:out'] = params['Ps:in'] * (1.0 + (unknowns['gamt'] - 1.0) / 2.0 * params['Mach:in'] ** 2) ** (unknowns['gamt'] / (unknowns['gamt'] - 1.0))
#            self.solve_totals_TP(params, unknowns)
#            set_total_TP(variables, data, unknowns['Tt:out'], unknowns['Pt:out'])
#        solve_statics_Mach(params, unknowns)
#        unknowns['area:out'] = params['W'] / (unknowns['rhos'] * unknowns['Vflow']) * 144.0
