from os.path import dirname, join
from Cantera import *
from scipy.optimize import newton

import pycycle

GAS_CONSTANT = 0.0685592 # BTU/lbm-R
SET_BY_NONE, SET_BY_Mach, SET_BY_area, SET_BY_Ps = 0, 1, 2, 3 # in general, target whichever variable was changed last

REACTANT_NAMES = []
REACTANT_SPLITS =[]
# add a reactant that can be mixed in
def add_reactant(reactants, splits):
    assert len(REACTANT_NAMES) == len(REACTANT_SPLITS)
    assert len(reactants) == len(splits)
    REACTANT_NAMES.append(reactants)
    REACTANT_SPLITS.append(splits)
# For right now, all FlowStations are Air/Fuel FlowStations
add_reactant(['N2', 'O2', 'AR', 'CO2'], [0.755184, 0.231416, 0.012916, 0.000485])
add_reactant(['H2O'], [1.0])
add_reactant(['CH2', 'CH'], [0.922189, 0.07781])
add_reactant(['C', 'H'], [0.86144, 0.13856])
add_reactant(['Jet-A(g)'], [1.0])
add_reactant(['H2'], [1.0])

class FlowStationData:
    def __init__(self, station_name, flow, species):
        self.station_name    = station_name
        self.flow            = flow
        self.species         = species

def init_fs_tree(add_var, station_name):
    '''Adds FlowStation variables to a data structure (e.g. a component's parameters or unknowns) and returns an object containing additional data.'''
    def add_all(all_vars):
        for var_name, args in all_vars.iteritems():
            default_value, desc, units = args[0], args[1], args[2]
            add_var('%s:%s' % (station_name, var_name), default_value, desc=desc, units=units)
    add_all({
        'ht':       (0.0, 'total enthalpy', 'Btu/lbm'),
        'Tt':       (0.0, 'total temperature', 'degR'),
        'Pt':       (0.0, 'total pressure', 'lbf/inch**2'),
        'rhot':     (0.0, 'total density', 'lbm/ft**3'),
        'gamt':     (0.0, 'total gamma', None),
        'Cp':       (0.0, 'specific heat at constant pressure', 'Btu/lbm*degR'),
        'Cv':       (0.0, 'specific heat at constant voluem', 'Btu/lbm*degR'),
        's':        (0.0, 'entropy', 'Btu/(lbm*R)'),
        'W':        (0.0, 'weight flow', 'lbm/s'),
        'FAR':      (0.0, 'fuel to air ratio', None),
        'WAR':      (0.0, 'weight to air ratio', None),
        'hs':       (0.0, 'static enthalpy', 'Btu/lbm'),
        'Ts':       (0.0, 'static temperature', 'degR'),
        'Ps':       (0.0, 'static pressure', 'lbf/inch**2'),
        'rhos':     (0.0, 'static density', 'lbm/ft**3'),
        'gams':     (0.0, 'static gamma', None),
        'Vflow':    (0.0, 'velocity', 'ft/s'),
        'Vsonic':   (0.0, 'speed of sound', 'ft/s'),
        'Mach':     (0.0, 'mach number', None),
        'area':     (0.0, 'flow area', 'inch**2'),
        'is_super': (False, 'selects preference for supersonic versus subsonic solution when setting area', None),
        'Wc':       (0.0, 'corrected weight flow', 'lbm/s')
    })
    #properties file path
    directory = dirname(pycycle.__file__)
    prop_file = join(directory, 'gri1000.cti')
    flow = importPhase(prop_file)
    species = [1.0, 0, 0, 0, 0, 0, 0, 0]
    return FlowStationData(station_name, flow, species)

def init_fs_standalone(variables={}, station_name=''):
    '''Returns a tuple containing a dictionary of FlowStation variables and an object containing additional data.'''
    variables = {}
    def add_var(name, default_value, *args, **kwargs):
        variables[name] = default_value
    fs_data = init_fs_tree(add_var, '')
    return (variables, fs_data)

def _set_comp(species, flow):
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

def set_dry_air(variables, data):
    '''Set the composition to dry air'''
    sn = data.station_name
    data.species[0] = 1
    variables['%s:WAR' % sn] = 0
    variables['%s:FAR' % sn] = 0
    _set_comp(data.species, data.flow)
        
def set_reactant(data, i):
    '''Set the composition to pure mixture of one of the reactants'''
    data.species = [0 for x in range(len(REACTANT_NAMES))]
    data.species[i - 1] = 1
              
def set_WAR(variables, data, WAR, set_statics_by):
    '''Set the compositon to air with water'''
    sn = data.station_name
    variables['%s:WAR' % sn] = WAR
    varialbes['%s:FAR' % sn] = 0
    data.species[0] = 1 / (1 + WAR)
    data.species[1] = WAR / (1 + WAR)
    _set_comp(data.species, data.flow)
    set_static(variables, data, set_statics_by)
        
def _total_calcs(variables, data, set_statics_by):
    sn = data.station_name
    variables['%s:ht' % sn]     = data.flow.enthalpy_mass() * 0.0004302099943161011
    variables['%s:s' % sn]      = data.flow.entropy_mass() * 0.000238845896627
    variables['%s:rhot' % sn]   = data.flow.density() * 0.0624
    variables['%s:Tt' % sn]     = data.flow.temperature() * 9.0 / 5.0
    variables['%s:Cp' % sn]     = data.flow.cp_mass() * 2.388459e-4
    variables['%s:Cv' % sn]     = data.flow.cv_mass() * 2.388459e-4
    variables['%s:gamt' % sn]   = variables['%s:Cp' % sn] / variables['%s:Cv' % sn]
    set_static(variables, data, set_statics_by)
    variables['%s:Vsonic' % sn] = math.sqrt(variables['%s:gams' % sn] * GasConstant * data.flow.temperature() / data.flow.meanMolecularWeight()) * 3.28084

def set_total_TP(variables, data, Tin, Pin, set_statics_by):
    '''Set total conditions based on T an P'''
    sn = data.station_name
    _set_comp(data.species, data.flow)
    variables['%s:Tt' % sn] = Tin
    variables['%s:Pt' % sn] = Pin
    data.flow.set(T=Tin * 5.0 / 9.0, P=Pin * 6894.75729)
    data.flow.equilibrate('TP')
    _total_calcs(variables, data, set_statics_by)

def set_total_hP(variables, data, hin, Pin, set_statics_by):
    '''Set total conditions based on h and P'''
    sn = data.station_name
    _set_comp(data.species, data.flow)
    variables['%s:ht' % sn] = hin
    variables['%s:Pt' % sn] = Pin
    data.flow.set(H=hin / 0.0004302099943161011, P=Pin * 6894.75729)
    data.flow.equilibrate('HP', loglevel=1)
    _total_calcs(variables, data, set_statics_by)

def set_total_sP(variables, data, Sin, Pin, set_statics_by):
    '''Set total condition based on S and P'''
    sn = data.station_name
    _set_comp(data.species, data.flow)
    variables['%s:s' % sn] = Sin
    variables['%s:Pt' % sn] = Pin
    data.flow.set(S=Sin / 0.000238845896627, P=Pin * 6894.75729)
    data.flow.equilibrate('SP', loglevel=1)
    _total_calcs(variables, data, set_statics_by)

def add(variables, data, variables2, data2):
    '''Add another station to this one; mix enthalpies and keep pressure and this stations value'''
    sn1, sn2 = data.station_name, data2.station_name
    vars1, vars2 = variables, variables2
    temp = ''
    for i in range(0, len(data._species)):
            data.species[i] = (vars1['%s:W' % sn1] * data.species[i] + vars2['%s:W' % sn2] * data2.species[i]) / (vars1['%s:W' % sn1] + vars2['%s:W' % sn2])
    _set_comp(data.species, data.flow)
    air1 = vars1['%s:W' % sn1] * (1.0 / (1.0 + vars1['%s:FAR' % sn1] + vars1['%s:WAR' % sn1]))
    air2 = vars2['%s:W' % sn2] * (1.0 / (1.0 + vars2['%s:FAR' % sn2] + vars2['%s:WAR' % sn2]))
    vars1['%s:FAR' % sn1] = (air1 * vars1['%s:FAR' % sn1] + air2 * vars2['%s:FAR' % sn2]) / (air1 + air2)
    vars1['%s:WAR' % sn1] = (air1 * vars1['%s:WAR' % sn1] + air2 * vars2['%s:WAR' % sn2]) / (air1 + air2)
    vars1['%s:ht' % sn1] = (vars1['%s:W' % sn1] * vars1['%s:ht' % sn1] + vars2['%s:W' % sn2] * vars2['%s:ht' % sn2]) / (vars1['%s:W' % sn1] + vars2['%s:W' % sn2])
    vars1['%s:W' % sn1] = vars1['%s:W' % sn1] + vars2['%s:W' % sn2]
    data.flow.set(T=vars1['%s:Tt' % sn1] * 5.0 / 9.0, P=vars1['%s:Pt' % sn1] * 6894.75729)
    data.flow.equilibrate('TP')
    data.flow.set(H=vars1['%s:ht' % sn1] / 0.0004302099943161011, P=vars1['%s:Pt' % sn1] * 6894.75729)
    data.flow.equilibrate('HP')
    vars1['%s:Tt' % sn1] = data.flow.temperature() * 9.0 / 5.0
    vars1['%s:s' % sn1] = data.flow.entropy_mass() * 0.000238845896627
    vars1['%s:rhot' % sn1] = data.flow.density() * 0.0624
    vars1['%s:gamt' % sn1] = data.flow.cp_mass() / data.flow.cv_mass()

def copy_from(variables, data, variables2, data2):
    '''Duplicates total properties from another flow station'''
    sn1, sn2 = data.station_name, data2.station_name
    vars1, vars2 = variables, variables2
    data.species = data2.species[:]
    vars1['%s:ht' % sn1]   = vars2['%s:ht' % sn2]
    vars1['%s:Tt' % sn1]   = vars2['%s:Tt' % sn2]
    vars1['%s:Pt' % sn1]   = vars2['%s:Pt' % sn2]
    vars1['%s:rhot' % sn1] = vars2['%s:rhot' % sn2]
    vars1['%s:gamt' % sn1] = vars2['%s:gamt' % sn2]
    vars1['%s:s' % sn1]    = vars2['%s:s' % sn2]
    vars1['%s:W' % sn1]    = vars2['%s:W' % sn2]
    vars1['%s:FAR' % sn1]  = vars2['%s:FAR' % sn2]
    vars1['%s:WAR' % sn1]  = vars2['%s:WAR' % sn2]
    temp =''
    _set_comp(data.species, data.flow)
    data.flow.set(T=vars1['%s:Tt' % sn1] * 5.0 / 9.0, P=vars1['%s:Pt' % sn1] * 6894.75729)
    data.flow.equilibrate('TP')
                    
def burn(variables, data, fuel, Wfuel, hfuel):
    '''Burn a fuel with this station'''
    sn = data.station_name
    flow_1 = variables['%s:W' % sn]
    variables['%s:W' % sn] = variables['%s:W' % sn] + Wfuel 
    for i in range(0, len(data.species)):
        if (fuel - 1) == i:
            data.species[i] = (flow_1 * data.species[i] + Wfuel) / variables['%s:W' % sn]
        else:
            data.species[i] = (flow_1 * data.species[i]) / variables['%s:W' % sn]
    variables['%s:ht' % sn] = (flow_1 * variables['%s:ht' % sn] + Wfuel * hfuel)/ variables['%s:W' % sn]
    air1 = flow_1 * (1.0 / (1.0 + variables['%s:FAR' % sn] + variables['%s:WAR' % sn]))
    variables['%s:FAR' % sn] = (air1 * variables['%s:FAR' % sn] + Wfuel) / air1
    _set_comp(data.species, data.flow)
    data.flow.set(T=2660 * 5 / 9, P=variables['%s:Pt' % sn] * 6894.75729)
    data.flow.equilibrate('TP')
    data.flow.set(H=variables['%s:ht' % sn] / 0.0004302099943161011, P=variables['%s:Pt' % sn] * 6894.75729)
    data.flow.equilibrate('HP')
    variables['%s:Tt' % sn] = data.flow.temperature() * 9.0 / 5.0
    variables['%s:s' % sn] = data.flow.entropy_mass() * 0.000238845896627
    variables['%s:rhot' % sn] = data.flow.density() * 0.0624
    variables['%s:gamt' % sn] = data.flow.cp_mass() / data.flow.cv_mass()

def set_static_Mach(variables, data):
    '''Set the statics based on Mach'''
    sn = data.station_name
    mach_target = variables['%s:Mach' % sn]
    def f(Ps):
        variables['%s:Ps' % sn] = Ps
        set_static_Ps(variables, data)
        print variables['%s:Mach' % sn], mach_target
        return variables['%s:Mach' % sn] - mach_target
    Ps_guess = variables['%s:Pt' % sn] * (1 + (variables['%s:gamt' % sn] - 1) / 2 * mach_target ** 2) ** (variables['%s:gamt' % sn] / (1 - variables['%s:gamt' % sn])) * 0.9
    newton(f, Ps_guess)
    print variables['%s:Mach' % sn]

def set_static_Ps(variables, data):
    '''Set the statics based on pressure'''
    sn = data.station_name
    _set_comp(data.species, data.flow)
    def f(Ts):
        data.flow.set(T=Ts * 5.0 / 9.0, P=variables['%s:Ps' % sn] * 6894.75729) # 6894.75729 Pa/psi
        data.flow.equilibrate('TP')
        return variables['%s:s' % sn] - data.flow.entropy_mass() * 0.000238845896627 # 0.0002... kCal/N-m
    newton(f, variables['%s:Ts' % sn])
    variables['%s:Ts' % sn] = data.flow.temperature() * 9.0 / 5.0
    variables['%s:rhos' % sn] = data.flow.density() * 0.0624
    variables['%s:gams' % sn] = data.flow.cp_mass() / data.flow.cv_mass()
    variables['%s:hs' % sn] = data.flow.enthalpy_mass() * 0.0004302099943161011
    variables['%s:Vflow' % sn] = math.sqrt((778.169 * 32.1740 * 2 * (variables['%s:ht' % sn] - variables['%s:hs' % sn]))) # 778.169 lb-f / J; 32.1740 ft/s^2 = g
    variables['%s:Vsonic' % sn] = math.sqrt(variables['%s:gams' % sn] * GasConstant * data.flow.temperature() / data.flow.meanMolecularWeight()) * 3.28084
    variables['%s:Mach' % sn] = variables['%s:Vflow' % sn] / variables['%s:Vsonic' % sn]
    variables['%s:area' % sn] = variables['%s:W' % sn] / (variables['%s:rhos' % sn] * variables['%s:Vflow' % sn]) * 144.0

def set_static_area(variables, data):
    '''Set the statics based on area'''
    sn = data.station_name
    target_area = variables['%s:area' % sn]
    Ps_guess = variables['%s:Pt' % sn] * (1 + (variables['%s:gamt' % sn] - 1) / 2) ** (variables['%s:gamt' % sn] / ( 1 - variables['%s:gamt' % sn])) # at mach 1
    def f(Ps):
        variables['%s:Ps' % sn] = Ps
        set_static_Ps(variables, data)
        return 1 - variables['%s:Mach' % sn]
    Ps_M1 = newton(f, Ps_guess)
    # find the subsonic solution first
    guess = (variables['%s:Pt' % sn] + Ps_M1) / 2
    def f(Ps):
        variables['%s:Ps' % sn] = Ps
        set_static_Ps(variables, data)
        return variables['%s:W' % sn] / (variables['%s:rhos' % sn] * variables['%s:Vflow' % sn]) * 144.0 - target_area
    newton(f, guess)
    # if you want the supersonic one, just keep going with a little lower initial guess    
    if variables['%s:is_super' % sn]:
        # jsg: wild guess of 1/M_subsonic
        mach_guess = 1 / variables['%s:Mach' % sn]
        Ps_guess = variables['%s:Pt' % sn] * (1 + (variables['%s:gamt' % sn] - 1) / 2 * mach_guess ** 2) ** (variables['%s:gamt' % sn] / (1 - variables['%s:gamt' % sn]))
        newton(f, Ps_guess)

def set_static(variables, data, set_by):
    '''Determine which static calc to use'''
    sn = data.station_name
    if variables['%s:Tt' % sn] and variables['%s:Pt' % sn]: # if non zero
        variables['%s:Wc' % sn] = math.sqrt(variables['%s:W' % sn] * (variables['%s:Tt' % sn] / 518.67)) / (variables['%s:Pt' % sn] / 14.696)
    if set_by == SET_BY_NONE:
        variables['%s:Ps' % sn]    = variables['%s:Pt' % sn]
        variables['%s:Ts' % sn]    = variables['%s:Tt' % sn]
        variables['%s:rhos' % sn]  = variables['%s:rhot' % sn]
        variables['%s:gams' % sn]  = variables['%s:gamt' % sn]
        variables['%s:hs' % sn]    = variables['%s:ht' % sn]
        variables['%s:Vflow' % sn] = 0.0
        variables['%s:Mach' % sn]  = 0.0
    elif set_by == SET_BY_Mach:
        set_static_Mach(variables, data)
    elif set_by == SET_BY_area:
        set_static_area(variables, data)
    elif set_by == SET_BY_Ps:
        set_static_Ps(variables, data)

def set_static_Ts_Ps_MN(variables, data, Ts, Ps, MN): 
    '''Set the statics based on Ts, Ps, and MN'''
    # UPDGRADED TO USE LOOPS
    sn = data.station_name
    # do this twice beacause gamt changes
    for n in range(2):
        variables['%s:Tt' % sn] = variables['%s:Ts' % sn] * (1 + (variables['%s:gamt' % sn] - 1) / 2.0 * MN ** 2)
        variables['%s:Pt' % sn] = Ps * (1 + (variables['%s:gamt' % sn] - 1) / 2.0 * MN ** 2) ** (variables['%s:gamt' % sn] / (variables['%s:gamt' % sn] - 1))
        set_total_TP(variables, data, variables['%s:Tt' % sn], variables['%s:Pt' % sn])
    variables['%s:Mach' % sn] = MN
    set_static_Mach(variables, data)
    variables['%s:area' % sn] = variables['%s:W' % sn] / (variables['%s:rhos' % sn] * variables['%s:Vflow' % sn]) * 144.0
