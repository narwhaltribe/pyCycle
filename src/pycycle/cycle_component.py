import math 

from openmdao.core.component import Component
from pycycle import flowstation

ALL_PARAMS = ('is_super', 'ht', 'Tt', 'Pt', 's', 'hs', 'Ts', 'Ps', 'Mach', 'area', 'W', 'FAR', 'WAR')
ALL_OUTPUTS = ('is_super', 'ht', 'Tt', 'Pt', 's', 'hs', 'Ts', 'Ps', 'Mach', 'area', 'W', 'FAR', 'WAR', 'Vsonic', 'Vflow', 'rhos', 'rhot', 'gams', 'gamt', 'Cp', 'Cv', 'Wc')

class CycleComponent(Component): 

    def __init__(self): 
        super(CycleComponent, self).__init__()
        self.add_param('design', False, desc='flag to indicate that the calculations are design conditions')

    @staticmethod
    def copy_from(comp1, name1, comp2, name2):
        '''Copies parameters from FlowStation 1 to FlowStation 2'''
        for var_name in ALL_PARAMS:
            comp2.params['%s:in:%s' % (name2, var_name)] = comp1.params['%s:in:%s' % (name1, var_name)]

    def _clear_unknowns(self, name, unknowns, var_names=ALL_OUTPUTS):
        '''Reset all of a FlowStation's unknowns to empty (-1.0).'''
        for var_name in var_names:
            unknowns['%s:out:%s' % (name, var_name)] = -1.0

    def _solve_flow_vars(self, name, params, unknowns):
        '''Solve a FlowStation's unknowns based on variables specified as parameters.'''
        def var(var_name):
            param_name = '%s:in:%s' % (name, var_name)
            output_name = '%s:out:%s' % (name, var_name)
            assert output_name in unknowns.keys() or param_name in params.keys()
            if not output_name in unknowns.keys():
                return params[param_name]
            if not param_name in params.keys():
                return unknowns[output_name]
            return unknowns[output_name] if unknowns[output_name] != -1 else params[param_name]
        out = flowstation.solve(ht=var('ht'), Tt=var('Tt'), Pt=var('Pt'), s=var('s'), hs=var('hs'), Ts=var('Ts'), Ps=var('Ps'), Mach=var('Mach'), area=var('area'), W=var('W'), is_super=var('is_super'))
        def set_vars(var_dict):
            for var_name, val in var_dict.iteritems():
                output_name = '%s:out:%s' % (name, var_name)
                unknowns[output_name] = val
        set_vars({'ht': out.ht,
                  'Tt': out.Tt,
                  'Pt': out.Pt,
                  's':  out.s,
                  'hs': out.hs,
                  'Ts': out.Ts,
                  'Ps': out.Ps,
                  'Mach': out.Mach,
                  'area': out.area,
                  'W': var('W'),
                  'FAR': var('FAR'),
                  'WAR': var('WAR'),
                  'Vsonic': out.Vsonic,
                  'Vflow': out.Vflow,
                  'rhos': out.rhos,
                  'rhot': out.rhot,
                  'gams': out.gams,
                  'gamt': out.gamt,
                  'Cp': out.Cp,
                  'Cv': out.Cv,
                  'Wc': out.Wc,
                  'is_super': var('is_super')})
    
    def _add_flowstation(self, name):
        '''Add a variable tree representing a FlowStation. Parameters are stored as self.parameters['FLOWSTATION NAME:in:VARIABLE NAME'] and outputs are stored as self.unknowns['FLOWSTATION NAME:out:VARIABLE NAME'].'''
        def add_output(var_name, desc, units=None, default_value=-1.0):
            self.add_output('%s:out:%s' % (name, var_name), -1.0, desc=desc, units=units)
        def add_duo(var_name, desc, units=None, default_value=-1.0):
            self.add_param('%s:in:%s' % (name, var_name), -1.0, desc=desc, units=units)
            add_output(var_name, desc, units)
        add_duo('is_super', 'selects preference for supersonic versus subsonic solution when setting area', default_value=False)
        add_duo('ht', 'total enthalpy', 'Btu/lbm')
        add_duo('Tt', 'total temperature', 'degR')
        add_duo('Pt', 'total pressure', 'lbf/inch**2')
        add_duo('s', 'entropy', 'Btu/(lbm*R)')
        add_duo('hs', 'static enthalpy', 'Btu/lbm')
        add_duo('Ts', 'static temperature', 'degR')
        add_duo('Ps', 'static pressure', 'lbf/inch**2')
        add_duo('Mach', 'Mach number')
        add_duo('area', 'flow area', 'inch**2')
        add_duo('W', 'weight flow', 'lbm/s', 0.0)
        add_duo('FAR', 'fuel-to-air ratio')
        add_duo('WAR', 'water-to-air ratio')
        add_output('Vsonic', 'speed of sound', 'ft/s')
        add_output('Vflow', 'velocity', 'ft/s')
        add_output('rhos', 'static density', 'lbm/ft**3')
        add_output('rhot', 'total density', 'lbm/ft**3')
        add_output('gams', 'static gamma')
        add_output('gamt', 'total gamma')
        add_output('Cp', 'specific heat at constant pressure', 'Btu/lbm*degR')
        add_output('Cv', 'specific heat at constant volume', 'Btu/lbm*degR')
        add_output('Wc', 'corrected weight flow', 'lbm/s')
