import math 

from openmdao.core.component import Component
from pycycle import flowstation

FLOW_VARS = {'ht':       ('total enthalpy', -1.0, 'Btu/lbm'),
             'Tt':       ('total temperature', -1.0, 'degR'),
             'Pt':       ('total pressure', -1.0, 'lbf/inch**2'),
             's':        ('entropy', -1.0, 'Btu/(lbm*R)'),
             'hs':       ('static enthalpy', -1.0, 'Btu/lbm'),
             'Ts':       ('static temperature', -1.0, 'degR'),
             'Ps':       ('static pressure', -1.0, 'lbf/inch**2'),
             'Mach':     ('Mach number', -1.0, ''),
             'area':     ('flow area', -1.0, 'inch**2'),
             'W':        ('weight flow', 0.0, 'lbm/s'),
             'FAR':      ('fuel-to-air ratio', -1.0, ''),
             'WAR':      ('water-to-air ratio', -1.0, ''),
             'Vsonic':   ('speed of sound', -1.0, 'ft/s'),
             'Vflow':    ('velocity', -1.0, 'ft/s'),
             'rhos':     ('static density', -1.0, 'lbm/ft**3'),
             'rhot':     ('total density', -1.0, 'lbm/ft**3'),
             'gams':     ('static gamma', -1.0, ''),
             'gamt':     ('total gamma', -1.0, ''),
             'Cp':       ('specific heat at constant pressure', -1.0, 'Btu/lbm*degR'),
             'Cv':       ('specific heat at constant volume', -1.0, 'Btu/lbm*degR'),
             'Wc':       ('corrected weight flow', -1.0, 'lbm/s'),
             'is_super': ('selects preference for supersonic versus subsonic solution when setting area', False, '')}

class CycleComponent(Component): 
    def __init__(self): 
        super(CycleComponent, self).__init__()
        self.add_param('design', False, desc='flag to indicate that the calculations are design conditions')

    @staticmethod
    def connect_flows(group, flow1, flow2):
        '''Connects flow variable trees. Both flow1 and flow2 are strings (e.g. 'component.flow_name')'''
        for var_name in FLOW_VARS:
            group.connect('%s:out:%s' % (flow1, var_name), '%s:in:%s' % (flow2, var_name))

    @staticmethod
    def copy_from(comp1, name1, comp2, name2):
        '''Copies parameters from FlowStation 1 to FlowStation 2'''
        for var_name in FLOW_VARS:
            comp2.params['%s:in:%s' % (name2, var_name)] = comp1.params['%s:in:%s' % (name1, var_name)]

    def _clear_unknowns(self, name, unknowns, var_names=FLOW_VARS):
        '''Reset all of a FlowStation's unknowns to empty (-1.0).'''
        for var_name in var_names:
            unknowns['%s:out:%s' % (name, var_name)] = FLOW_VARS[var_name][1]

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
            return unknowns[output_name] if unknowns[output_name] != FLOW_VARS[var_name][1] else params[param_name]
        def set_vars(var_dict):
            for var_name, val in var_dict.iteritems():
                output_name = '%s:out:%s' % (name, var_name)
                unknowns[output_name] = val
        try:
            out = flowstation.solve(ht=var('ht'), Tt=var('Tt'), Pt=var('Pt'), s=var('s'), hs=var('hs'), Ts=var('Ts'), Ps=var('Ps'), Mach=var('Mach'), area=var('area'), W=var('W'), is_super=var('is_super'))
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
        except flowstation.ArgumentError:
            for var_name in FLOW_VARS.keys():
                set_vars({var_name: var(var_name)})
    
    def _add_flowstation(self, name):
        '''Add a variable tree representing a FlowStation. Parameters are stored as self.parameters['FLOWSTATION NAME:in:VARIABLE NAME'] and outputs are stored as self.unknowns['FLOWSTATION NAME:out:VARIABLE NAME'].'''
        for var_name, props in FLOW_VARS.iteritems():
            param_name = '%s:in:%s' % (name, var_name)
            output_name = '%s:out:%s' % (name, var_name)
            if props[2] == '':
                self.add_param(param_name, props[1], desc=props[0])
                self.add_output(output_name, props[1], desc=props[0])
            else:
                self.add_param(param_name, props[1], desc=props[0], units=props[2])
                self.add_output(output_name, props[1], desc=props[0], units=props[2])
