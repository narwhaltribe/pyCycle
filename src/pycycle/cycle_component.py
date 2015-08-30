import math 

from openmdao.core.component import Component
from pycycle import flowstation

class CycleComponent(Component): 

    def __init__(self): 
        super(CycleComponent, self).__init__()
        self.add_param('design', False, desc='flag to indicate that the calculations are design conditions')

    def _clear_unknowns(self, name, unknowns):
        '''Reset all of a FlowStation's unknowns to empty (-1.0).'''
        for var_name in ('ht', 'Tt', 'Pt', 's', 'hs', 'Ts', 'Ps', 'Mach', 'area', 'W', 'FAR', 'WAR', 'Vsonic', 'Vflow', 'rhos', 'rhot', 'gams', 'gamt', 'Cp', 'Cv', 'Wc'):
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
                  'Wc': out.Wc})
    
    def _add_flowstation(self, name):
        '''Add a variable tree representing a FlowStation. Parameters are stored as self.parameters['FLOWSTATION NAME:in:VARIABLE NAME'] and outputs are stored as self.unknowns['FLOWSTATION NAME:out:VARIABLE NAME'].'''
        def add_output(var_name, desc, units=None):
            self.add_output('%s:out:%s' % (name, var_name), -1.0, desc=desc, units=units)
        def add_duo(var_name, desc, units=None):
            self.add_param('%s:in:%s' % (name, var_name), -1.0, desc=desc, units=units)
            add_output(var_name, desc, units)
        self.add_param('%s:in:is_super' % name, False, desc='selects preference for supersonic versus subsonic solution when setting area')
        add_duo('ht', 'total enthalpy', 'Btu/lbm')
        add_duo('Tt', 'total temperature', 'degR')
        add_duo('Pt', 'total pressure', 'lbf/inch**2')
        add_duo('s', 'entropy', 'Btu/(lbm*R)')
        add_duo('hs', 'static enthalpy', 'Btu/lbm')
        add_duo('Ts', 'static temperature', 'degR')
        add_duo('Ps', 'static pressure', 'lbf/inch**2')
        add_duo('Mach', 'Mach number')
        add_duo('area', 'flow area', 'inch**2')
        add_duo('W', 'weight flow', 'lbm/s')
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
