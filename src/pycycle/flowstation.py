__all__ = ['CanteraFlowStation']

from os.path import dirname, join

from openmdao.main.api import VariableTree
from openmdao.lib.datatypes.api import Float, VarTree, Enum

from Cantera import *

import pycycle #used to find file paths

GAS_CONSTANT = 0.0685592 #BTU/lbm-R

#secant solver with a limit on overall step size
def secant(func, x0, TOL=1e-5, x_min=1e15, x_max=1e15 ):
    if x0 >= 0:
        x1 = x0*(1 + 1e-4) + 1e-4
    else:
        x1 = x0*(1 + 1e-4) - 1e-4
    f1, f = func(x1), func(x0)
    if (abs(f) > abs(f1)):
        x1, x0 = x0, x1 
        f1, f = f, f1
    dx = f * (x0 - x1) / float(f - f1)  
    count = 0  
    while 1:
        if abs(dx) < TOL * (1 + abs(x0)): 
        #if abs((f1-f)/(f+1e-10)) < TOL: 
            return x0 -dx
        dx = f * (x0 - x1) / float(f - f1)  
        df = abs((f1-f)/(f+1e-10))      
        if x0-dx < x_min: 
            #x1, x0 = x0, x0*(1+.01*abs(dx)/dx)
            x1, x0 = x0, (x_max+x0)/2
        elif x0-dx > x_max: 
            x1, x0 = x0, (x_max+x0)/2
        else:    
            x1, x0 = x0, x0 - dx
        f1, f = f, func(x0) 
        count = count + 1


class CanteraFlowStation(VariableTree):

    reactants=["Air", "H2O"]

    ht=Float(0.0, desc='total enthalpy', unit='Btu/lbm')
    Tt=Float(0.0, desc='total temperature', unit='R')
    Pt=Float(0.0, desc='total pressure', unit='lbf/in**2')
    rhot=Float(0.0, desc='total density', unit='lbm/ft**3') 
    gamt=Float(0.0, desc='total gamma', unit='') 
    s =Float(0.0, desc='entropy', unit='Btu/lbm-R')
    W =Float(0.0, desc='weight flow', unit='lbm/s') 
    FAR =Float(0.0, desc='fuel to air ratio', unit='') 
    WAR =Float(0.0, desc='water to air ratio', unit='') 
    hs=Float(0.0, desc='static enthalpy', unit='Btu/lbm')
    Ts=Float(0.0, desc='static temperature', unit='R')
    Ps=Float(0.0, desc='static pressure', unit='lbf/in**2')
    rhos=Float(0.0, desc='static density', unit='lbm/ft**3')
    gams=Float(0.0, desc='static gamma', unit='')    
    Vflow =Float(0.0, desc='Velocity', unit='ft/s')   
    Mach=Float(0.0, desc='Mach number', unit='')
    area =Float(0.0, desc='flow area', unit='in**2') 
    sub_or_super = Enum(('sub','super'), desc="selects preference for subsonic or supersonice solution when setting area")
    
    Wc = Float(0.0, desc='corrected weight flow', unit='lbm/s') 


    #intialize station        
    def __init__(self,*args,**kwargs): 
        super(CanteraFlowStation, self).__init__(*args,**kwargs)

        #properties file path
        _dir = dirname(pycycle.__file__)
        _prop_file = join(_dir,'gri1000.cti')

        self._trigger = 0
        
        self._species=[1.0, 0, 0, 0, 0, 0, 0, 0]
        self._mach_or_area=0    
        self._flow=importPhase(_prop_file)
        self._flowS=importPhase(_prop_file)
        self.setDryAir()
        self.setTotalTP(518,15)
    
    #add a reactant that can be mixed in
    def add_reactant(self, reactant):
            for i in range(0, len(self.reactants)):
                if self.reactants[i] == reactant:
                    return
            self.reactants.append(reactant)

    def _W_changed(self): 
        if self._trigger == 0:
            self._trigger=1
            self.setStatic()
            self._trigger=0

    #trigger action on Mach
    def _Mach_changed(self):
        if self._trigger == 0:
            self._trigger=1
            self._mach_or_area=1
            self.setStatic()
            self._trigger=0
                    
    #trigger action on area        
    def _area_changed(self):
        if self._trigger == 0:
            self._trigger=1
            self._mach_or_area=2
            self.setStatic()
            self._trigger=0
           
    #trigger action on static pressure       
    def _Ps_changed(self):
        if self._trigger == 0:
            self._trigger=1
            self._mach_or_area=3
            self.setStatic()
            self._trigger=0 


    
    #set the composition to dry air
    def setDryAir(self):
        self._flow.setMassFractions("Air:1.00")             
        self._species=(len(self.reactants)+1)*[0, ]
        self._species[0]=1 
        self.WAR=0
        self.FAR=0
        self.setStatic()
        self._trigger=0
    
    #set the compositon to air with water
    def setWAR(self, WAR):
        self._trigger=1
        self.WAR=WAR
        self.FAR=0
        self._flow.setMassFractions("Air:"+str((1-WAR)/(1+WAR))+" H2O:"+str((WAR)/(1+WAR)))
        self._species=len(self.reactants)*[0, ]
        self._species[0]=(1-WAR)/(1+WAR)
        self._species[1]=(WAR)/(1+WAR)
        self.setStatic()
        self._trigger=0

        #set total conditions based on T an P
    def setTotalTP(self, Tin, Pin):
        self._trigger=1
        self.Tt=Tin
        self.Pt=Pin                
        self._flow.set(T=Tin*5./9., P=Pin*6894.75729)
        self._flow.equilibrate('TP')
        self.ht=self._flow.enthalpy_mass()*0.0004302099943161011
        self.s=self._flow.entropy_mass()*0.000238845896627
        self.rhot=self._flow.density()*.0624
        self.Tt=self._flow.temperature()*9./5.
        self.gamt=self._flow.cp_mass()/self._flow.cv_mass()
        self._flowS=self._flow 
        self.setStatic()
        self._calculated_properties()
        self._trigger=0

    #set total conditions based on h and P
    def setTotal_hP(self, hin, Pin):
        self._trigger=1 
        self.ht=hin
        self.Pt=Pin
        self._flow.set(H=hin/.0004302099943161011, P=Pin*6894.75729)
        self._flow.equilibrate('HP')
        self.Tt=self._flow.temperature()*9./5.
        self.s=self._flow.entropy_mass()*0.000238845896627    
        self.rhot=self._flow.density()*.0624
        self.gamt=self._flow.cp_mass()/self._flow.cv_mass()
        self.setStatic()
        self._calculated_properties()
        self._trigger=0

    #set total condition based on S and P
    def setTotalSP(self, sin, Pin):
        self._trigger=1
        self.s=sin
        self.Pt=Pin             
        self._flow.set(S=sin/0.000238845896627, P=Pin*6894.75729)
        self._flow.equilibrate('SP', loglevel=1)
        self.Tt=self._flow.temperature()*9./5.
        self.ht=self._flow.enthalpy_mass()*0.0004302099943161011
        self.rhot=self._flow.density()*.0624
        self.gamt=self._flow.cp_mass()/self._flow.cv_mass()
        self.setStatic()
        self._calculated_properties()
        self._trigger=0

    #add another station to this one
    #mix enthalpies and keep pressure and this stations value
    def add(self, FS2):
        temp =""
        for i in range(0, len(self._species)):
                self._species[i]=(self.W*self._species[i]+FS2.W*FS2._species[i])/(self.W + FS2.W)
                temp=temp+self.reactants[i]+":"+str(self._species[i])+" "
        self._flow.setMassFractions(temp)      
        air1 = self.W * ( 1. / ( 1. + self.FAR + self.WAR ));
        air2 = FS2.W *( 1. / ( 1 + FS2.WAR + FS2.FAR ));
        self.FAR = ( air1 * self.FAR + air2*FS2.FAR )/( air1 + air2 );
        self.WAR = ( air1 * self.WAR + air2*FS2.WAR )/( air1 + air2 );
        self.ht=(self.W*self.ht+FS2.W+FS2.ht)/(self.W+FS2.W)
        self._flow.setMassFractions(temp)
        self.W=self.W +(FS2.W)
        self._flow.set(H=self.ht/0.0004302099943161011, P=self.Pt*6894.75729)
        self._flow.equilibrate('HP')
        self.Tt=self._flow.temperature()*9./5.
        self.s=self._flow.entropy_mass()* 0.000238845896627
        self.rhot=self._flow.density()*.0624
        self.gamt=self._flow.cp_mass()/self._flow.cv_mass()          
                    
    def copy_from(self, FS2):
        """duplicates total properties from another flow station""" 
        self.ht=FS2.ht
        self.Tt=FS2.Tt
        self.Pt=FS2.Pt
        self.rhot=FS2.rhot
        self.gamt=FS2.gamt
        self.s =FS2.s
        self.W =FS2.W
        self.FAR =FS2.FAR
        self.WAR =FS2.WAR
        temp =""
        for i in range(0, len(self.reactants)):
                self._species[i]=FS2._species[i]
                temp=temp+self.reactants[i]+":"+str(self._species[i])+" "
        self._flow.setMassFractions(temp)
        self._flow.set(T=self.Tt*5./9., P=self.Pt*6894.75729)
        self._flow.equilibrate('TP')
                    
    #burn a fuel with this station        
    def burn(self, fuel, Wfuel, hfuel):
        flow_1=self.W
        ht=self.ht
        self.W=self.W + Wfuel 
        for i in range(0, len(self.reactants)):
            if fuel == self.reactants[i]:
                self._species[i]=(flow_1*self._species[i]+Wfuel)/ self.W
            else:
                self._species[i]=(flow_1*self._species[i])/ self.W
        ht=(flow_1 * ht + Wfuel * hfuel)/ self.W
        self.ht=ht
        air1=flow_1 * (1. / (1. + self.FAR + self.WAR))
        self.FAR=(air1 * self.FAR + Wfuel)/(air1)
        temp =""
        for i in range(0, len(self.reactants)):
            temp=temp+self.reactants[i]+":"+str(self._species[i])+" "
        self._flow.setMassFractions(temp)        
        self._flow.set(H=ht/0.0004302099943161011, P=self.Pt*6894.75729)
        self._flow.equilibrate('HP')
        self.Tt=self._flow.temperature()*9./5.
        self.s=self._flow.entropy_mass()*0.000238845896627  
        self.rhot=self._flow.density()*.0624
        self.gamt=self._flow.cp_mass()/self._flow.cv_mass() 
       
    #set the statics based on Mach
    def setStaticMach(self):
        mach_target = self.Mach
        def f(Ps):
            self.Ps=Ps
            self.setStaticPs()
            return self.Mach - mach_target

        Ps_guess = self.Pt*(1 + (self.gamt-1)/2*mach_target**2)**(self.gamt/(1-self.gamt))
        secant(f, Ps_guess, x_min=0, x_max=self.Pt)


    #set the statics based on pressure
    def setStaticPs(self):
        self._flowS=self._flow 
        self._flowS.set(S=self.s/0.000238845896627, P=self.Ps*6894.75729) 
        self._flowS.equilibrate('SP')
        self.Ts=self._flowS.temperature()*9./5.
        self.rhos=self._flowS.density()*.0624
        self.gams=self._flowS.cp_mass()/self._flowS.cv_mass() 
        self.hs=self._flowS.enthalpy_mass()*0.0004302099943161011                   
        Vson=math.sqrt(self.gams*GasConstant*self._flowS.temperature()/self._flowS.meanMolecularWeight())*3.28084
        self.Vflow=(778.169*32.1740*2*(self.ht-self.hs))**.5
        self.Mach=self.Vflow / Vson
        self.area= self.W / (self.rhos*self.Vflow)*144. 

    def setStaticArea(self): 
        target_area = self.area
        Ps_guess=self.Pt*(1 + (self.gamt-1)/2)**(self.gamt/(1-self.gamt)) #at mach 1
        def f(Ps):
            self.Ps = Ps
            self.setStaticPs()
            return 1-self.Mach
        
        Ps_M1 = secant(f,Ps_guess,x_min=0,x_max=self.Pt)

        #find the subsonic solution first
        guess = (self.Pt+Ps_M1)/2
        def f(Ps):
            self.Ps = Ps
            self.setStaticPs()
            #print "TEST", self.Ps, self.Pt
            return self.W/(self.rhos*self.Vflow)*144.-target_area
        secant(f,  guess, x_min=0, x_max=self.Pt)

        #if you want the supersonic one, just keep going with a little lower initial guess    
        if self.sub_or_super == "super":
            #jsg: wild guess of 1/M_subsonic
            mach_guess = 1/self.Mach
            Ps_guess=self.Pt*(1 + (self.gamt-1)/2*mach_guess**2)**(self.gamt/(1-self.gamt))
            secant(f, Ps_guess, x_min=0, x_max=Ps_M1)

    #determine which static calc to use
    def setStatic(self):
        if (self.Tt and self.Pt): # if non zero
            self.Wc = self.W*(self.Tt/518.67)**.5/(self.Pt/14.696)

        if self._mach_or_area == 0:
            self.Ps = self.Pt 
            self.Ts = self.Tt
            self.rhos = self.rhot
            self.gams = self.gamt
            self.hs = self.ht 
            self.Vflow = 0
            self.Mach = 0

        elif self._mach_or_area == 1:
            self.setStaticMach()

        elif self._mach_or_area ==2:
            self.setStaticArea()
            
        elif self._mach_or_area == 3:
            self.setStaticPs()

    def _calculated_properties(self): 
        self.Wc = self.W*(self.Tt/518.67)**.5/(self.Pt/14.696)


    #set the statics based on Ts, Ps, and MN
    #UPDGRAEDE TO USE LOOPS
    def setStaticTsPsMN(self, Ts, Ps, MN): 
        self._trigger=1 
        self.Tt=Ts*(1+(self.gamt - 1) /2.* MN**2)
        self.Pt=Ps*(1+(self.gamt - 1) /2.* MN**2)**(self.gamt /(self.gamt -1))
        self.setTotalTP(self.Tt, self.Pt)
        self._trigger=1
        self.Mach=MN 
        self.setStaticMach()
        self.area= self.W / (self.rhos * self.Vflow)*144. 
        self._trigger=0

#variable class used in components
class FlowStation(VarTree): 
    def __init__(self,*args,**metadata): 
        super(FlowStation,self).__init__(CanteraFlowStation(), *args, **metadata)
