**This code is under development and not yet suitable for general use.**

### Guide to pyCycle for OpenMDAO 1.0

**Files updated for OpenMDAO 1.0**

- flowstation.py
- cyclecomponent.py
- compressor.py
- test/test_compressor.py
- test/test_util.py (new in this fork)

**Files not ready for OpenMDAO 1.0**

- duct.py
- heat_exchanger.py
- inlet.py
- nozzle.py
- splitter.py
- start.py
- test/test_duct.py
- test/test_flowstation.py
- test/test_heat_exchanger.py
- test/test_inlet.py
- test/test_nozzle.py
- test/test_splitter.py
- test/test_start.py

**Changes**

OpenMDAO no longer treats variable trees as distinct classes, so flow_station.py has been restructured as a library of stateless functions that calculate thermodynamics values directly from parameters. Passing known parameters to `solve()` will be sufficient for nearly all applications. FlowStation variable trees must be stored manually by components.

The `CycleComponent` class has been given new helper functions for initializing, reseting, and solving FlowStation variable trees. The easiest way to use FlowStations is to subclass `CycleComponent` and use `self._add_flowstation()` and `self._solve_flow_vars()`.

### pyCycle

This plugin provides basic thermodynamic cycle modeling tools for OpenMDAO. It's main feature is the inclusion of a FlowStation variable that provides access to isentropic flow relationships.

# Pre-Reqs


### OpenMDAO 
This is an OpenMDAO plugin, so we assume you have already installed a version of OpenMDAO. 

#### Mac OS-X Tip
On OS-X it is strongly suggested that you setup Python, and other OpenMDAO pre-reqs with 
homebrew. You can follow these [detailed instructions](http://www.lowindata.com/2013/installing-scientific-python-on-mac-os-x/)
but once you have homebrew installed and setup, here is the short version: 


    brew install git
    brew install python
    brew install gfortran
    pip install numpy
    pip install scipy
    brew install freetype
    pip install matplotlib


### Cantera
In addition, this plugin requires [Cantera](https://code.google.com/p/cantera/) 
and the python wrapper for it. You can [compile cantera from scratch](http://cantera.github.io/docs/sphinx/html/compiling.html), 
or follow the instructions below for a bit easier route. 



#### Windows
Cantera provides [pre-compiled](https://code.google.com/p/cantera/downloads/list) binaries 
that are by far the easiest choice for windows. You should follow their instructions 
for the best way to install it and make sure it works. 

But note that you need to install version 2.0.2 and **NOT** the 2.1 beta version. Get only the official 
release version of it! 

https://code.google.com/p/cantera/wiki/WindowsInstallation

#### Mac OS-X
Assuming you've used homebrew to get OpenMDAO setup, then just use it to install Cantera too! 


    brew install cantera


#### Linux
You're best bet is to follow the [official instructions](http://cantera.github.io/docs/sphinx/html/compiling.html). 


# Installation
Make sure you're in an activated OpenMDAO environment. Then you have two options: 

## Regular install
This will install the plugin, but you won't have access to the source in order to mess around with it. 

    plugin install --github pycycle


## Source install
If you want to have access to the source and make changes to it while you're working then you should do 
an install from source. Don't be scared by the cantera source install though, compared to that this is pretty 
easy! It's only two steps. 

1)Clone the repository to your local machine. If you want to clone ours: 
    
    git clone https://github.com/OpenMDAO-Plugins/pyCycle.git
    

Or you can fork it and clone that to your local machine. Either way... 

2) cd into the pyCycle directory
    
    python setup.py develop
    

