import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
from OpenPNM.Fluids.__GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r'''
    Creates Fluid object with a default name 'water' and preset values for water
    
    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this fluid object will be attached.  
        
    Notes
    -----
    This explicit association is necessary so the Fluid object can initialize
    data arrays of the correct size to store network data.
    
    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> water = OpenPNM.Fluids.Water(network=pn)
    '''
    def __init__(self,name=None,**kwargs):
        super(Water,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        # Set default T and P since most propery models require it
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0
        self['pore.Tc'] = 647.096
        self['pore.Pc'] = 22.06e6
        self['pore.MW'] = 0.0180
        self.add_method(model=OpenPNM.Fluids.models.misc.constant,
                        propname='pore.diffusivity',
                        static=True,
                        value=2e-9)
        self.add_method(model=OpenPNM.Fluids.models.misc.constant,
                        propname='pore.viscosity',
                        static=True,
                        value=0.001)
        self.add_method(model=OpenPNM.Fluids.models.misc.constant,
                        propname='pore.molar_density',
                        static=True,
                        R=44445)
        self.add_method(model=OpenPNM.Fluids.models.misc.constant,
                        propname='pore.surface_tension',
                        static=True,
                        R=0.072)
        self.add_method(model=OpenPNM.Fluids.models.misc.constant,
                        propname='pore.contact_angle',
                        static=True,
                        R=110)

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Fluids.Water(network=pn)
