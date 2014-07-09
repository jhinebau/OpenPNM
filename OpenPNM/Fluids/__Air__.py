import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
from OpenPNM.Fluids.__GenericFluid__ import GenericFluid

class Air(GenericFluid):
    r"""
    Creates Fluid object with a default name 'air' and preset values for air
    
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
    >>> air = OpenPNM.Fluids.Air(network=pn)
    """
    def __init__(self,name=None,**kwargs):
        super(Air,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        # Set default T and P since most propery models require it
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0
        self['pore.Tc'] = 132.65
        self['pore.Pc'] = 3.771e6
        self['pore.MW'] = 0.0291
        self.add_method(propname='pore.diffusivity',
                        model=OpenPNM.Fluids.models.diffusivity.Fuller,
                        static=False,
                        MA=0.03199,
                        MB=0.0291,
                        vA=16.3,
                        vB=19.7)
        self.add_method(propname='pore.viscosity',
                        model=OpenPNM.Fluids.models.misc.constant,
                        static=True,
                        value=1.75e-5)
        self.add_method(propname='pore.molar_density',
                        model=OpenPNM.Fluids.models.molar_density.ideal_gas,
                        static=False,
                        R=8.314)


if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Fluids.Air(network=pn)