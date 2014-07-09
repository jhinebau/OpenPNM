from OpenPNM.Geometry import models as gm

def BoundaryPores(network,pores):
    r"""
    Custom function for generating geometry for Toray090 GDL materials
    """
   
    network.add_method(propname='pore.seed',
                       pores=pores,
                       static=True,
                       model=gm.pore_misc.constant,
                       value=1)
                  
    network.add_method(propname='pore.diameter',
                       pores=pores,
                       static=True,
                       model=gm.pore_diameter.constant,
                       value=0)
                       
    network.add_method(propname='pore.volume',
                       pores=pores,
                       static=True,
                       model=gm.pore_misc.constant,
                       value=0)
                       