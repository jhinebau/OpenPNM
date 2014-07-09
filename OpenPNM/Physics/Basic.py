from OpenPNM.Physics import models as pm

def Phys(fluid,network,pores,throats):
    r"""
    Custom function for generating geometry for Toray090 GDL materials
    """
    fluid.add_method(propname='throat.capillary_pressure',
                     throats=throats,
                     static=True,
                     model=pm.capillary_pressure.washburn,
                     seed=None)



