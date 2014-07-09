from OpenPNM.Physics import models as pm

def Phys(fluid,network,pores,throats):
    r"""
    A suite of basic physics models typically used on pore network modeling
    """
    fluid.add_method(propname='throat.capillary_pressure',
                     throats=throats,
                     static=True,
                     model=pm.capillary_pressure.washburn)



