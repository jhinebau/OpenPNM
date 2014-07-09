from OpenPNM.Geometry import models as gm

def Toray090(network,pores,throats):
    r"""
    Custom function for generating geometry for Toray090 GDL materials
    """
    network.add_method(propname='pore.seed',
                       pores=pores,
                       static=False,
                       model=gm.pore_seed.random,
                       seed=None)
    network.add_method(propname='throat.seed',
                       throats=throats,
                       static=True,
                       model=gm.throat_seed.neighbor,
                       mode='min')
    network.add_method(propname='pore.diameter',
                       pores=pores,
                       static=False,
                       model=gm.pore_diameter.sphere,
                       psd_name='weibull_min',
                       psd_shape=2.5,
                       psd_loc=5e-4,
                       psd_scale=4e-6)
    network.add_method(propname='throat.diameter',
                       throats=throats,
                       static=True,
                       model=gm.throat_diameter.cylinder,
                       tsd_name='weibull_min',
                       tsd_shape=2.5,
                       tsd_loc=5e-4,
                       tsd_scale=4e-6)
    network.add_method(propname='pore.volume',
                       pores=pores,
                       static=True,
                       model=gm.pore_volume.sphere)
    network.add_method(propname='throat.length',
                       throats=throats,
                       static=True,
                       model=gm.throat_length.straight)
    network.add_method(propname='throat.volume',
                       throats=throats,
                       static=True,
                       model=gm.throat_volume.cylinder)
                       
                       

