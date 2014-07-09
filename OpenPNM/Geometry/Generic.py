from OpenPNM.Geometry import models as gm

def StickAndBall(network,pores,throats):
    r"""
    Custom function for generating a generic StickAndBall geometry for testing
    purposes.
    """
   
    network.add_method(propname='pore.seed',
                       pores=pores,
                       static=True,
                       model=gm.pore_seed.random,
                       seed=1)
                  
    network.add_method(propname='throat.seed',
                       throats=throats,
                       static=True,
                       model=gm.throat_seed.neighbor,
                       mode='min')
                       
    network.add_method(propname='pore.diameter',
                       pores=pores,
                       static=True,
                       model=gm.pore_diameter.sphere,
                       psd_name='weibull_min',
                       psd_shape=1,
                       psd_loc=1,
                       psd_scale=1)
                       
    network.add_method(propname='throat.diameter',
                       throats=throats,
                       static=True,
                       model=gm.throat_diameter.cylinder,
                       tsd_name='weibull_min',
                       tsd_shape=1,
                       tsd_loc=1,
                       tsd_scale=1)
                       
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
                       
                       
                       

