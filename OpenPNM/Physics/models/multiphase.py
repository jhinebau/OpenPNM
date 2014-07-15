r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as sp

def conduit_conductance(network,
                        fluid,
                        throats,
                        throat_conductance,
                        throat_occupancy='throat.occupancy',
                        pore_occupancy='pore.occupancy',
                        mode='strict',
                        factor=1e-5,
                        **kwargs):
    r"""
    Add a new multiphase conductance property to the conduits of network, where a 
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas. 
    
    This method "closes" conduits that are not sufficiently filled with the 
    specified fluid by multiplying the original conductance by a very small *factor*.

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
        The fluid of interest
        
    occupied_condition : 'occupancy'
        The name of the pore and throat property that dictates whether conduit is "closed" or not
        
    mode : 'strict' or 'medium' or 'loose'
        How agressive the method should be in "closing" conduits. 
        'strict' implies that if any pore or throat in the conduit is unoccupied by the given fluid, the conduit is closed.
        'medium' implies that if either the throat or both pores are unoccupied, the conduit is closed
        'loose' will only close the conduit if the throat is unoccupied.
        
    factor : 1/1e3
        The "closing factor" which becomes multiplied to the original conduit's conductance to severely limit transport.

    Notes
    -----
    This function requires that all the necessary fluid properties already be 
    calculated.

    """
    if (mode == 'loose'):
        closed_conduits = -fluid[throat_occupancy]
    else:
        throats_closed = -fluid[throat_occupancy]
        connected_pores = network.find_connected_pores(throats)
        pores_1 = connected_pores[:,0]
        pores_2 = connected_pores[:,1]
        pores_1_closed = -fluid[pore_occupancy][pores_1]
        pores_2_closed = -fluid[pore_occupancy][pores_2]
        if(mode == 'medium'):
            closed_conduits = throats_closed | (pores_1_closed & pores_2_closed)
            
        if(mode == 'strict'):
            closed_conduits = pores_1_closed | throats_closed | pores_2_closed
    open_conduits = -closed_conduits
    throat_value = fluid[throat_conductance]
    value = throat_value*open_conduits + throat_value*closed_conduits*factor
    return value

def late_throat_filling(network,fluid,throats,Pc_star,eta):
    pass

def late_pore_filling(network,
                      fluid,
                      pores,
                      Pc,
                      Swp_star=0.2,
                      eta=3,
                      wetting_fluid=False,
                      pore_occupancy='pore.occupancy',
                      throat_capillary_pressure='throat.capillary_pressure',
                      **kwargs):
    r'''
    Applies a late pore filling model to calculate fractional pore filling as 
    a function of applied capillary pressure.
    
    Parameters
    ----------
    Pc : float
        The capillary pressure in the non-wetting phase (Pc > 0)
    Swp_star : float
        The residual wetting phase in an invaded pore immediately after
        nonwetting phase invasion
    eta : float
        Exponent to control the rate at which wetting phase is displaced
    wetting_fluid : boolean
        Indicates whether supplied fluid is the wetting or non-wetting phase
    
    
    '''
    
    prop = fluid[throat_capillary_pressure]
    neighborTs = network.find_neighbor_throats(pores,flatten=False)
    Pc_star = sp.array([sp.amin(prop[row]) for row in neighborTs])
    Swp = Swp_star*(Pc_star/Pc)**eta
    if wetting_fluid:
        values = Swp*fluid[pore_occupancy]*(Pc_star<Pc)
    else:
        values = (1-Swp)*(~fluid[pore_occupancy])*(Pc_star<Pc)
    return values
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    