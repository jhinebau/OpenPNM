r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as sp

def sphere(network,pores):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    pore_volumes = sp.pi/6*network['pore.diameter']**3
    return pore_volumes
    
