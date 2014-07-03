r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as sp

def sphere(pore_diameter):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    pore_volumes = sp.pi/6*pore_diameter**3
    return pore_volumes
    
