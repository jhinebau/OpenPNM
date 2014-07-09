r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as sp
import scipy.stats as spst



def cylinder(throat_diameter,**kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = throat_diameter
    throat_area = sp.constants.pi/4*(D)**2
    return throat_area



