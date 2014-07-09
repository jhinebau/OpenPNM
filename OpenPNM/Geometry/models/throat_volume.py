r"""
===============================================================================
Submodule -- throat_volume
===============================================================================

"""
import scipy as sp


def cylinder(network,throats,**kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    throat_length = network['throat.length'][throats]
    throat_diameter = network['throat.diameter'][throats]
    throat_volume = sp.pi/4*throat_length*throat_diameter**2
    return throat_volume


