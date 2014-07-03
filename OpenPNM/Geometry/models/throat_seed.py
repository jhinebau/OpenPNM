r"""
===============================================================================
Submodule -- throat_seeds
===============================================================================

"""
import scipy as sp


def neighbor_min(throats,network):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    Pn = network.find_connected_pores(throats)
    pore_seeds = network['pore.seed'][Pn]
    throat_seed = sp.amin(pore_seeds,axis=1)
    return throat_seed





