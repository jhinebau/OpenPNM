r"""
===============================================================================
Submodule -- pore_seeds
===============================================================================

"""
import scipy as _sp

def random(geometry,seed=None):
    r"""
    Assign random number to pore bodies for later use in pore size distributions
    """
    pores = geometry.pores()
    _sp.random.seed(seed)  # Set seed for generator
    pore_seeds = _sp.random.rand(_sp.shape(pores)[0])
    return pore_seeds