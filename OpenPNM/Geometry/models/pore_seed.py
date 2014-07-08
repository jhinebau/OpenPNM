r"""
===============================================================================
Submodule -- pore_seeds
===============================================================================

"""
import scipy as _sp

def constant(network,pores,value=0):
    r"""
    Assign specified number to pore bodies for later use in pore size distributions
    """
    pore_seeds = _sp.ones(_sp.shape(pores)[0])*value
    return pore_seeds

def random(network,pores,seed=None):
    r"""
    Assign random number to pore bodies for later use in pore size distributions
    """
    _sp.random.seed(seed)  # Set seed for generator
    pore_seeds = _sp.random.rand(_sp.shape(pores)[0])
    return pore_seeds