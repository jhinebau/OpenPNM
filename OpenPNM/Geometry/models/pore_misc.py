r"""
===============================================================================
Submodule -- pore_misc
===============================================================================
Miscillaneous and generic methods for pores
"""
import scipy as _sp


def constant(network,pores,value=0):
    r"""
    Assign specified number to pores
    """
    values = _sp.ones(_sp.shape(pores)[0])*value
    return values

def random(network,pores,seed=None):
    r"""
    Assign random number to pores
    """
    _sp.random.seed(seed)  # Set seed for generator
    values = _sp.random.rand(_sp.shape(pores)[0])
    return values
    