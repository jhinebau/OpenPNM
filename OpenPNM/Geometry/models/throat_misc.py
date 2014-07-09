r"""
===============================================================================
Submodule -- pore_misc
===============================================================================
Miscillaneous and generic methods for pores
"""
import scipy as _sp


def constant(network,throats,value=0,**kwargs):
    r"""
    Assign specified number to throats
    """
    pore_seeds = _sp.ones(_sp.shape(throats)[0])*value
    return pore_seeds

def random(network,throats,seed=None,**kwargs):
    r"""
    Assign random number to throats
    """
    _sp.random.seed(seed)  # Set seed for generator
    pore_seeds = _sp.random.rand(_sp.shape(throats)[0])
    return pore_seeds
    