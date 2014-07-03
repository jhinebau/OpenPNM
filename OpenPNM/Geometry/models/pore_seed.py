r"""
===============================================================================
Submodule -- pore_seeds
===============================================================================

"""
import scipy as sp

def random(pores,seed=None):
    r"""
    Assign random number to pore bodies for later use in pore size distributions
    """
    if type(pores) != int:
        pores = sp.shape(pores)[0]
    sp.random.seed(seed)  # Set seed for generator
    pore_seeds = sp.random.rand(pores)
    return pore_seeds