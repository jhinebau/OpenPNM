r"""
===============================================================================
Submodule -- contact_angle
===============================================================================

"""
import scipy as _sp

def constant(fluid,value=0,**kwargs):
    r"""
    Assign specified number to pores
    """
    values = _sp.ones((fluid.Np,))*value
    return values

def random(fluid,seed=None,**kwargs):
    r"""
    Assign random number to pores
    """
    _sp.random.seed(seed)  # Set seed for generator
    values = _sp.random.rand(fluid.Np,)
    return values