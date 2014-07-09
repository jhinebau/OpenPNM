r"""
===============================================================================
Submodule -- throat_seeds
===============================================================================

"""
import scipy as _sp

def constant(network,throats,value=0,**kwargs):
    r"""
    Assign specified number to throat for later use in throat size distributions
    """
    pore_seeds = _sp.ones(_sp.shape(pores)[0])*value
    return pore_seeds

def neighbor(network,throats,mode='min',**kwargs):
    r"""
    Adopt the seed value from neighboring pores
    
    Parameters
    ----------
    mode : String
        Controls how to determine the throat seed value.  Options are:
        
        - min : (default) Takes the minimum value of it's two neighbors
        - max : Takes the maximum value of it's two neighbors
        - mean : Takes the arithmetic mean value of it's two neighbors
        
    """
    Pn = network.find_connected_pores(throats)
    pore_seeds = network['pore.seed'][Pn]
    if mode=='min':
        throat_seed = _sp.amin(pore_seeds,axis=1)
    elif mode=='max':
        throat_seed = _sp.amax(pore_seeds,axis=1)
    elif mode=='mean':
        throat_seed = _sp.mean(pore_seeds,axis=1)
    return throat_seed





