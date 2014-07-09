r"""
===============================================================================
Submodule -- throat_length
===============================================================================
The functions in this module provide ways to calculate throat lengths

"""
import scipy as sp


def straight(network,throats):
    r"""
    Calculate throat length 
    """
    #Initialize throat_property['length']
    pores = network['throat.conns'][throats,:]
    C1 = network['pore.coords'][pores[:,0],:]
    C2 = network['pore.coords'][pores[:,1],:]
    E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
    D1 = network['pore.diameter'][pores[:,0]]
    D2 = network['pore.diameter'][pores[:,1]]
    throat_length = E-(D1+D2)/2
    return throat_length
    
