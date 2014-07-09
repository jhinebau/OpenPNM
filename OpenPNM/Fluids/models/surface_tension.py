r"""
===============================================================================
Submodule -- surface_tension
===============================================================================

Some text here?

"""

import scipy as sp

def Eotvos(fluid,**kwargs):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    k = 2.1e-7  # The Eotvos constant
    Tc = fluid['pore.Tc']
    T = fluid['pore.temperature']
    Vm = 1/fluid['pore.molar_density']
    value = k*(Tc-6-T)/(Vm**(2/3))
    return value

def GuggenheimKatayama(fluid,network,propname,K2,n,**kwargs):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    T = fluid.get_pore_data(prop='temperature')
    Pc = fluid.get_pore_data(prop='Pc')
    Tc = fluid.get_pore_data(prop='Tc')
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    value = sigma_o*(1-T/Tc)**n
    fluid.set_pore_data(prop=propname,data=value)

def BrockBird_scaling(fluid,network,propname,sigma_o,To,**kwargs):
    r"""
    Uses Brock_Bird model to adjust surface tension from it's value at a given reference temperature to temperature of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    sigma_o : float
        Surface tension at reference temperature (N/m)

    To : float
        Temperature at reference conditions (K)
    """
    Tc = fluid.get_pore_data(prop='Tc')
    Ti = fluid.get_pore_data(prop='temperature')
    Tro = To/Tc
    Tri = Ti/Tc
    value = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    fluid.set_pore_data(prop=propname,data=value)



