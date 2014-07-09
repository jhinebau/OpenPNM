r"""
===============================================================================
Submodule -- thermal_conductivity
===============================================================================

"""
import scipy as sp

def Chung(fluid,Cv,MW,acentric,**kwargs):
    r"""
    Uses Chung et al. model to estimate thermal conductivity for gases with 
    ow pressure(<10 bar) from first principles at conditions of interest

    Parameters
    ----------
    Cv :  float, array_like
        Heat capacity at constant volume (J/(mol.K))
    MW : float, array_like
        Molecular weight of the component (kg/mol)
    acentric : float, array_like
        Acentric factor of the component

    """
    R = 8.314
    T = fluid['pore.temperature']
    mu = fluid['pore.viscosity']
    Tc = fluid['pore.Tc']
    Tr = T/Tc
    z = 2.0 + 10.5*Tr**2
    beta = 0.7862 - 0.7109*acentric + 1.3168*acentric**2
    alpha = Cv/R -3/2
    s = 1 + alpha*((0.215+0.28288*alpha-1.061*beta+0.26665*z)/(0.6366+beta*z+1.061*alpha*beta))
    value = 3.75*s*(mu)*R/(MW)
    return value

def Sato(fluid,network,propname,Tb,MW,**kwargs):
    r"""
    Uses Sato et al. model to estimate thermal conductivity for pure liquids 
    from first principles at conditions of interest

    Parameters
    ----------
    Tb :  float, array_like
        Boiling temperature of the component (K)
    MW : float, array_like
        Molecular weight of the component (kg/mol)

    """
    T = fluid['pore.temperature']
    Tc = fluid['pore.Tc']
    Tbr = Tb/Tc
    Tr = T/Tc
    value = (1.11/((MW*1e3)**0.5))*(3+20*(1-Tr)**(2/3))/(3+20*(1-Tbr)**(2/3))
    return value
