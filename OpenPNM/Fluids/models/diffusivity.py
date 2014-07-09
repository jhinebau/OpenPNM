r"""
===============================================================================
Submodule -- diffusivity
===============================================================================

"""
import scipy as sp

def constant(fluid,network,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    fluid.set_pore_data(prop=propname,data=value)

def na(fluid,network,propname,**params):
    r"""
    Assigns nonsensical, but numerical value of -1.  
    This ensurse stability of other methods 
    but introduces the possibility of being misused.
    """
    value = -1
    fluid.set_pore_data(prop=propname,data=value)

def Fuller(fluid, network, propname, MA, MB, vA, vB, **params): 
    r"""
    Uses Fuller model to estimate diffusion coefficient for gases from first 
    principles at conditions of interest

    Parameters
    ----------
    T :  float, array_like
        Temperature of interest [K]
    P :  float, array_like
        Pressure of interest [Pa]
    MA : float, array_like
        Molecular weight of component A [kg/mol]
    MB : float, array_like
        Molecular weight of component B [kg/mol]
    VA:  float, array_like
        Sum of atomic diffusion volumes for component A
    VB:  float, array_like
        Sum of atomic diffusion volumes for component B
    """

    T = fluid.get_pore_data(prop='temperature')
    P = fluid.get_pore_data(prop='pressure')
    MAB = 2*(1/MA+1/MB)**(-1)
    MAB = MAB*1e3
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3)+vB**(1./3))**2)*1e-4
    fluid.set_pore_data(prop=propname,data=value)

def Fuller_scaling(fluid,network,propname,DABo, To, Po,**params):
    r"""
    Uses Fuller model to adjust a diffusion coefficient for gases from reference conditions to conditions of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    Po, To : float, array_like
        Pressure & temperature at reference conditions, respectively
    """
    Ti = fluid.get_pore_data(prop='temperature')
    Pi = fluid.get_pore_data(prop='pressure')
    value = DABo*(Ti/To)**1.75*(Po/Pi)
    fluid.set_pore_data(prop=propname,data=value)

def TynCalus(fluid,network,propname,VA ,VB ,sigma_A ,sigma_B,viscosity='viscosity',**params):
    r"""
    Uses Tyn_Calus model to estimate diffusion coefficient in a dilute liquid solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    T :  float, array_like
        Temperature of interest (K)
    viscosity :  float, array_like
        Viscosity of solvent (Pa.s)
    VA : float, array_like
        Molar volume of component A at boiling temperature (m3/mol)
    VB : float, array_like
        Molar volume of component B at boiling temperature (m3/mol)
    sigmaA:  float, array_like
        Surface tension of component A at boiling temperature (N/m)
    sigmaB:  float, array_like
        Surface tension of component B at boiling temperature (N/m)

    """
    T = fluid.get_pore_data(prop='temperature')
    mu = fluid.get_pore_data(prop=viscosity)
    value = 8.93e-8*(VB*1e6)**0.267/(VA*1e6)**0.433*T*(sigma_B/sigma_A)**0.15/(mu*1e3)
    fluid.set_pore_data(prop=propname,data=value)

def TynCalus_Scaling(fluid,network,propname,DABo ,To ,mu_o ,viscosity='viscosity', **params):
    r"""
    Uses Tyn_Calus model to adjust a diffusion coeffciient for liquids from reference conditions to conditions of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    mu_o, To : float, array_like
        Viscosity & temperature at reference conditions, respectively
    """
    Ti = fluid.get_pore_data(prop='temperature')
    mu_i = fluid.get_pore_data(prop=viscosity)
    value = DABo*(Ti/To)*(mu_o/mu_i)
    fluid.set_pore_data(prop=propname,data=value)
