r"""
===============================================================================
Submodule -- viscosity
===============================================================================

Models for predicting fluid viscosity

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

def Reynolds(fluid,network,propname,uo,b,**params):
    r"""
    Uses Reynolds method for the temperature dependance of shear viscosity

    Parameters
    ----------
    u0, b :  float, array_like
            Coefficients of Reynolds method

    """
    T = fluid.get_pore_data(prop='temperature')
    value = uo*sp.exp(-1*b*T)
    fluid.set_pore_data(prop=propname,data=value)

def Chung(fluid,network,propname,Vc,MW,acentric,kappa,dipole,**params):
    r"""
    Uses Chung et al. model to estimate viscosity for gases with low pressure(not near the critical pressure) from first principles at conditions of interest

    Parameters
    ----------
    Vc :  float, array_like
        Critical volume of the component (m3/mol)
    MW : float, array_like
        Molecular weight of the component (kg/mol)
    acentric : float, array_like
        Acentric factor of the component
    kappa : float, array_like
        Special correction for highly polar substances
    dipole :float, array_like
        Dipole moment (C.m)

    """
    T = fluid.get_pore_data(prop='temperature')
    Tc = fluid.get_pore_data(prop = 'Tc')
    Tr= T/Tc
    Tstar = 1.2593*Tr
    A = 1.161415
    B = 0.14874
    C = 0.52487
    D = 0.77320
    E = 2.16178
    F = 2.43787
    sigma = (A*(Tstar)**(-B)) + C*(sp.exp(-D*Tstar)) + E*(sp.exp(-F*Tstar))
    dipole_r = 131.3*(dipole*2.997e29)/((Vc*1e6)*Tc)**0.5
    f = 1-0.2756*acentric + 0.059035*(dipole_r**4) + kappa
    value = 40.785*f*(MW*1e3*T)**0.5/((Vc*1e6)**(2/3)*sigma)*1e-7
    fluid.set_pore_data(prop=propname,data=value)

