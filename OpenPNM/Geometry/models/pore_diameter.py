r"""
===============================================================================
Submodule -- pore_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst


def sphere(pore_seed, psd_name='weibull_max',psd_shape=1,psd_loc=1,psd_scale=1,pdf=None):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    if pdf == None:  # Fetch pdf from scipy.stats
        prob_fn = getattr(spst,psd_name)
        pdf = prob_fn(psd_shape,loc=psd_loc,scale=psd_scale)
    pore_diameter = pdf.ppf(pore_seed)
    return pore_diameter
  