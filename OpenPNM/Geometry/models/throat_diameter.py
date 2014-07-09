r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
import scipy as sp
import scipy.stats as spst



def cylinder(network,throats,tsd_name='weibull_max',tsd_shape=1,tsd_loc=1,tsd_scale=1,pdf=None):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    if pdf == None:  # Fetch pdf from scipy.stats
        prob_fn = getattr(spst,tsd_name)
        pdf = prob_fn(tsd_shape,loc=tsd_loc,scale=tsd_scale)
    throat_diameter = pdf.ppf(network['throat.seed'])
    return throat_diameter

