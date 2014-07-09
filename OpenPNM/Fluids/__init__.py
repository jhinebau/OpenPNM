r"""
*******************************************************************************
:mod:`OpenPNM.Fluids` -- Fluid Property Estimation Methods
*******************************************************************************

.. module:: OpenPNM.Fluids

Contents
--------
This module contains methods for estimating fluid properties

.. autoclass:: GenericFluid
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: Water
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: Air
   :members:
   :undoc-members:
   :show-inheritance:

"""

from .__GenericFluid__ import GenericFluid
from .__Water__ import Water
from .__Air__ import Air
from . import models
