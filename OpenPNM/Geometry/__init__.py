r"""
*******************************************************************************
:mod:`OpenPNM.Geometry` -- Classes related to the creation of pore and throat geometry
*******************************************************************************

.. module:: OpenPNM.Geometry

Contents
--------
Contains methods for applying pore and throat geometry

  
"""
#Import every file in the directory, including both classes and methods
#import os,sys
#dir = os.path.dirname(os.path.abspath(__file__))
#for item in os.listdir(dir):
#    if item.split('.')[-1] == 'py':
#        if item == '__init__.py':
#            pass
#        elif item[0:2] == '__':
#            exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
#        else:
#            exec('from . import ' + format(item.split('.')[0]))
        

from . import models
from . import GDL
from . import Boundary
from . import Generic
from . import Sandstone
