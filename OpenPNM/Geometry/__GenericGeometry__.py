#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import scipy.sparse as sprs
import scipy.stats as spst

class GenericGeometry(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericGeometry - Base class to construct pore networks
    
    This class contains the interface definition for the construction of networks
    
    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
        
    """

    def __init__(self):
        
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__()
        self._logger.debug("Method: Constructor")
        
    def generate(self, **params):
        r"""
        Generate the network
        """
        self._logger.debug("self.generate()")
        self.generate_setup(**params)
        self.generate_pores()
        self.generate_throats()
        self._net.update()
#        self.add_boundaries()
        self.generate_pore_seeds()
        self.generate_throat_seeds()
        self.generate_pore_diameters(params['psd_info'])
        self.generate_throat_diameters(params['tsd_info'])
        self.calc_pore_volumes()
        self.calc_throat_lengths()
        self.calc_throat_volumes()
        self._net.update()
        self._logger.debug("\t end of self.generate()")
        return self._net  

    def generate_setup(self,**params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology        
        """
        self._logger.error("generation_setup: not implemented")

    def generate_pores(self):
        r"""
        Generate the pores (numbering, types and coordinates)
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology        
        """
        self._logger.error("generate_pores: not implemented")
        
    def generate_throats(self):        
        r"""
        Generate the throats (numbering and types)
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("generate_throats: not implemented")
        
    def add_boundaries(self):
        r"""
        Add boundary pores around network (numbering and types)
        
        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("add_boundaries: not implemented")        
        
    def generate_throat_type(self):
        self._logger.info("update_throat_type: Start of method")
        self._net.throat_properties["type"] = self._net._incmatrix.transpose() * self._net.pore_properties["type"]
        self._logger.debug("update_throat_type: End of method")

    def generate_pore_seeds(self):
        r"""
        
        Notes
        -----
        Assigns random seed to pores
        """
        self._logger.info("generate_pore_seeds: Assign each pore a random seed")
        Np = self._net.get_num_pores()
        self._net.pore_properties['seed'] = sp.random.rand(Np)
        #Set boundary pore to 1 (max possible value) so throats adopt seed from interal pore
        self._net.pore_properties['seed'][self._net.pore_properties['type']>0] = 1
        self._logger.debug("generate_pore_seeds: End of method")

    def generate_throat_seeds(self):
        r"""
        
        Notes
        -----
        Assigns random seeds to throats based on smaller of neighboring pores
        """
        self._logger.info("generate_throat_seeds: Assign each throat its smaller neighboring pore seed")
        self._net.throat_properties['seed'] = sp.amin(self._net.pore_properties['seed'][self._net.throat_properties['connections']],1)
        self._logger.debug("generate_throat_seeds: End of method")
        
    def generate_pore_diameters(self,psd_info):
        r"""
        Calculates pore diameter from given statisical distribution using the random seeds provided by generate_pore_seeds()
        
        Notes
        -----
                
        """
        self._logger.info("generate_pore_diameters: Generate pore diameter from "+psd_info['name']+" distribution")
        prob_fn = getattr(spst,psd_info['name'])
        P = prob_fn(psd_info['shape'],loc=psd_info['loc'],scale=psd_info['scale'])
        self._net.pore_properties['diameter'] = P.ppf(self._net.pore_properties['seed'])
        #Set boundadry pores to size 0
        self._net.pore_properties['diameter'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("generate_pore_diameters: End of method")

    def generate_throat_diameters(self,tsd_info):
        r"""
        
        Notes
        -----
        Calculates throat diameter from given statisical distribution using the random seeds provided by generate_throat_seeds()
        """
        self._logger.info("generate_throat_diameters: Generate throat diameter from "+tsd_info['name']+" distribution")
        prob_fn = getattr(spst,tsd_info['name'])
        P = prob_fn(tsd_info['shape'],loc=tsd_info['loc'],scale=tsd_info['scale'])
        self._net.throat_properties['diameter'] = P.ppf(self._net.throat_properties['seed'])
        self._logger.debug("generate_throat_diameters: End of method")
        
    def calc_pore_volumes(self):
        r"""
        Calculates pore volume from diameter assuming a spherical pore
        """
        self._logger.info("calc_pore_volumes: Setting pore volumes assuming cubic bodies")
        #Set internal pore volumes to 1
        self._net.pore_properties['volume'] = self._net.pore_properties['diameter']**3
        #Set boundary pore volumes to 0
        self._net.pore_properties['volume'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("calc_pore_volumes: End of method")
        
    def calc_throat_volumes(self):
        r"""
        Calculates throat volume from diameter and length assuming a cylindrical pore of constant cross-section
        
        Notes
        -----
        
        """
        self._logger.info("calc_throat_volumes: Setting throat volumes assuming square cross-section")
        #Set internal pore volumes to 1
        self._net.throat_properties['volume'] = self._net.throat_properties['length']*self._net.throat_properties['diameter']**2
        self._logger.debug("calc_throat_volumes: End of method")
        
    def calc_throat_lengths(self):
        r"""
        Determine throat length from distance between pores
        
        Notes
        -----
        
        """
        self._logger.info("calc_throat_lengths: Determine throat length from distance between pores")
        #Initialize throat_property['length']
        self._net.throat_properties['length'] = sp.zeros_like(self._net.throat_properties['type'])
        C1 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,0]]
        C2 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,1]]
        E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
        D1 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,0]]
        D2 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,1]]
        self._net.throat_properties['length'] = E - (D1 + D2)/2
        #Perform check for unphysical throat lengths
        if sp.sum(self._net.throat_properties['length']<0):
            self._logger.warning("calc_throat_lengths: Some negative throat lengths exist, some pores overlap!")
        self._logger.debug("calc_throat_lengths: End of method")     
        
    def translate_coordinates(self,net,displacement=[0,0,0]):
        net.pore_properties['coords'] = net.pore_properties['coords'] + displacement
        
    def scale_coordinates(self,net,scale=[1,1,1]):
        net.pore_properties['coords'] = net.pore_properties['coords']*scale 
        
    def get_net(self):
        r"""
        Return the produced network
        """
        return self._net

if __name__ == '__main__':
    test=GenericGeometry(loggername="TestGenerator")        

