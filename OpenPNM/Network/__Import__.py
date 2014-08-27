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
import OpenPNM.Utilities.misc as misc
import scipy as sp
import scipy.io as spio
import os
from .__GenericNetwork__ import GenericNetwork

class MatFile(GenericNetwork):
    r"""
    MatFile - constructs a pore network from a perfectly formatted .mat file (MATLAB)
    
    This class contains the interface definition for the construction
    of networks
    
    Parameters
    ----------

    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)
    all other parameters are in the generate() command

    
    """
    def __init__(self, **kwargs):
        
        r"""
        Initialize
        """
        super(MatFile,self).__init__(**kwargs)
        
    def generate(self,filename='standard_cubic_5x5x5.mat', path='LocalFiles', xtra_pore_data=None, xtra_throat_data=None):
        '''
        Create network from Matlab file. Returns OpenPNM.Network.GenericNetwork() 
        object. The extra data of 'type' will trigger internal and boundary pores

        Parameters
        ----------

        Critical\n
        filename : string
            filename = 'standard_cubic_5x5x5.mat' (default)\n
            Name of mat file\n
        path : string
            path='LocalFiles' (default)\n
            the location of the mat file on your computer \n
        xtra_pore_data : list of strings
            xtra_pore_data = ['type','shape','material']
            any additional props to look for in the dictionary
        xtra_throat_data : list of strings
            xtra_throat_data = ['type','shape','material']
            any additional props to look for in the dictionary

        Examples:
        ---------

        generate network using example mat file

        >>> import OpenPNM as PNM
        >>> pn=PNM.Geometry.MatFile(name='matfile')
        >>> pn.generate(filename='standard_cubic_5x5x5.mat', path='LocalFiles')
        '''
        if path == 'LocalFiles':
            long_path = os.path.abspath(__file__)
            short_path, fname = os.path.split(long_path)
            short_path, foldername = os.path.split(short_path)  
            path, foldername = os.path.split(short_path)  
            path = os.path.join(path,'LocalFiles')
        self._path = path
        filepath = os.path.join(self._path,filename)
        self._xtra_pore_data=xtra_pore_data
        self._xtra_throat_data=xtra_throat_data
        self._dictionary=spio.loadmat(filepath)
        
        self._Np=sp.size(self._dictionary['pnumbering'])
        self._Nt=sp.size(self._dictionary['tnumbering'])
        
        #Run through generation steps
        self._add_pores()
        self._add_throats()
        self._remove_disconnected_clusters()
        self._add_xtra_pore_data()
        self._add_xtra_throat_data()
        self._add_geometry()
        
    def _add_pores(self):
        Pind = sp.arange(0,self._Np)
        self['pore.all'] = sp.ones_like(Pind,dtype=bool)
        self._logger.info('Writing pore data')
        self['pore.coords']=self._dictionary['pcoords']
        
    def _add_throats(self):
        Tind = sp.arange(0,self._Nt)
        self['throat.all']=sp.ones_like(Tind,dtype=bool)
        self._logger.info('Writing throat data')
        self['throat.conns']=self._dictionary['tconnections']
        
    def _remove_disconnected_clusters(self):
        bad_pores = sp.array([],dtype=int)
        self._pore_map = self.pores()
        self._throat_map = self.throats()
        health = self.check_network_health()        
        Np = self.num_pores()
        Nt = self.num_throats()
        cluster_sizes = [sp.shape(x)[0] for x in health['disconnected_clusters']]
        acceptable_size = min([min([50,Np/2]),max(cluster_sizes)]) # 50 or less, if it's a really small network.
        #step through each cluster of pores. If its a small cluster, add it to the list
        for cluster in health['disconnected_clusters']:
            if sp.shape(cluster)[0] < acceptable_size:
                bad_pores = sp.append(bad_pores,sp.ravel(cluster))
        bad_throats = sp.unique(self.find_neighbor_throats(bad_pores))
        #Create map for pores
        if sp.shape(bad_pores)[0] > 0:
            i = 0
            self._pore_map = sp.zeros((Np-sp.shape(bad_pores)[0],),dtype=int)
            for pore in self.pores():
                if pore not in bad_pores:
                    self._pore_map[i] = pore 
                    i += 1
        #Create map for throats
        if sp.shape(bad_throats)[0] > 0:
            i = 0
            self._throat_map = sp.zeros((Nt-sp.shape(bad_throats)[0],),dtype=int)
            for throat in self.throats():
                if throat not in bad_throats:
                    self._throat_map[i] = throat                    
                    i += 1
        self.trim(pores=bad_pores)
        #Fix the pore transformer
        try:        
            if sp.shape(bad_pores)[0] > 0:
                i = 0
                old_transform = self._dictionary['pname_transform']
                self._dictionary['pname_transform'] = sp.zeros((Np-sp.shape(bad_pores)[0],),dtype=int)
                for pore in self.pores():
                    if pore not in bad_pores:
                        self._dictionary['pname_transform'][i] = old_transform[pore]
                        i += 1
        except:
            self._logger.info('Could not update pname_transform. Imported network may not have had it.')

    def _add_geometry(self):
        try: 
            boundary_pores = sp.where(self['pore.type']!=0)[0]
            boundary_throats = sp.where(self['throat.type']!=0)[0]
            add_boundaries = True
        except: 
            boundary_pores = sp.array([])
            boundary_throats = sp.array([])
            self._logger.info('No boundary pores added.')
        Ps = sp.where([pore not in boundary_pores for pore in self.pores()])[0]
        Ts = sp.where([throat not in boundary_throats for throat in self.throats()])[0]
        geom = OpenPNM.Geometry.GenericGeometry(network=self,name='internal',pores=Ps,throats=Ts)
        geom['pore.volume'] = sp.ravel(sp.array(self._dictionary['pvolume'][self._pore_map[Ps]]))
        geom['pore.diameter'] = sp.ravel(sp.array(self._dictionary['pdiameter'][self._pore_map[Ps]]))
        geom['throat.diameter'] = self._dictionary['tdiameter'][self._throat_map[Ts]]
        geom.add_model(propname='pore.area',model=OpenPNM.Geometry.models.pore_area.spherical)
        geom.add_model(propname='throat.area',model=OpenPNM.Geometry.models.throat_area.cylinder)
        
        if add_boundaries:
            boun = OpenPNM.Geometry.Boundary(network=self,pores=boundary_pores,throats=boundary_throats,name='boundary')
    
    def _add_xtra_pore_data(self):
        xpdata = self._xtra_pore_data
        if xpdata is not None:
            if type(xpdata) is type([]):
                for pdata in xpdata:
                    try:
                        self['pore.'+pdata]=self._dictionary['p'+pdata][self._pore_map]
                    except:
                        self._logger.warning('Could not add pore data: '+pdata+' to network')
            else:
                try:
                    self['pore.'+xpdata]=self._dictionary['p'+xpdata][self._pore_map]
                except:
                    self._logger.warning('Could not add pore data: '+xpdata+' to network')

    def _add_xtra_throat_data(self):
        xtdata = self._xtra_throat_data
        if xtdata is not None:
            if type(xtdata) is type([]):
                for tdata in xtdata:
                    try:
                        self['throat.'+tdata]=self._dictionary['t'+tdata][self._throat_map]
                    except:
                        self._logger.warning('Could not add throat data: '+tdata+' to network')
            else:
                try:
                    self['throat.'+xtdata]=self._dictionary['t'+xtdata][self._throat_map]
                except:
                    self._logger.warning('Could not add throat data: '+xtdata+' to network')

    def domain_length(self,face_1,face_2):
        r'''
        Calculate the distance between two faces of the network
        
        Parameters
        ----------
        face_1 and face_2 : array_like
            Lists of pores belonging to opposite faces of the network
            
        Returns
        -------
        The length of the domain in the specified direction
        
        Notes
        -----
        - Does not yet check if input faces are perpendicular to each other
        '''
        #Ensure given points are coplanar before proceeding
        if misc.iscoplanar(self['pore.coords'][face_1]) and misc.iscoplanar(self['pore.coords'][face_2]):
            #Find distance between given faces
            x = self['pore.coords'][face_1]
            y = self['pore.coords'][face_2]
            Ds = misc.dist(x,y)
            L = sp.median(sp.amin(Ds,axis=0))
        else:
            self._logger.warning('The supplied pores are not coplanar. Length will be approximate.')
            f1 = self['pore.coords'][face_1]
            f2 = self['pore.coords'][face_2]
            distavg = [0,0,0]
            distavg[0] = sp.absolute(sp.average(f1[:,0]) - sp.average(f2[:,0]))
            distavg[1] = sp.absolute(sp.average(f1[:,1]) - sp.average(f2[:,1]))
            distavg[2] = sp.absolute(sp.average(f1[:,2]) - sp.average(f2[:,2]))
            L = max(distavg)
        return L

        
    def domain_area(self,face):
        r'''
        Calculate the area of a given network face
        
        Parameters
        ----------
        face : array_like
            List of pores of pore defining the face of interest
            
        Returns
        -------
        The area of the specified face
        '''
        coords = self['pore.coords'][face]
        rads = self['pore.diameter'][face]/2.
        # calculate the area of the 3 principle faces of the bounding cuboid
        dx = max(coords[:,0]+rads) - min(coords[:,0]-rads)
        dy = max(coords[:,1]+rads) - min(coords[:,1]-rads)
        dz = max(coords[:,2]+rads) - min(coords[:,2]-rads)
        yz = dy*dz # x normal
        xz = dx*dz # y normal
        xy = dx*dy # z normal
        # find the directions parallel to the plane
        directions = sp.where([yz,xz,xy]!=max([yz,xz,xy]))[0] 
        try:
            # now, use the whole network to do the area calculation
            coords = self['pore.coords']
            rads = self['pore.diameter']/2.
            d0 = (max(coords[:,directions[0]]+rads) - min(coords[:,directions[0]]-rads))
            d1 = (max(coords[:,directions[1]]+rads) - min(coords[:,directions[1]]-rads))
            A = d0*d1        
        except:
            # if that fails, use the max face area of the bounding cuboid
            A = max([yz,xz,xy])
        if not misc.iscoplanar(self['pore.coords'][face]):
            self._logger.warning('The supplied pores are not coplanar. Area will be approximate')
        return A

        