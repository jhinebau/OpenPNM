"""
module __Cubic2__: Generate cubic networks based on Delaunay Tessellations
===============================================================================

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM
import OpenPNM.Utilities.misc as misc
import sys
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.ndimage as spim
from OpenPNM.Network.__GenericNetwork__ import GenericNetwork

class Cubic2(GenericNetwork):
    r"""

    """

    def __init__(self,**kwargs):
        '''
        Create Cubic2 network object
        '''
        super(Cubic2,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        
    def generate(self,**params):
        r'''
        Method to trigger the generation of the network
        
        Parameters
        ----------
        domain_size : list of floats, [Lx,Ly,Lz]
            Bounding cube for internal pore positions
        num_pores : int
            Number of pores to place randomly within domain

        '''
        self._logger.info(sys._getframe().f_code.co_name+": Start of network topology generation")
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
#        self._add_boundaries()
#        self._add_labels()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")

    def _generate_setup(self, domains):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")
        self._domains = domains

    def _generate_pores(self):
        r"""
        Generate the pores
        """
        domains = self._domains
        pts = sp.ndarray((0,3))
        labels = {}
        ind_t = -1
        for item in domains.keys():
            Lc = domains[item]['Lc']
            [Lx,Ly,Lz] = sp.array(domains[item]['region'][1]) - sp.array(domains[item]['region'][0])
            [Nx,Ny,Nz] = sp.rint(sp.array([Lx,Ly,Lz])/Lc).astype(dtype=int)
            ind = sp.arange(Nx*Ny*Nz)
            temp = Lc/2+Lc*sp.array(sp.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=float).T
            temp = temp + sp.array(domains[item]['region'][0])
            pts = sp.vstack((pts,temp))
            labels[item] = ind_t + 1 + ind
            ind_t = sp.amax(ind) + ind_t + 1
        self['pore.coords'] = pts

        #Convert labels
        self['pore.all'] = np.ones_like(pts[:,0],dtype=bool)
        for item in labels:
            self['pore.'+item] = sp.zeros_like(self['pore.all'],dtype=bool)
            self['pore.'+item][labels[item]] = True



    def _generate_throats(self):
        r"""
        Generate the throats connections
        """
        domains = self._domains
        self._logger.info(sys._getframe().f_code.co_name+": Define connections between pores")
        Np = self.num_pores()
        pts = self['pore.coords']
        #Perform tessellation
        self._logger.debug(sys._getframe().f_code.co_name+": Beginning tessellation")
        Tri = sptl.Delaunay(pts)
        self._logger.debug(sys._getframe().f_code.co_name+": Converting tessellation to adjacency matrix")
        Np = sp.shape(pts)[0]
        a = Tri.simplices
        b = []
        b.append(sp.vstack((a[:,0],a[:,1])).T)
        b.append(sp.vstack((a[:,0],a[:,2])).T)
        b.append(sp.vstack((a[:,0],a[:,3])).T)
        b.append(sp.vstack((a[:,1],a[:,0])).T)
        b.append(sp.vstack((a[:,1],a[:,2])).T)
        b.append(sp.vstack((a[:,1],a[:,3])).T)
        b.append(sp.vstack((a[:,2],a[:,0])).T)
        b.append(sp.vstack((a[:,2],a[:,1])).T)
        b.append(sp.vstack((a[:,2],a[:,3])).T)
        b.append(sp.vstack((a[:,3],a[:,0])).T)
        b.append(sp.vstack((a[:,3],a[:,1])).T)
        b.append(sp.vstack((a[:,3],a[:,2])).T)
        c = sp.vstack((b))
        d = sp.ones_like(c[:,0])
        adjmat = sprs.coo.coo_matrix((d,(c[:,0],c[:,1])),(Np,Np))
        #Convert to csr and back to remove duplicate entries
        adjmat = adjmat.tocsr()
        adjmat = adjmat.tocoo()
        adjmat = sprs.triu(adjmat,k=1,format="coo")
        self._logger.debug(sys._getframe().f_code.co_name+": Conversion to adjacency matrix complete")
        
        #Enter current throats into network
        conns = sp.vstack((adjmat.row, adjmat.col)).T
        self['throat.conns'] = conns 
        self._logger.debug(sys._getframe().f_code.co_name+": Propigating pore labels to new throats")
        self['throat.all'] = np.ones_like(conns[:,0],dtype=bool)
        labels = self.labels('pore')
        labels.remove('pore.all')
        for item in labels:
            self.propagate_labels(pores=self.pores(item))
        
        self._logger.debug(sys._getframe().f_code.co_name+": Identify and trim non-cubic connections")
        dist = misc.dist(self['pore.coords'],self['pore.coords'])
        tlen = dist[self['throat.conns'][:,0],self['throat.conns'][:,1]]
        bad_Ts = sp.zeros_like(self.throats(),dtype=bool)
        for item in domains.keys():
            t = self.throats(item)
            temp = tlen[t]>(1.1*domains[item]['Lc'])
            bad_Ts[t[temp]] = True
        
        #Find throats connecting domains before trimming them
        Ts1 = self.throats(labels=['CLa','GDLa'],mode='intersection')
        Ps1 = self.find_connected_pores(throats=Ts1)
        Ts2 = self.throats(labels=['CLc','GDLc'],mode='intersection')
        Ps2 = self.find_connected_pores(throats=Ts2)
        self.trim(throats=bad_Ts)
        
        self._logger.debug(sys._getframe().f_code.co_name+": Replace connections between domains")
        #Now re-add throats that are not too long
        Ds1 = misc.dist(self['pore.coords'][Ps1[:,0]],self['pore.coords'][Ps1[:,1]])
        Ds2 = misc.dist(self['pore.coords'][Ps2[:,0]],self['pore.coords'][Ps2[:,1]])
        
        #Remap DsN onto PsN
        thresh = domains['GDLa']['Lc']/1.5
        tpore1 = Ps1[sp.where(Ds1<thresh)[0],0]
        tpore2 = Ps1[sp.where(Ds1<thresh)[1],1]
        conns = sp.vstack((tpore1,tpore2)).T
        self.extend(throat_conns=conns)
        #Remap DsN onto PsN
        tpore1 = Ps2[sp.where(Ds2<thresh)[0],0]
        tpore2 = Ps2[sp.where(Ds2<thresh)[1],1]
        conns = sp.vstack((tpore1,tpore2)).T
        self.extend(throat_conns=conns)
        
        self._logger.debug(sys._getframe().f_code.co_name+": End of method")
        
    def _add_labels(self):
        r'''
        This finds surface pors simply by proximity to the domain boundaries.
        A better approach is necessary 
        '''
        coords = self.get_pore_data(prop='coords')
        self.set_pore_info(label='front',locations=(coords[:,0]<(0.1*self._Lx)))
        self.set_pore_info(label='back',locations=(coords[:,0]>(0.9*self._Lx)))
        self.set_pore_info(label='left',locations=(coords[:,1]<(0.1*self._Ly)))
        self.set_pore_info(label='right',locations=(coords[:,1]>(0.9*self._Ly)))
        self.set_pore_info(label='bottom',locations=(coords[:,2]<(0.1*self._Lz)))
        self.set_pore_info(label='top',locations=(coords[:,2]>(0.9*self._Lz)))
        bnds = self.get_pore_indices(labels=['front','back','left','right','bottom','top'])
        self.set_pore_info(label='boundary',locations=bnds)
        self.set_pore_info(label='internal',locations='all')
        
    def _add_boundaries(self):
        r"""
        This is an alternative means of adding boundaries
        """
        self._logger.info("add_boundaries: start of method")

        import scipy.spatial as sptl
        import scipy.sparse as sprs
        Lx = self._Lx
        Ly = self._Ly
        Lz = self._Lz
        Np = self.num_pores()
        btype = self._btype
        boffset = 0.05

        #Translate internal pores to each face of domain
        poffset = np.zeros((7,3))
        poffset[[2,5],0] = [-Lx, Lx]
        poffset[[3,4],1] = [-Ly, Ly]
        poffset[[1,6],2] = [-Lz, Lz]
        pcoords = pcoords0 = self.get_pore_data(prop='coords')
        for i in np.r_[1:7]:
            pcoords = np.concatenate((pcoords,pcoords0 + poffset[i,:]),axis=0)

        #Use some twisted logic to get bval list of + for boundary and - for periodic faces
        bval = [0, 1, 2, 3, 4, 5, 6]*(np.array([0, btype[2], btype[0], btype[1], btype[1], btype[0], btype[2]])*-2+1)
        ptype = np.zeros((Np,),dtype=int)
        for i in np.r_[1:7]:
            ptype = np.concatenate((ptype,np.ones((Np,),dtype=int)*bval[i]),axis=0)

        #pnum contains the internal ID number of the boundary pores (for connecting periodic points)
        pnum = self.get_pore_indices()
        pnum = np.tile(pnum,7)

        Tri = sptl.Delaunay(pcoords)
        adjmat = sprs.lil_matrix((np.shape(pcoords)[0],np.shape(pcoords)[0]),dtype=int)
        for i in np.arange(0,np.shape(Tri.simplices)[0]):
            #Keep only simplices that are fully in real domain
            adjmat[Tri.simplices[i],Tri.simplices[i]] = 1
        adjmat = sprs.triu(adjmat,k=1,format="lil")
        for i in np.arange(0,Np):
            #Add periodic throats to the netowrk (if any)
            tpore2 = pnum[adjmat.rows[i]][ptype[adjmat.rows[i]]<0]
            tpore1 = np.ones_like(tpore2,dtype=int)*i
            conns = self.get_throat_data(prop='conns')
            conns = np.concatenate((conns,np.vstack((tpore1,tpore2)).T),axis=0)
            #Add boundary pores and throats to the network
            newporetyps = np.unique(ptype[adjmat.rows[i]][ptype[adjmat.rows[i]]>0])
            newporenums = np.r_[self.num_pores():self.num_pores()+np.size(newporetyps)]
            tpore2 = newporenums
            tpore1 = np.ones_like(tpore2,dtype=int)*i
            conns = np.concatenate((conns,np.vstack((tpore1,tpore2)).T),axis=0)
            self.set_throat_data(prop='conns',data=conns)
            bcoords = np.zeros((7,3),dtype=float)
            coords = self.get_pore_data(prop='coords')
            bcoords[1,:] = [coords[i,0], coords[i,1], 0-Lz*boffset]
            bcoords[2,:] = [0-Lx*boffset, coords[i,1], coords[i,2]]
            bcoords[3,:] = [coords[i,0], -Ly*boffset, coords[i,2]]
            bcoords[4,:] = [coords[i,0], Ly+Ly*boffset, coords[i,2]]
            bcoords[5,:] = [Lx+Lx*boffset, coords[i,1], coords[i,2]]
            bcoords[6,:] = [coords[i,0], coords[i,1], Lz+Lz*boffset]
            newporecoords = bcoords[newporetyps,:]
            coords = np.concatenate((coords,newporecoords),axis=0)
            self.set_pore_data(prop='coords',data=coords)
        #Reset number of pores and throats (easier than tracking it)
        nums = np.r_[0:np.shape(coords)[0]]
        self.set_pore_data(prop='numbering',data=nums)
        self.set_pore_info(label='numbering',locations=np.ones((nums[-1]+1,),dtype=bool))
        nums = np.r_[0:np.shape(conns)[0]]
        self.set_throat_data(prop='numbering',data=nums)
        self.set_throat_info(label='numbering',locations=np.ones((nums[-1]+1,),dtype=bool))
        self._logger.debug("add_boundaries: end of method")

    def _add_boundaries_old(self):
        self._logger.info("add_boundaries_old: Start of method")

        self.add_opposing_boundaries(btype=[2,5])
        self.add_opposing_boundaries(btype=[3,4])
        self.add_opposing_boundaries(btype=[1,6])

    def _add_opposing_boundaries(self,btype=[1,6]):
        r"""
        btype indicates which two boundaries are being added by type
        """
        self._logger.info("add_opposing_boundaries: start of method")

        if btype==[2,5]:
            D=0
            W=1
            H=2
        elif btype==[3,4]:
            D=1
            W=0
            H=2
        elif btype==[1,6]:
            D=2
            W=1
            H=0

        Lx = self.domain_size[D]
        Ly = self.domain_size[W]
        Lz = self.domain_size[H]
        #Rotate pore coordinates (use only internal pores)
        pnum = self._net.pore_data['numbering'][self._net.pore_data['type']==0]
        pcoords = np.zeros_like(self._net.pore_data['coords'][pnum,:])
        pcoords[:,0] = self._net.pore_data['coords'][pnum,D]
        pcoords[:,1] = self._net.pore_data['coords'][pnum,W]
        pcoords[:,2] = self._net.pore_data['coords'][pnum,H]

        #Determine dimensions of image from dimensions of domain
        f = 100 #minimum image dimension
        im_dim = [0,0,0]
        im_dim[0] = np.floor(f*Lx/np.min([Lx,Ly,Lz]))
        im_dim[1] = np.floor(f*Ly/np.min([Lx,Ly,Lz]))
        im_dim[2] = np.floor(f*Lz/np.min([Lx,Ly,Lz]))
        im_dim = np.array(im_dim,dtype=int)

        #Convert pore coordinates into image subscripts
        im_subs = np.zeros_like(pcoords,dtype=int)
        im_subs[:,0] = pcoords[:,0]*im_dim[0]/Lx
        im_subs[:,1] = pcoords[:,1]*im_dim[1]/Ly
        im_subs[:,2] = pcoords[:,2]*im_dim[2]/Lz
        #Find linear indices of each pore in the new image
        im_inds = np.ravel_multi_index((im_subs[:,0], im_subs[:,1], im_subs[:,2]), dims=(im_dim), order='F')

        #Generate 3D image of points (place pore numbers at each site for use later)
        img = np.zeros(im_dim,dtype=int)
        img[im_subs[:,0],im_subs[:,1],im_subs[:,2]] = pnum

        #Perform distance transform on points and also get 'indicies' of each point
        img_dt, ind_dt = spim.distance_transform_edt(img==0)

        #Project all* internal points to x face
        #*Note that it's possible/likely that mutliple internal points map to the same boundary point
        img_bd0 = np.zeros([im_dim[1],im_dim[2]],dtype=int)
        img_bd1 = np.zeros([im_dim[1],im_dim[2]],dtype=int)
        img_bd0[im_subs[:,1],im_subs[:,2]] = im_inds
        img_bd1[im_subs[:,1],im_subs[:,2]] = im_inds

        #Create 2D array of distance transform indices for 0 and end faces
        dt_D0 = ind_dt[0,0,:,:]*(img_bd0>0)  #0 face
        dt_D1 = ind_dt[0,-1,:,:]*(img_bd1>0) #end face

        #Create a 2D mask containing x coordinates of internal points and -1 elsewhere
        img_D0 = -np.ones([im_dim[1],im_dim[2]],dtype=int)
        img_D1 = -np.ones([im_dim[1],im_dim[2]],dtype=int)
        img_D0[im_subs[:,1],im_subs[:,2]] = im_subs[:,0]
        img_D1[im_subs[:,1],im_subs[:,2]] = im_subs[:,0]

        #Find where x value of internal points corresponds to x value of distance transform indices
        img_bd0 = (img_D0 == dt_D0)*img_bd0
        img_bd1 = (img_D1 == dt_D1)*img_bd1

        #Convert boundary sites to linear indices
        inds_bd0 = img_bd0[np.nonzero(img_bd0)]
        inds_bd1 = img_bd1[np.nonzero(img_bd1)]

        #Use linear indices to find pore ID nums
        nums_bd0 = img[np.unravel_index(inds_bd0, dims=(im_dim), order='F')]
        nums_bd1 = img[np.unravel_index(inds_bd1, dims=(im_dim), order='F')]
        nums_bd = np.append(nums_bd0,nums_bd1)
        types_bd = np.append(np.zeros_like(nums_bd0),np.ones_like(nums_bd1))

        #Add new boundary pores and throats to the network
        Np = self._net.num_pores() #Get all pores including previously added boundaries
        bp_numbering = np.r_[Np:Np+np.size(nums_bd)]
        bp_type = (types_bd==0)*btype[0] + (types_bd==1)*btype[1]
        bp_coords = np.zeros([np.size(nums_bd),3])
        bp_coords[types_bd==0,D] = np.zeros_like(nums_bd0)-.0001
        bp_coords[types_bd==0,W] = pcoords[nums_bd0,1]
        bp_coords[types_bd==0,H] = pcoords[nums_bd0,2]
        bp_coords[types_bd==1,D] = np.ones_like(nums_bd1)*Lx+0.0001
        bp_coords[types_bd==1,W] = pcoords[nums_bd1,1]
        bp_coords[types_bd==1,H] = pcoords[nums_bd1,2]
        self._net.pore_data['numbering'] = np.append(self._net.pore_data['numbering'],bp_numbering)
        self._net.pore_data['type'] = np.append(self._net.pore_data['type'],bp_type)
        self._net.pore_data['coords'] = np.concatenate((self._net.pore_data['coords'],bp_coords))
        Nt = self._net.num_throats()
        bt_numbering = np.r_[Nt:Nt+np.size(nums_bd)]
        bt_type = np.ones(np.size(nums_bd),dtype=int)*2
        bt_connections = np.zeros([np.size(nums_bd),2],dtype=int)
        bt_connections[:,0] = nums_bd
        bt_connections[:,1] = bp_numbering
        self._net.throat_data['numbering'] = np.append(self._net.throat_data['numbering'],bt_numbering)
        self._net.throat_data['type'] = np.append(self._net.throat_data['type'],bt_type)
        self._net.throat_data['conns'] =  np.concatenate((self._net.throat_data['conns'],bt_connections))

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
