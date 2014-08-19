"""
module __Core__: Base class to construct pore network tools
==================================================================

.. warning:: The classes of this module should be loaded through the 'Base.__init__.py' file.

"""

import pprint,collections
from functools import partial
import scipy as sp
import OpenPNM
from OpenPNM.Base import Base
from OpenPNM.Utilities import misc


class Core(Base):
    r'''
    This is the core class from which OpenPNM objects derive.  
    '''
    def __init__(self, **kwargs):
        r'''
        Initialize
        '''
        super(Core,self).__init__(**kwargs)
        self._logger.info("Construct Core subclass from Base")

        #Initialize ordered dict for storing property models
        self._models = collections.OrderedDict()

        self._logger.debug("Construction of Core class complete")
        
    def __setitem__(self,key,value):
        r'''
        This is a subclass of the default __setitem__ behavior.  The main aim
        is to limit what type and shape of data can be written to protect
        the integrity of the network.
        '''
        #Enforce correct dict naming
        element = key.split('.')[0]
        if (element != 'pore') and (element != 'throat'):
            self._logger.error('Array name must begin with \'pore\' or \'throat\'')
            return
        #Convert value to an ndarray
        value = sp.array(value,ndmin=1)
        #Skip checks for 'coords', 'conns'
        if (key == 'pore.coords') or (key == 'throat.conns'):
            super(Base, self).__setitem__(key,value)
            return
        #Skip checks for protected props, and prevent changes if defined
        if key.split('.')[1] in ['all','map']:
            if key in self.keys():
                self._logger.error(key+' is already defined.')
            else:
                self._logger.debug(key+' is being defined.')
                super(Base, self).__setitem__(key,value)
            return
        #Write value to dictionary
        if sp.shape(value)[0] == 1:  # If value is scalar
            self._logger.debug('Broadcasting scalar value into vector: '+key)
            value = sp.ones((self._count(element),),dtype=value.dtype)*value
            super(Base, self).__setitem__(key,value)
        elif sp.shape(value)[0] == self._count(element):
            self._logger.debug('Updating vector: '+key)
            super(Base, self).__setitem__(key,value)
        else:
            self._logger.error('Cannot write vector with an array of the wrong length: '+key)
            
    def add_model(self,propname,model,regen_mode='static',**kwargs):
        r'''
        Add specified property estimation model to the object.
        
        Parameters
        ----------
        propname : string
            The name of the property to use as dictionary key, such as
            'pore.diameter' or 'throat.length'
            
        model : function
            The property estimation function to use
            
        regen_mode : string
            Controls when and if the property is regenerated. Options are:
            
            * 'static' : The property is stored as static data and is only regenerated when the object's `regenerate` is called
            * 'constant' : The property is calculated once when this method is first run, but always maintains the same value

        Notes
        -----
        This method is inherited by all net/geom/phys/phase objects.  It takes
        the received model and stores it on the object under private dictionary
        called _models.  This dict is an 'OrderedDict', so that the models can
        be run in the same order they are added.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn)
        >>> import OpenPNM.Geometry.models as gm
        >>> f = gm.pore_misc.random  # Get model from Geometry library
        >>> geom.add_model(propname='pore.seed',model=f)
        >>> list(geom._models.keys())  # Look in private dict to verify model was added
        ['pore.seed']
        
        '''
        #Determine object type, and assign associated objects
        self_type = self.__module__.split('.')[1]
        network = self._net
        phase = None
        geometry = None
        physics = None
        if self_type == 'Geometry':
            geometry = self
        elif self_type == 'Phases':
            phase = self
        elif self_type == 'Physics':
            phase = self._phases[0]
            physics = self
        #Build partial function from given kwargs
        fn = partial(model,network=network,phase=phase,geometry=geometry,physics=physics,**kwargs)
        if regen_mode == 'static':
            self[propname] = fn()  # Generate data and store it locally
            self._models[propname] = fn  # Store model in a private attribute
        if regen_mode == 'constant':
             self[propname] = fn()  # Generate data and store it locally
             
    def remove_model(self,propname):
        r'''
        Remove pore scale property models from current object.
        
        Parameters
        ----------
        propname : string
            The property name for which the model should be removed
            
        Notes
        -----
        This does not delete the data associated with the model, only the model.
        With this model removed it is now possible to write constants to the
        propname dictionary and they won't be overwritten when the object is 
        regenerated.
        '''
        self._models.pop(propname,None)
        
    def regenerate(self, props='',mode='inclusive'):
        r'''
        This updates properties using any models on the object that were 
        assigned using ``add_model``

        Parameters
        ----------
        props : string or list of strings
            The names of the properties that should be updated, defaults to 'all'
        mode : string
            This controls which props are regenerated and how.  Options are:
            
            * 'inclusive': (default) This regenerates all given properties
            * 'exclude': This generates all given properties EXCEPT the given ones
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
        >>> import OpenPNM.Geometry.models as gm  # Import Geometry model library
        >>> f = gm.pore_misc.random  # Get random seed generator model
        >>> geom.add_model(propname='pore.seed',model=f,seed=None)  # Add model to Geometry object
        >>> geom['pore.seed'][0]  # Look at seed value in pore 0
        ...
        >>> geom.regenerate()  # Regenerate all models
        >>> geom['pore.seed'][0]  # Look at same seed value again
        ...
        
        '''
        if props == '':
            props = self._models.keys()
        elif type(props) == str:
            props = [props]
        if mode == 'exclude':
            temp = list(self._models.keys())
            for item in props:
                temp.remove(item)
            props = temp
        for item in props:
            if item in self._models.keys():
                self[item] = self._models[item]()
            else:
                self._logger.warning('Requested proptery is not a dynamic model: '+item)
    

    #--------------------------------------------------------------------------
    '''Setter and Getter Methods'''
    #--------------------------------------------------------------------------
    def set_data(self,pores=None,throats=None,prop='',data='',mode='merge'):
        r'''
        Write data according to input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve.
        pores : array_like, or string 'all'
            List of pore indices to which data will be written.
        throats : array_like, or string 'all'
            List of throat indices to which data will be written.

        data: array_like
        
        mode : string
            Set the mode to be used for writing data to the specified indices.  Options are:
            
            * 'merge' : (default) Adds data to specified locations while 
            maintaining pre-existing data
            
            * 'overwrite' : Adds data to specified locations while 
            removing all pre-existing data
            
            * 'remove_data' : Removes data from specified locations.
            
            * 'remove_prop' : Removes property key from the dictionary of the object

        '''        

        if mode in ['merge','overwrite','remove_data','remove_prop']:
            #Checking for the indices
            data= sp.array(data,ndmin=1)
            if pores != None:
                if pores == 'all':  pores = self.pores()
                else:   pores = sp.array(pores,ndmin=1)
                element='pore'
                locations = pores
            elif throats != None:
                if throats == 'all':    throats = self.throats()
                else:   throats = sp.array(throats,ndmin=1)
                element='throat'
                locations = throats 
            elif throats==None and pores==None:
                self._logger.error('No pores or throats have been received!')
                raise Exception()
            if prop.split('.')[0] == element:
                prop = prop.split('.')[1]
            if mode=='remove_data':
                if data=='':
                    try: 
                        self[element+'.'+prop]   
                        self[element+'.'+prop][locations] = sp.nan
                        self._logger.debug('For the '+element+' property '+prop+', the specified data have been deleted in '+self.name)                        
                    except: self._logger.error(element+' property '+prop+' in '+self.name+' has not been defined. Therefore, no data can be removed!')                              
                else:   
                    self._logger.error('For the '+element+' property '+prop+' in '+self.name+': The (remove_data) mode, will remove data from the specified locations. No data should be sent!')        
                    raise Exception()             
            elif mode=='remove_prop':
                if data=='':
                    del self[element+'.'+prop]
                    self._logger.debug(element+' property '+prop+' has been deleted from the dictionary in '+self.name)
                else:   
                    self._logger.error('For the property '+prop+' in '+self.name+': The (remove_prop) mode, will remove the property from the dictionary. No data should be sent!')        
                    raise Exception() 
            else:
                if data.ndim > 1: data = data.squeeze()
                try:
                    self[element+'.'+prop]
                    temp_word = 'updated for '
                except: temp_word = 'added to '            
                if sp.shape(data)[0]==1 or sp.shape(locations)[0]==sp.shape(data)[0]:          
                    try: 
                        self[element+'.'+prop]
                        if mode=='overwrite':
                            self[element+'.'+prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                    except: self[element+'.'+prop] = sp.zeros((getattr(self,'num_'+element+'s')(),))*sp.nan
                    self[element+'.'+prop][locations] = data
                    self._logger.debug(element+' property '+prop+' has been '+temp_word+self.name)
    
                else:         
                    self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', data can be scalar or should have the same size as the locations!')
                    raise Exception()
        else:  self._logger.error('For adding '+element+' property '+prop+' to '+self.name+', the mode '+mode+' cannot be applied!')

    
    def get_data(self,pores=None,throats=None,prop='',mode=''):
        r'''
        Retrieves data from the object to which it is associated according to 
        input arguments.
        
        Parameters
        ----------
        prop : string
            Name of property to retrieve
        pores : array_like, or string 'all'
            List of pore indices from which to retrieve data.  If set to 'all',
            then ALL values are returned
        throats : array_like, or string 'all'
            List of throat indies from which to retrieve data.  If set to 'all'
            , then ALL values are returned.
        mode : string
            Set the mode to be used for getting data from the specified indices.  Options are:
            
            * '' : (default) gets data from the specified locations.
            
            * 'interpolate': Determines a pore (or throat) property as the average of it's neighboring 
            throats (or pores)

        '''   

        if mode in ['','interpolate']:
            if pores != None:
                if pores == 'all':  pores = self.pores()
                else:   pores = sp.array(pores,ndmin=1)
                element='pore'
                locations = pores
            elif throats != None:
                if throats == 'all':    throats = self.throats()
                else:   throats = sp.array(throats,ndmin=1)
                element='throat'
                locations = throats 
            elif throats==None and pores==None:
                self._logger.error('No pores or throats have been received!')
                raise Exception('No pores/throats have been sent!') 
            if prop.split('.')[0] in ['pore','throat']:
                if prop.split('.')[0]==element:
                    prop = prop.split('.')[1]
                else:
                    raise Exception('Wrong property has been sent for the locations! Pore property can only be returned for pore indices. The same applies to the throats.')                    
            
            if mode=='interpolate':
                if element == 'pore':
                    try:
                        self['throat.'+prop]
                    except:
                        self._logger.error('For getting pore property using interpolate mode in '+self.name+' : throat property '+prop+' has not been found!')
                        raise Exception('throat data for the property not found!')
                    return self.interpolate_data(self['throat.'+prop])[locations]
                elif element=='throat':
                    try:
                        self['pore.'+prop]
                    except:
                        self._logger.error('For getting throat property using interpolate mode in '+self.name+' : pore property '+prop+' has not been found!')
                        raise Exception('pore data for the property not found!')
                    return self.interpolate_data(self['pore.'+prop])[locations]
            else:
                try:    self[element+'.'+prop]
                except: 
                    self._logger.error(element+' property '+prop+' does not exist for '+self.name)
                    raise Exception('No data found for the property!')
                try:    return self[element+'.'+prop][locations]
                except: self._logger.error('In '+self.name+', data for '+element+' property '+prop+' in these locations cannot be returned.')
 
        else: self._logger.error('For getting property: '+prop+' in '+self.name+', the mode '+mode+' cannot be applied!')    
 
    

    def set_info(self,pores=None,throats=None,label='',mode='merge'):
        r'''
        Apply a label to a selection of pores or throats
        
        Parameters
        ----------
        label : string
            The name of the pore labels you wish to apply (e.g. 'top')
        pores or throats : array_like
            An array containing the pore (or throat) indices where the labels 
            should be applied.  Can be either a boolean mask of Np length with 
            True at labels locations (default), a list of indices where labels 
            should be applied.
        mode : string
            Set the mode to be used for writing labels.  Options are:
            
            * 'merge' : (default) Adds label to specified locations while 
            maintaining pre-existing labels
            
            * 'overwrite' : Adds label to specified locations while 
            removing all pre-existing labels
            
            * 'remove_info' : Removes label from the specified locations.
            
            * 'remove_label' : Removes the entire label from the dictionary of the object
        '''
        if mode in ['merge','overwrite','remove_info','remove_label']:
            if pores != None:
                if pores == 'all':  pores = self.pores()
                else:   pores = sp.array(pores,ndmin=1)
                element='pore'
                locations = pores
            elif throats != None:
                if throats == 'all':    throats = self.throats()
                else:   throats = sp.array(throats,ndmin=1)
                element='throat'
                locations = throats 
            elif throats==None and pores==None:
                self._logger.error('No pores or throats have been received!')
                raise Exception()
            locations = sp.array(locations,ndmin=1)
            if label:
                if label.split('.')[0] == element:
                    label = label.split('.')[1]
                if label=='all':
                    try: 
                        old_label = self[element+'.'+label]
                        if sp.shape(old_label)[0]<sp.shape(locations)[0]:
                            self[element+'.'+label] = sp.ones_like(locations,dtype=bool)
                            self._logger.info('label=all for '+element+'has been updated to a bigger size!')
                            for info_labels in getattr(self,'labels')():
                                if info_labels.split('.')[0] == element:
                                    info_labels = info_labels.split('.')[1]
                                if info_labels!=label:
                                    temp = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                                    temp[old_label] = self[element+'.'+info_labels]
                                    self[element+'.'+info_labels] = temp
                        elif sp.shape(old_label)[0]>sp.shape(locations)[0]: 
                            self._logger.error('To apply a new numbering label (label=all) to '+element+'s, size of the locations cannot be less than total number of '+element+'s!!')
                    except: self[element+'.'+label] = sp.ones_like(locations,dtype=bool)
                else:    
                    try: self[element+'.'+label]
                    except: self[element+'.'+label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                    if mode=='overwrite':
                        self[element+'.'+label] = sp.zeros((getattr(self,'num_'+element+'s')(),),dtype=bool)
                        self[element+'.'+label][locations] = True
                    elif mode=='remove_info':  self[element+'.'+label][locations] = False
                    elif mode=='merge':  self[element+'.'+label][locations] = True
                    elif mode=='remove_label':  del self[element+'.'+label]
            else: self._logger.error('No label has been defined for these locations.')                
    
        else: self._logger.error('For setting '+element+' '+label+' to '+self.name+', the mode '+mode+' cannot be applied!')    

    
    def get_info(self,pores=None,throats=None,label=''):
        r'''
        Retrieves the locations where the specified label is applied
        
        Parameters
        ----------
        label : string
            Label of interest
        pores (or throats) : array_like
            List of pores or throats
            
        See Also
        --------
        `get_labels`
            
        '''
        #Clean up label argument to remove leading 'pore' or 'throat'
        if pores != None:
            if pores == 'all':  pores = self.pores()
            else:   pores = sp.array(pores,ndmin=1)
            element='pore'
            locations = pores 
        elif throats != None:
            if throats == 'all':    throats = self.throats()
            else:   throats = sp.array(throats,ndmin=1)
            element='throat'
            locations = throats 

        elif throats==None and pores==None:
            self._logger.error('No pores or throats have been received!')
            raise Exception()
        if label.split('.')[0] == element:
            label = label.split('.')[1]
            
        temp = self[element+'.'+label]
        return temp[locations]
            

           
    def amalgamate_data(self,objs=[]):
        r"""
        Returns a dictionary containing ALL pore data from all netowrk and/or
        phase objects received as arguments
        """
        if type(objs) != list:
            objs = list(objs)
        data_amalgamated = {}
        exclusion_list = ['pore.centroid','pore.vertices','throat.centroid','throat.offset_verts','throat.verts','throat.normals','throat.perimeter']
        for item in objs:
#            try:
            if item.__module__[0:15] == 'OpenPNM.Network': #if network object, combine geometry and network keys
                keys = []
                for key in item.keys():
                    keys.append(key)
                for geom in item._geometries:
                    for key in geom.keys():
                        if key not in keys:
                            keys.append(key)
            else:
                if item.__module__[0:14] == 'OpenPNM.Phases':
                    keys = []
                    for key in item.keys():
                        keys.append(key)
                    for physics in item._physics:
                        for key in physics.keys():
                            if key not in keys:
                                keys.append(key)
            keys.sort()
            for key in keys:
                if key not in exclusion_list:
                    if sp.amax(item[key]) < sp.inf:
                        dict_name = item.name+'.'+key
                        data_amalgamated.update({dict_name : item[key]})
#            except: 
#                self._logger.error('Only network and phase items contain data')
        return data_amalgamated
            
    def props(self,element='',mode='all'):
        r'''
        Returns a list containing the names of all defined pore or throat
        properties. 
        
        Parameters
        ----------
        element : string, optional
            Can be either 'pore' or 'throat' to specify what properties are
            returned.  If no element is given, both are returned
            
        mode : string, optional
            Controls what type of properties are returned.  Options are:
            
            - 'all' : Returns all properties on the object
            - 'models' : Returns only properties that are associated with a model 
            - 'constants' : Returns only properties that are set as constant values

        Returns
        -------
        A an alphabetically sorted list containing the string name of all 
        pore or throat properties currently defined.  This list is an iterable,
        so is useful for scanning through properties.

        See Also
        --------
        `labels`
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.props()
        ['pore.coords', 'throat.conns']
        >>> pn.props('pore')
        ['pore.coords']
        '''
        
        props = []
        for item in self.keys():
            if self[item].dtype != bool:
                props.append(item)
        models = list(self._models.keys())
        constants = props.copy()
        for item in models:
            constants.remove(item)
        
        if element in ['pore','pores']:
            element = 'pore'
        elif element in ['throat','throats']:
            element = 'throat'
            
        temp = []
        if mode == 'all':
            if element == '': temp = props
            else: temp = [item for item in props if item.split('.')[0]==element]
        elif mode == 'models':
            if element == '': temp = models
            else: temp = [item for item in models if item.split('.')[0]==element]
        elif mode == 'constants':
            if element == '': temp = constants
            else: temp = [item for item in constants if item.split('.')[0]==element]
        return misc.PrintableList(temp)

            
    def _get_labels(self,element='',locations=[],mode='union'):
        r'''
        This is the actual label getter method, but it should not be called directly.  
        Wrapper methods have been created, use get_labels().
        '''
        labels = []
        for item in self.keys():
            if item.split('.')[0] == element:
                if self[item].dtype == bool:
                    labels.append(item)
        labels.sort()
        if locations == []:
            return misc.PrintableList(labels)
        else:
            labels = sp.array(labels)
            arr = sp.zeros((sp.shape(locations)[0],len(labels)),dtype=bool)
            col = 0
            for item in labels:
                arr[:,col] = self[item][locations]
                col = col + 1
            if mode == 'count':
                return sp.sum(arr,axis=1)
            if mode == 'union':
                temp = labels[sp.sum(arr,axis=0)>0]
                temp.tolist()
                return misc.PrintableList(temp)
            if mode == 'intersection':
                temp = labels[sp.sum(arr,axis=0)==sp.shape(locations,)[0]]
                temp.tolist()
                return misc.PrintableList(temp)
            if mode == 'difference':
                temp = labels[sp.sum(arr,axis=0)!=sp.shape(locations,)[0]]
                temp.tolist()
                return misc.PrintableList(temp)
            if mode == 'mask':
                return arr
            if mode == 'none':
                temp = sp.ndarray((sp.shape(locations,)[0],),dtype=object)
                for i in sp.arange(0,sp.shape(locations,)[0]):
                    temp[i] = list(labels[arr[i,:]])
                return temp
            else:
                print('unrecognized mode')
                
    def labels(self,element='',pores=[],throats=[],mode='union'):
        r'''
        Returns the labels applied to specified pore or throat locations
        
        Parameters
        ----------
        pores (or throats) : array_like
            The pores (or throats) whose labels are sought.  If left empty a 
            dictionary containing all pore and throat labels is returned.
        mode : string, optional
            Controls how the query should be performed
            
            * 'none' : An N x Li list of all labels applied to each input pore (or throats). Li can vary betwen pores (and throats)
            
            * 'union' : A list of labels applied to ANY of the given pores (or throats)
            
            * 'intersection' : Label applied to ALL of the given pores (or throats)
            
            * 'difference' : Labels NOT applied to ALL pores (or throats)
            
            * 'count' : The number of labels on each pores (or throats)
            
            * 'mask' : returns an N x Lt array, where each row corresponds to a pore (or throat) location, and each column contains the truth value for the existance of labels as returned from labels(pores='all',mode='union')).
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.labels(pores=[0,1,5,6])
        ['pore.all', 'pore.bottom', 'pore.front', 'pore.left']
        >>> pn.labels(pores=[0,1,5,6],mode='intersection')
        ['pore.all', 'pore.bottom']
        '''
        if (pores == []) and (throats == []):
            if element == '':
                temp = []
                temp = self._get_labels(element='pore')
                temp = temp + self._get_labels(element='throat')
            elif element in ['pore','pores']:
                temp = self._get_labels(element='pore',locations=[], mode=mode)
            elif element in ['throat','throats']:
                temp = self._get_labels(element='throat',locations=[], mode=mode)
            else:
                self._logger.error('Unrecognized element')
                return
        elif pores != []:
            if pores == 'all':
                pores = self.pores()
            pores = sp.array(pores,ndmin=1)
            temp = self._get_labels(element='pore',locations=pores, mode=mode)
        elif throats != []:
            if throats == 'all':
                throats = self.throats()
            throats = sp.array(throats,ndmin=1)
            temp = self._get_labels(element='throat',locations=throats,mode=mode)
        return temp
            
    def filter_by_label(self,pores=[],throats=[],label=''):
        r'''
        Returns which of the supplied pores (or throats) has the specified label
        
        Parameters
        ----------
        pores, or throats : array_like
            List of pores or throats to be filtered
            
        label : string
            The label to apply as a filter
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.filter_by_label(pores=[0,1,5,6],label='left')
        array([0, 1])
        '''
        if pores != []:
            label = 'pore.'+label.split('.')[-1]
            all_labels = self.labels('pore')
            mask = self.labels(pores=pores,mode='mask')
            ind = all_labels.index(label)
            temp = mask[:,ind]
            pores = sp.array(pores,ndmin=1)
            return pores[temp]
        elif throats != []:
            label = 'throat.'+label.split('.')[-1]
            all_labels = self.labels('throat')
            mask = self.labels(throats=throats,mode='mask')
            ind = all_labels.index(label)
            temp = mask[:,ind]
            throats = sp.array(throats,ndmin=1)
            return throats[temp]            
        
    def _get_indices(self,element,labels,mode):
        r'''
        This is the actual method for getting indices, but should not be called
        directly.  Use pores or throats instead.
        '''
        if element+'.all' not in self.keys():
            raise Exception('Cannot proceed without {}.all'.format(element))
        if type(labels) == str: labels = [labels] #convert string to list, if necessary
        if mode == 'union':
            union = sp.zeros_like(self[element+'.all'],dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    union = union + self[element+'.'+item.split('.')[-1]]
            ind = union
        elif mode == 'intersection':
            intersect = sp.ones_like(self[element+'.all'],dtype=bool)
            for item in labels: #iterate over labels list and collect all indices
                    intersect = intersect*self[element+'.'+item.split('.')[-1]]
            ind = intersect
        elif mode == 'not_intersection':
            not_intersect = sp.zeros_like(self[element+'.all'],dtype=int)
            for item in labels: #iterate over labels list and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                not_intersect = not_intersect + sp.int8(info)
            ind = (not_intersect == 1)
        elif mode == 'difference':
            none = sp.zeros_like(self[element+'.all'],dtype=int)
            for item in labels: #iterate over labels list and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                none = none - sp.int8(info)
            ind = (none == 0)
        #Extract indices from boolean mask
        ind = sp.where(ind==True)[0]
        return ind
        
    def pores(self,labels='all',mode='union'):
        r'''
        Returns pore locations where given labels exist.
        
        Parameters
        ----------
        labels : list of strings, optional
            The pore label(s) whose locations are requested.
            If omitted, all pore inidices are returned.
        mode : string, optional
            Specifies how the query should be performed.  The options are:    

            * 'union' : (default) All pores with ANY of the given labels are 
            returned.

            * 'intersection' : Only pore with ALL the given labels are 
            returned.

            * 'not_intersection' : Only pores with exactly one of the given 
            labels are returned.
            
            * 'difference' : Only pores with none of the given labels are returned.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pind = pn.get_pore_indices(labels=['top','front'],mode='union')
        >>> pind[[0,1,2,-3,-2,-1]]
        array([  0,   5,  10, 122, 123, 124], dtype=int64)
        >>> pn.get_pore_indices(labels=['top','front'],mode='intersection')
        array([100, 105, 110, 115, 120], dtype=int64)
        '''
        if labels == 'all':
            Np = sp.shape(self['pore.all'])[0]
            ind = sp.arange(0,Np)
        else:
            ind = self._get_indices(element='pore',labels=labels,mode=mode)
        return ind
        
    def throats(self,labels='all',mode='union'):
        r'''
        Returns throat locations where given labels exist.
        
        Parameters
        ----------
        labels : list of strings, optional
            The throat label(s) whose locations are requested.
            If omitted, 'all' throat inidices are returned.
        mode : string, optional
            Specifies how the query should be performed.  The options are: 

            * 'union' : (default) All throats with ANY of the given labels are 
            returned.

            * 'intersection' : Only throats with ALL the given labels are 
            counted.

            * 'not_intersection' : Only throats with exactly one of the given 
            labels are counted.
            
            * 'difference' : Only throats with none of the given labels are returned.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> Tind = pn.throats()
        >>> Tind[0:5]
        array([0, 1, 2, 3, 4], dtype=int64)
        
        '''
        if labels == 'all':
            Nt = sp.shape(self['throat.all'])[0]
            ind = sp.arange(0,Nt)
        else:
            ind = self._get_indices(element='throat',labels=labels,mode=mode)
        return ind
        
    def tomask(self,pores=None,throats=None):
        r'''
        Convert a list of pore or throat indices into a boolean mask of the 
        correct length
        
        Parameters
        ----------
        pores or throats : array_like
            List of pore or throat indices
            
        Returns
        -------
        mask : array_like
            A boolean mask of length Np or Nt with True in the locations of
            pores or throats received.  
        
        '''
        if pores != None:
            Np = sp.shape(self['pore.all'])[0]
            pores = sp.array(pores,ndmin=1)
            mask = sp.zeros((Np,),dtype=bool)
            mask[pores] = True
        if throats != None:
            Nt = sp.shape(self['throat.all'])[0]
            throats = sp.array(throats,ndmin=1)
            mask = sp.zeros((Nt,),dtype=bool)
            mask[throats] = True
        return mask
            
    def toindices(self,mask):
        r'''
        Convert a boolean mask a list of pore or throat indices
        
        Parameters
        ----------
        mask : array_like booleans
            A boolean array with True at locations where indices are desired.
            The appropriate indices are returned based an the length of mask, 
            which must be either Np or Nt long.  
            
        Returns
        -------
        indices : array_like
            A list of pore or throat indices corresponding the locations where
            the received mask was True.
            
        Notes
        -----
        This behavior could just as easily be accomplished by using the mask
        in pn.pores()[mask] or pn.throats()[mask].  This method is just a thin
        convenience function and is a compliment to tomask().
        
        '''
        if sp.shape(mask)[0] == self.num_pores():
            indices = self.pores()[mask]
        elif sp.shape(mask)[0] == self.num_throats():
            indices = self.throats()[mask]
        else:
            raise Exception('Mask received was neither Np nor Nt long')
        return indices
        
    def interpolate_data(self,data):
        r"""
        Determines a pore (or throat) property as the average of it's neighboring 
        throats (or pores)

        Parameters
        ----------
        data : array_like
            A list of specific values to be interpolated.  List MUST be either
            Np or Nt long
        
        Returns
        -------
        An array containing interpolated pore (or throat) data

        Notes
        -----
        - This uses an unweighted average, without attempting to account for 
        distances or sizes of pores and throats.
        - Only one of pores, throats OR data are accepted

        """
        if self.__module__.split('.')[1] == 'Network': 
            net = self
            Ts = net.throats()
            Ps = net.pores()
            label = 'all'
        elif self.__module__.split('.')[1] in ['Phases','Algorithms']: 
            net = self._net
            Ts = net.throats()
            Ps = net.pores()
            label = 'all'          
        elif self.__module__.split('.')[1] in ['Geometry','Physics']:
            net = self._net
            Ts = net.throats(self.name)
            Ps = net.pores(self.name)
            label = self.name
        if sp.shape(data)[0] == self.Nt:
            #Upcast data to full network size
            temp = sp.ones((net.Nt,))*sp.nan
            temp[Ts] = data
            data = temp
            temp = sp.ones((net.Np,))*sp.nan
            for pore in Ps:
                neighborTs = net.find_neighbor_throats(pore)
                neighborTs = net.filter_by_label(throats=neighborTs,label=label)
                temp[pore] = sp.mean(data[neighborTs])
            values = temp[Ps]
        elif sp.shape(data)[0] == self.Np:
            #Upcast data to full network size
            temp = sp.ones((net.Np,))*sp.nan
            temp[Ps] = data
            data = temp
            Ps12 = net.find_connected_pores(throats=Ts,flatten=False)
            values = sp.mean(data[Ps12],axis=1)
        else:
            self._logger.error('Received data was an ambiguous length')
            raise Exception()
        return values
        
    def _interleave_data(self,prop,sources):
        r'''
        Retrieves requested property from associated objects, to produce a full 
        Np or Nt length array.
        
        Parameters
        ----------
        prop : string
            The property name to be retrieved
            
        sources : list
            List of object names OR objects from which data is retrieved
            
        Returns
        -------
        A full length (Np or Nt) array of requested property values.  
        
        Notes
        -----
        This makes an effort to maintain the data 'type' when possible; however
        when data is missing this can be tricky.  Float and boolean data is
        fine, but missing ints are converted to float when nans are inserted.
        '''
        element = prop.split('.')[0]
        #Utilize a pre-existing dummy 'temp' variable on the object to save time
        try:
            temp = self._temp[element]
        except:
            Np = sp.shape(self['pore.all'])[0]
            Nt = sp.shape(self['throat.all'])[0]
            self._temp = {}
            self._temp['pore'] = sp.empty((Np,))
            self._temp['throat'] = sp.empty((Nt,))
            temp = self._temp[element]
        dtypes = []
        prop_found = False  #Flag to indicate if prop was found on a sub-object
        for item in sources:
            #Check sources were given as list of objects OR names
            try: item.name
            except: item = self.find_object(obj_name=item)
            locations = self._get_indices(element=element,labels=item.name,mode='union')
            if prop not in item.keys():
                values = sp.ones_like(locations)*sp.nan
                dtypes.append('nan')
            else:
                prop_found = True
                values = item[prop]
                dtypes.append(values.dtype.name)
            temp[locations] = values  #Assign values
        #Check if requested prop was found on any sub-objects
        if prop_found == False:
            raise KeyError(prop)
        #Analyze and assign data type
        if sp.all([t in ['bool','nan'] for t in dtypes]):  # If all entries are 'bool' (or 'nan')
            temp = sp.array(temp,dtype='bool')*~sp.isnan(temp)
        elif sp.all([t == dtypes[0] for t in dtypes]) :  # If all entries are same type
            temp = sp.array(temp,dtype=dtypes[0])
        else:
            self._logger.warning('Retrieved data is of different type...converting to float')
        return temp
        
    def num_pores(self,labels='all',mode='union'):
        r'''
        Returns the number of pores of the specified labels

        Parameters
        ----------
        labels : list of strings, optional
            The pore labels that should be included in the count.  
            If not supplied, all pores are counted.
        labels : list of strings
            Label of pores to be returned
        mode : string, optional
            Specifies how the count should be performed.  The options are: 
            
            * 'union' : (default) All pores with ANY of the given labels are counted.

            * 'intersection' : Only pores with ALL the given labels are counted.
            
            * 'not_intersection' : Only pores with exactly one of the given labels are counted.
            
            * 'difference' : Only pores with none of the given labels are counted.
            
        Returns
        -------
        Np : int
            Number of pores with the specified labels 
            
        See Also
        --------
        `num_throats`, `count`
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.num_pores()
        125
        >>> pn.num_pores(labels=['top'])
        25
        >>> pn.num_pores(labels=['top','front'],mode='union') #'union' is default
        45
        >>> pn.num_pores(labels=['top','front'],mode='intersection')
        5
        >>> pn.num_pores(labels=['top','front'],mode='not_intersection')
        40
        
        '''
        if labels == 'all':
            Np = sp.shape(self['pore.all'])[0]
        else:
            #convert string to list, if necessary
            if type(labels) == str: 
                labels = [labels]
            #Count number of pores of specified type
            Ps = self.pores(labels=labels,mode=mode)
            Np = sp.shape(Ps)[0] 
        return Np
        
    @property
    def Np(self):
        return self.num_pores()  
            
    def num_throats(self,labels='all',mode='union'):
        r'''
        Return the number of throats of the specified labels

        Parameters
        ----------
        labels : list of strings, optional
            The throat labels that should be included in the count.  
            If not supplied, all throats are counted.
        mode : string, optional
            Specifies how the count should be performed.  The options are: 

            * 'union' : (default) All throats with ANY of the given labels are counted.

            * 'intersection' : Only throats with ALL the given labels are counted.

            * 'not_intersection' : Only throats with exactly one of the given labels are counted.
            
            * 'difference' : Only throats with none of the given labels are counted.

        Returns
        -------
        Nt : int
            Number of throats with the specified labels
            
        See Also
        --------
        `num_pores`, `count`

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.num_throats()
        300
        >>> pn.num_throats(labels=['top'])
        40
        >>> pn.num_throats(labels=['top','front'],mode='union') #'union' is default
        76
        >>> pn.num_throats(labels=['top','front'],mode='intersection')
        4
        >>> pn.num_throats(labels=['top','front'],mode='not_intersection')
        72
        
        '''
        if labels == 'all':
            Nt = sp.shape(self['throat.all'])[0]
        else:
            #convert string to list, if necessary
            if type(labels) == str: labels = [labels]
            #Count number of pores of specified type
            Ts = self.throats(labels=labels,mode=mode)
            Nt = sp.shape(Ts)[0]
        return Nt
        
    @property
    def Nt(self):
        return self.num_throats()    
        
    def _count(self,element=None):
        r'''
        Returns a dictionary containing the number of pores and throats in 
        the network, stored under the keys 'pore' or 'throat'
        
        Parameters
        ----------
        element : string, optional
            Can be either 'pore' , 'pores', 'throat' or 'throats', which 
            specifies which count to return.
            
        Returns
        -------
        A dictionary containing the number of pores and throats under the 
        'pore' and 'throat' key respectively.  
        
        See Also
        --------
        `num_pores`, `num_throats`
        
        Notes
        -----
        The ability to send plurals is useful for some types of 'programmatic'
        access.  For instance, the standard argument for locations is pores
        or throats.  If these are bundled up in a **kwargs dict then you can 
        just use the dict key in count() without removing the 's'.
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn._count()
        {'pore': 125, 'throat': 300}
        >>> pn._count('pore')
        125
        '''
        if element in ['pore','pores']:
            temp = self.num_pores()
        elif element in ['throat','throats']:
            temp = self.num_throats()
        elif element == None:
            temp = {}
            temp['pore'] = self.num_pores()
            temp['throat'] = self.num_throats()
        return temp
        
    def check_data_health(self,props=[],element='',quiet=False):
        r'''
        Check the health of pore and throat data arrays.  
        
        Parameters
        ----------
        element : string, optional
            Can be either 'pore' or 'throat', which will limit the checks to 
            only those data arrays.
        props : list of pore (or throat) properties, optional
            If given, will limit the health checks to only the specfied
            properties.  Also useful for checking existance.
        quiet : bool, optional
            By default this method will output a summary of the health check.
            This can be disabled by setting quiet to False.
            
        Returns
        -------
        Returns a True if all check pass, and False if any checks fail.  This
        is ideal for programatically checking data integrity prior to running
        an algorithm.
        
        '''
        health = {}
        flag = True
        if props == []:
            props = self.props(element)
        else:
            if type(props) == str:
                props = [props]
            if props[0].split('.')[0] not in ['pore','throat']:
                self._logger.error('Properties must be either pore or throat')
        for item in props:
            try: 
                if sp.sum(sp.isnan(self[item])) > 0:
                    health[item] = 'Has NaNs'
                    flag = False
                elif sp.shape(self[item])[0] == 1:
                    health[item] = 'Healthy Scalar'
                elif sp.shape(self[item])[0] == self._count(item.split('.')[0]):
                    health[item] = 'Healthy Vector'
                else:
                    health[item] = 'Wrong Length'
                    flag = False
            except: 
                health[item] = 'Does not exist'
                flag = False
        if quiet == False:
            pprint.pprint(health)
        return flag
        
    def __str__(self):
        header = '-'*60
        print(header)
        print(self.__module__.replace('__','')+': \t'+self.name)
        print(header)
        print("{a:<5s} {b:<35s} {c:<10s}".format(a='#', b='Properties', c='Valid Values'))
        print(header)
        count = 0
        props = self.props()
        props.sort()
        for item in props:
            count = count + 1
            prop=item
            if len(prop)>35:
                prop = prop[0:32]+'...'
            required = self._count(item.split('.')[0])
            defined = required - sp.sum(sp.isnan(self[item]))
            print("{a:<5d} {b:<35s} {c:>5d} / {d:<5d}".format(a=count, b=prop, c=defined, d=required))
        print(header)
        print("{a:<5s} {b:<35s} {c:<10s}".format(a='#', b='Labels', c='Assigned Locations'))
        print(header)
        count = 0
        labels = self.labels()
        labels.sort()
        for item in labels:
            count = count + 1
            prop=item
            if len(prop)>35:
                prop = prop[0:32]+'...'
            print("{a:<5d} {b:<35s} {c:<10d}".format(a=count, b=prop, c=sp.sum(self[item])))
        print(header)
        return ''

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)

