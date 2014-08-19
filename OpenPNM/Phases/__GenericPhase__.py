"""
module Physics
===============================================================================

"""
import sys, os
import OpenPNM
from OpenPNM.Base import Core
import scipy as sp

class GenericPhase(Core):
    r"""
    Base class to generate a generic phase object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common phases are included with OpenPNM and can be found under OpenPNM.Phases.

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this Phase should be attached
    components : list of OpenPNM Phase objects
        These Phase objects are ficticious or virtual phases that are the pure
        components from which the mixture is made.  These are used to calculate
        and store any pure component data.  If none are supplied then this 
        object will act like either a pure component, a mixture whose properties
        are well known (like air) and need not be found from consideration of
        the pure component properties.
    name : str, optional
        A unique string name to identify the Phase object, typically same as 
        instance name but can be anything.

    """
    def __init__(self,network,components=[],name=None,**kwargs):
        super(GenericPhase,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        # Attach objects to self for internal access
        self._net = network
        self.name = name
        [self._phases.append(comp) for comp in components]
        
        # Append this Phase to the Network
        network._phases.append(self)
        
        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = self._net['pore.all']
        self['throat.all'] = self._net['throat.all']        
        
        # Set default T and P since most propery models require it
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0
        
    def __setitem__(self,prop,value):
        for phys in self._physics:
            if prop in phys.keys():
                self._logger.error(prop+' is already defined in at least one associated Geometry object')
                return
        super(GenericPhase,self).__setitem__(prop,value)
        
    def __getitem__(self,key):
        if key not in self.keys():
            self._logger.debug(key+' not on Phase, constructing data from Physics')
            return self._interleave_data(key,sources=self._physics)
        else:
            return super(GenericPhase,self).__getitem__(key)
            
    def set_component(self,phase,mode='add'):
        r'''
        Pa
        '''
        if type(phase) == str:
            phase = self.find_object(obj_name=phase.name)
        if mode == 'add':
            if phase.name in self.phases():
                self._logger.error('Phase already present')
            else:
                self._phases.append(phase)
        elif mode == 'remove':
            if phase.name in self.phases():
                self._phases.remove(phase)
            else:
                self._logger.error('Phase not found')
                
    def regenerate(self,**kwargs):
        for item in self._phases:
            item.regenerate(**kwargs)
        super().regenerate(**kwargs)
    
    #Pull in doc string for the Core regenerate method
    regenerate.__doc__ = Core.regenerate.__doc__
            
    def check_mixture_health(self):
        r'''
        '''
        mole_sum = sp.zeros((self.Np,))
        for comp in self._phases:
            try:
                mole_sum = mole_sum + comp['pore.mole_fraction']
            except:
                pass
        return mole_sum
    


