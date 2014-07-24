import OpenPNM

def geom(pn):
    #==============================================================================
    '''Build Geometry'''
    #==============================================================================   
    #### GDL pores and throats for anode side
    ###----------------------------------------------------------------------------
    geom_GDL_anode = OpenPNM.Geometry.Toray090(name='GDL_anode',network=pn,loglevel=10)
    geom_GDL_anode.set_locations(pores=pn.pores('GDLa'),throats=pn.throats('GDLa'))
    
    #### GDL pores and throats for cathode side
    ###----------------------------------------------------------------------------
    geom_GDL_cathode = OpenPNM.Geometry.Toray090(name='GDL_cathode',network=pn,loglevel=10)
    geom_GDL_cathode.set_locations(pores=pn.pores('GDLc'),throats=pn.throats('GDLc'))
    
    #### CL pores and throats for anode side
    ###----------------------------------------------------------------------------
    geom_CL_anode = OpenPNM.Geometry.Toray090(name='CL_anode',network=pn,loglevel=10)
    geom_CL_anode.set_locations(pores=pn.pores('CLa'),throats=pn.throats('CLa'))
    
    #### CL pores and throats for cathode side
    ###----------------------------------------------------------------------------
    geom_CL_cathode = OpenPNM.Geometry.Toray090(name='CL_cathode',network=pn,loglevel=10)
    geom_CL_cathode.set_locations(pores=pn.pores('CLc'),throats=pn.throats('CLc'))
    
    #### ME pores and throats
    ###----------------------------------------------------------------------------
    geom_PEM = OpenPNM.Geometry.Toray090(name='PEM',network=pn,loglevel=10)
    geom_PEM.set_locations(pores=pn.pores('PEM'),throats=pn.throats('PEM'))
    
    
    T_GDL_CL_anode = pn.find_interface_throats(['GDLa','CLa'])
    T_GDL_CL_cathode = pn.find_interface_throats(['GDLc','CLc'])
    T_PEM_CL_anode = pn.find_interface_throats(['PEM','CLa'])
    T_PEM_CL_cathode = pn.find_interface_throats(['PEM','CLc'])
    
    #### GDL and CL, interfacial throats, anode side 
    ###----------------------------------------------------------------------------
    geom_GDL_CL_anode = OpenPNM.Geometry.Toray090(name='GDL_CL_anode_interface',network=pn,loglevel=10)
    geom_GDL_CL_anode.set_locations(throats=T_GDL_CL_anode)
    
    #### GDL and CL, interfacial throats, cathode side 
    ###----------------------------------------------------------------------------
    geom_GDL_CL_cathode = OpenPNM.Geometry.Toray090(name='GDL_CL_cathode_interface',network=pn,loglevel=10)
    geom_GDL_CL_cathode.set_locations(throats=T_GDL_CL_cathode)
    
    #### ME and CL, interfacial throats, anode side 
    ###----------------------------------------------------------------------------
    geom_PEM_CL_anode = OpenPNM.Geometry.Toray090(name='PEM_CL_anode_interface',network=pn,loglevel=10)
    geom_PEM_CL_anode.set_locations(throats= T_PEM_CL_anode)
    
    #### ME and CL, interfacial throats, cathode side 
    ###----------------------------------------------------------------------------
    geom_PEM_CL_cathode = OpenPNM.Geometry.Toray090(name='PEM_CL_cathode_interface',network=pn,loglevel=10)
    geom_PEM_CL_cathode.set_locations(throats= T_PEM_CL_cathode)
    
    ###----------------------------------------------------------------------------
    pn.regenerate_geometries()
    
    
    
def fluid(pn):
    #==============================================================================
    '''Build Fluids'''
    #==============================================================================
    air = OpenPNM.Fluids.Air(name='Air',network=pn, loglevel=50)
    air.apply_conditions(temperature=350, pressure=200000)
    air.add_property(prop='electrical_conductivity',model='constant',value=5e-12)
    ###----------------------------------------------------------------------------
    solid = OpenPNM.Fluids.GenericFluid(name='Solid',network=pn, loglevel=10)
    solid.apply_conditions(temperature=350)
    solid.add_property(prop='electrical_conductivity',model='constant',value=2e5)
    ###----------------------------------------------------------------------------
    Nafion = OpenPNM.Fluids.GenericFluid(name='Nafion',network=pn, loglevel=50)
    Nafion.add_property(prop='diffusivity',model='constant',value=2e-10)
    Nafion.add_property(prop='molar_density',model='constant',value=44445)
    Nafion.add_property(prop='electrical_conductivity',model='constant',value=1e3)
    ###----------------------------------------------------------------------------
    hydrogen = OpenPNM.Fluids.GenericFluid(name='Hydrogen',network=pn, loglevel=50)
    hydrogen.set_pore_data(prop='Tc',data=33.49)
    hydrogen.set_pore_data(prop='Pc',data=1.296e6)
    hydrogen.set_pore_data(prop='MW',data=0.002)
    hydrogen.add_property(prop='diffusivity',model='constant',value=6.1e-5)
    hydrogen.add_property(prop='viscosity',model='constant',value=8.3e-5)
    hydrogen.add_property(prop='molar_density',model='ideal_gas',R=8.314)
    hydrogen.add_property(prop='electrical_conductivity',model='constant',value=5e-12)    
    
    #Use Network's Fluid regeneration method
    pn.regenerate_fluids()
    pn._air = air
    pn._hydrogen = hydrogen
    pn._Nafion = Nafion
    pn._solid = solid
    
    return(air,solid,Nafion,hydrogen)
    
def physics(pn,PEM_Lc,PEM_sigma,C_ionomer_fraction,CL_Lc,porosity_CL,tort_CL):    
        
    #==============================================================================
    '''Build Physics Objects'''
    #==============================================================================
    air = pn._air
    hydrogen = pn._hydrogen
    solid = pn._solid
    Nafion = pn._Nafion
    #### Physics for Nafion in CL_anode
    sigma_val_CL = PEM_sigma*C_ionomer_fraction**1.5*(PEM_Lc**2/PEM_Lc)
    ###----------------------------------------------------------------------------
    phys_Nafion_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_CL_anode',network=pn, fluid=Nafion,geometry=pn.find_object(obj_name='CL_anode'),loglevel=50)
    phys_Nafion_CL_anode.add_property(prop='electronic_conductance',  model='constant',value=sigma_val_CL)
    
    #### Physics for Nafion in CL_cathode
    ###----------------------------------------------------------------------------
    phys_Nafion_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_CL_cathode',network=pn, fluid=Nafion,geometry=pn.find_object(obj_name='CL_cathode'),loglevel=50)
    phys_Nafion_CL_cathode.add_property(prop='electronic_conductance',  model='constant',value=sigma_val_CL)
    
    
    #### Physics for Nafion in ME
    sigma_val_PEM = PEM_sigma*(PEM_Lc**2/PEM_Lc)
    ###----------------------------------------------------------------------------
    phys_Nafion_PEM = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_ME',network=pn, fluid=Nafion,geometry=pn.find_object(obj_name='PEM'),loglevel=50)
    phys_Nafion_PEM.add_property(prop='electronic_conductance',  model='constant',value=sigma_val_PEM)
    
        
    #### Physics for Nafion in GDL_CL_cathode interface
    ###----------------------------------------------------------------------------
    phys_Nafion_GDL_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_GDL_CL_cathode',network=pn, fluid=Nafion,geometry=pn.find_object(obj_name='GDL_CL_cathode_interface'),loglevel=50)
    phys_Nafion_GDL_CL_cathode.add_property(prop='electronic_conductance', model='constant',value=1e-100)

    #### Physics for Nafion in GDL_CL_anode interface
    ###----------------------------------------------------------------------------
    phys_Nafion_GDL_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_GDL_CL_anode',network=pn, fluid=Nafion,geometry=pn.find_object(obj_name='GDL_CL_anode_interface'),loglevel=50)
    phys_Nafion_GDL_CL_anode.add_property(prop='electronic_conductance', model='constant',value=1e-100)


    #### Physics for Nafion in PEM_CL_anode interface
    ###----------------------------------------------------------------------------
    phys_Nafion_PEM_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_PEM_CL_anode',network=pn, fluid=Nafion,geometry=pn.find_object(obj_name='PEM_CL_anode_interface'),loglevel=50)
    phys_Nafion_PEM_CL_anode.add_property(prop='electronic_conductance', model='constant',value=sigma_val_PEM)
    
    #### Physics for Nafion in PEM_CL_cathode interface
    ###----------------------------------------------------------------------------
    phys_Nafion_PEM_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_PEM_CL_cathode',network=pn, fluid=Nafion,geometry=pn.find_object(obj_name='PEM_CL_cathode_interface'),loglevel=50)
    phys_Nafion_PEM_CL_cathode.add_property(prop='electronic_conductance',model='constant',value=sigma_val_PEM)    
    
    #==============================================================================================================================
    #==============================================================================================================================
    
    #### Physics for hydrogen in GDL_anode
    ###----------------------------------------------------------------------------
    phys_hydrogen_GDL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_GDL_anode',network=pn, fluid=hydrogen,geometry=pn.find_object(obj_name='GDL_anode'),loglevel=50)
    phys_hydrogen_GDL_anode.add_property(prop='diffusive_conductance', model='bulk_diffusion')
    
    #### Physics for hydrogen in CL_anode
    ###----------------------------------------------------------------------------
    phys_hydrogen_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_CL_anode',network=pn, fluid=hydrogen,geometry=pn.find_object(obj_name='CL_anode'),loglevel=50)
    DH_val_CL = hydrogen['pore.molar_density']*hydrogen['pore.diffusivity']*(porosity_CL/tort_CL)*(CL_Lc**2/CL_Lc)

    phys_hydrogen_CL_anode.add_property(prop='diffusive_conductance', model= 'constant', value = DH_val_CL)  
    
    
    #### Physics for hydrogen in GDL_CL_anode interface
    ###----------------------------------------------------------------------------
    phys_hydrogen_GDL_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_GDL_CL_anode',network=pn, fluid=hydrogen,geometry=pn.find_object(obj_name='GDL_CL_anode_interface'),loglevel=50)
    phys_hydrogen_GDL_CL_anode.add_property(prop='diffusive_conductance', model= 'constant', value = DH_val_CL)

    #### Physics for hydrogen in GDL_CL_anode interface
    ###----------------------------------------------------------------------------
    phys_hydrogen_PEM_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_PEM_CL_anode',network=pn, fluid=hydrogen,geometry=pn.find_object(obj_name='PEM_CL_anode_interface'),loglevel=50)
    phys_hydrogen_PEM_CL_anode.add_property(prop='diffusive_conductance', model='constant',value=1e-100)    
   
    #==============================================================================================================================
    #==============================================================================================================================   
   
    #### Physics for air in GDL_cathode
    ###----------------------------------------------------------------------------
    phys_air_GDL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_GDL_cathode',network=pn, fluid=air,geometry=pn.find_object(obj_name='GDL_cathode'),loglevel=50)
    phys_air_GDL_cathode.add_property(prop='diffusive_conductance', model='bulk_diffusion')
    
    
    #### Physics for air in CL_cathode
    ###----------------------------------------------------------------------------
    phys_air_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_CL_cathode',network=pn, fluid=air,geometry=pn.find_object(obj_name='CL_cathode'),loglevel=50)
    DO_val_CL = air['pore.molar_density']*air['pore.diffusivity']*(porosity_CL/tort_CL)*(CL_Lc**2/CL_Lc)
    phys_air_CL_cathode.add_property(prop='diffusive_conductance', model='constant', value = DO_val_CL)
    
    #### Physics for air in GDL_CL_cathode interface
    ###----------------------------------------------------------------------------
    phys_air_GDL_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_GDL_CL_cathode',network=pn, fluid=air,geometry=pn.find_object(obj_name='GDL_CL_cathode_interface'),loglevel=50)
    phys_air_GDL_CL_cathode.add_property(prop='diffusive_conductance', model='constant', value = DO_val_CL)

    #### Physics for air in GDL_CL_cathode interface
    ###----------------------------------------------------------------------------
    phys_air_PEM_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_PEM_CL_cathode',network=pn, fluid=air,geometry=pn.find_object(obj_name='PEM_CL_cathode_interface'),loglevel=50)
    phys_air_PEM_CL_cathode.add_property(prop='diffusive_conductance', model='constant',value=1e-100)    
     
    #==============================================================================================================================
    #==============================================================================================================================      
     
    #### Physics for solid in GDL_anode
    ###----------------------------------------------------------------------------
    phys_solid_GDL_anode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_anode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='GDL_anode'),loglevel=50)
    phys_solid_GDL_anode.add_property(prop='electronic_conductance', model='series_resistors')
    
    #### Physics for solid in GDL_cathode
    ###----------------------------------------------------------------------------
    phys_solid_GDL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_cathode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='GDL_cathode'),loglevel=50)
    phys_solid_GDL_cathode.add_property(prop='electronic_conductance', model='series_resistors')
    
    #### Physics for solid in CL_anode
    ###----------------------------------------------------------------------------
    phys_solid_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_solid_CL_anode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='CL_anode'),loglevel=50)
    phys_solid_CL_anode.add_property(prop='electronic_conductance', model='series_resistors')
    
    #### Physics for solid in CL_cathode
    ###----------------------------------------------------------------------------
    phys_solid_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_solid_CL_cathode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='CL_cathode'),loglevel=50)
    phys_solid_CL_cathode.add_property(prop='electronic_conductance', model='series_resistors')
    

    #### Physics for solid in GDL_CL_anode interface
    ###----------------------------------------------------------------------------
    phys_solid_GDL_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_CL_anode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='GDL_CL_anode_interface'),loglevel=50)
    phys_solid_GDL_CL_anode.add_property(prop='electronic_conductance', model='series_resistors')

    #### Physics for solid in GDL_CL_cathode interface
    ###----------------------------------------------------------------------------
    phys_solid_GDL_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_CL_cathode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='GDL_CL_cathode_interface'),loglevel=50)
    phys_solid_GDL_CL_cathode.add_property(prop='electronic_conductance', model='series_resistors')


    #### Physics for solid in GDL_CL_anode interface
    ###----------------------------------------------------------------------------
    phys_solid_PEM_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_solid_PEM_CL_anode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='PEM_CL_anode_interface'),loglevel=50)
    phys_solid_PEM_CL_anode.add_property(prop='electronic_conductance', model='series_resistors')

    #### Physics for solid in GDL_CL_cathode interface
    ###----------------------------------------------------------------------------
    phys_solid_PEM_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_solid_PEM_CL_cathode',network=pn, fluid=solid,geometry=pn.find_object(obj_name='PEM_CL_cathode_interface'),loglevel=50)
    phys_solid_PEM_CL_cathode.add_property(prop='electronic_conductance', model='series_resistors')

    #==============================================================================================================================
    #==============================================================================================================================
    #Use Network's Physics regeneration method
    pn.regenerate_physics()
    
