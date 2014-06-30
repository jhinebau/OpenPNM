import OpenPNM
import scipy as sp
import OpenPNM.Utilities.misc as misc


#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic2(loglevel=10,name='net')

domains = {}
domains['GDLa'] = {'region': [[0,0,0],[0.0005,0.0005,0.0003]], 'Lc': 0.00005}
domains['CLa']  = {'region': [[0,0,0.0003],[0.0005,0.0005,0.00035]], 'Lc': 0.00001}
domains['PEM']  = {'region': [[0,0,0.00035],[0.0005,0.0005,0.0004]], 'Lc': 0.00001}
domains['CLc']  = {'region': [[0,0,0.0004],[0.0005,0.0005,0.00045]], 'Lc': 0.00001}
domains['GDLc'] = {'region': [[0,0,0.00045],[0.0005,0.0005,0.00075]], 'Lc': 0.00005}

pn.generate(domains=domains)

pn.save_object('big_net.npz')

###==============================================================================
##'''Build Geometry'''
###==============================================================================

X = pn.get_data(prop='coords',pores='all')[:,0]
Y = pn.get_data(prop='coords',pores='all')[:,1]
Z = pn.get_data(prop='coords',pores='all')[:,2]

P_GDL_a = pn.tomask(pores=pn.pores('GDLa'))
P_GDL_c = pn.tomask(pores=pn.pores('GDLc'))
P_CL_a = pn.tomask(pores=pn.pores('CLa'))
P_CL_c = pn.tomask(pores=pn.pores('CLc'))
P_PEM = pn.tomask(pores=pn.pores('PEM'))

T_GDL_a = pn.throats('GDLa')
T_GDL_c = pn.throats('GDLc')
T_CL_a = pn.throats('CLa')
T_CL_c = pn.throats('CLc')
T_PEM = pn.throats('PEM')
T_GDL_CL_a = pn.find_interface_throats(['GDLa','CLa'])
T_GDL_CL_c = pn.find_interface_throats(['GDLc','CLc'])
T_PEM_CL_a = pn.find_interface_throats(['PEM','CLa'])
T_PEM_CL_c = pn.find_interface_throats(['PEM','CLc'])



#==============================================================================
'''Build Geometry'''
#==============================================================================

#### GDL pores and throats for anode side
###----------------------------------------------------------------------------
geom_GDL_anode = OpenPNM.Geometry.Toray090(name='GDL_anode',network=pn,loglevel=10)
T_GDL_a_not_interface = pn.find_neighbor_throats(P_GDL_a,flatten=True,mode='intersection')
geom_GDL_anode.set_locations(pores=P_GDL_a, throats=T_GDL_a_not_interface)

#### GDL pores and throats for cathode side
###----------------------------------------------------------------------------
geom_GDL_cathode = OpenPNM.Geometry.Toray090(name='GDL_cathode',network=pn,loglevel=10)
T_GDL_c_not_interface = pn.find_neighbor_throats(P_GDL_c,flatten=True,mode='intersection')
geom_GDL_cathode.set_locations(pores=P_GDL_c, throats=T_GDL_c_not_interface)

#### CL pores and throats for anode side
###----------------------------------------------------------------------------
geom_CL_anode = OpenPNM.Geometry.Toray090(name='CL_anode',network=pn,loglevel=10)
T_CL_a_not_interface = pn.find_neighbor_throats(P_CL_a,flatten=True,mode='intersection')
geom_CL_anode.set_locations(pores=P_CL_a,throats=T_CL_a_not_interface)

#### CL pores and throats for cathode side
###----------------------------------------------------------------------------
geom_CL_cathode = OpenPNM.Geometry.Toray090(name='CL_cathode',network=pn,loglevel=10)
T_CL_c_not_interface = pn.find_neighbor_throats(P_CL_c,flatten=True,mode='intersection')
geom_CL_cathode.set_locations(pores=P_CL_c,throats=T_CL_c_not_interface)

#### ME pores and throats
###----------------------------------------------------------------------------
geom_PEM = OpenPNM.Geometry.Toray090(name='PEM',network=pn,loglevel=10)
T_PEM_not_interface = pn.find_neighbor_throats(P_PEM,flatten=True,mode='intersection')
geom_PEM.set_locations(pores=P_PEM,throats=T_PEM_not_interface)


#### GDL and CL, interfacial throats, anode side 
###----------------------------------------------------------------------------
geom_GDL_CL_anode = OpenPNM.Geometry.Toray090(name='GDL_CL_anode_interface',network=pn,loglevel=10)
geom_GDL_CL_anode.set_locations(throats=T_GDL_CL_a)

#### GDL and CL, interfacial throats, cathode side 
###----------------------------------------------------------------------------
geom_GDL_CL_cathode = OpenPNM.Geometry.Toray090(name='GDL_CL_cathode_interface',network=pn,loglevel=10)
geom_GDL_CL_cathode.set_locations(throats=T_GDL_CL_c)

#### ME and CL, interfacial throats, anode side 
###----------------------------------------------------------------------------
geom_PEM_CL_anode = OpenPNM.Geometry.Toray090(name='PEM_CL_anode_interface',network=pn,loglevel=10)
geom_PEM_CL_anode.set_locations(throats= T_PEM_CL_a)

#### ME and CL, interfacial throats, cathode side 
###----------------------------------------------------------------------------
geom_PEM_CL_cathode = OpenPNM.Geometry.Toray090(name='PEM_CL_cathode_interface',network=pn,loglevel=10)
geom_PEM_CL_cathode.set_locations(throats= T_PEM_CL_c)

###----------------------------------------------------------------------------
pn.regenerate_geometries()

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
hydrogen.set_pore_data(prop='Tc',data=132.65)
hydrogen.set_pore_data(prop='Pc',data=3.771e6)
hydrogen.set_pore_data(prop='MW',data=0.002)
hydrogen.add_property(prop='diffusivity',model='constant',value=2e-5)
hydrogen.add_property(prop='viscosity',model='constant',value=1.9e-5)
hydrogen.add_property(prop='molar_density',model='ideal_gas',R=8.314)
hydrogen.add_property(prop='electrical_conductivity',model='constant',value=5e-12)

#Use Network's Fluid regeneration method
pn.regenerate_fluids()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================

#### Physics for Nafion in CL_anode
###----------------------------------------------------------------------------
phys_Nafion_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_CL_anode',network=pn, fluid=Nafion,geometry=geom_CL_anode,loglevel=50)
phys_Nafion_CL_anode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for Nafion in CL_cathode
###----------------------------------------------------------------------------
phys_Nafion_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_CL_cathode',network=pn, fluid=Nafion,geometry=geom_CL_cathode,loglevel=50)
phys_Nafion_CL_cathode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for Nafion in GDL_CL_cathode interface
###----------------------------------------------------------------------------
phys_Nafion_GDL_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_GDL_CL_cathode',network=pn, fluid=Nafion,geometry=geom_GDL_CL_cathode,loglevel=50)
phys_Nafion_GDL_CL_cathode.add_property(prop='electronic_conductance', model='constant',value=1e-100)

#### Physics for Nafion in GDL_CL_anode interface
###----------------------------------------------------------------------------
phys_Nafion_GDL_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_GDL_CL_anode',network=pn, fluid=Nafion,geometry=geom_GDL_CL_anode,loglevel=50)
phys_Nafion_GDL_CL_anode.add_property(prop='electronic_conductance', model='constant',value=1e-100)

#### Physics for Nafion in PEM
###----------------------------------------------------------------------------
phys_Nafion_PEM = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_ME',network=pn, fluid=Nafion,geometry=geom_PEM,loglevel=50)
phys_Nafion_PEM.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for Nafion in PEM_CL_anode interface
###----------------------------------------------------------------------------
phys_Nafion_PEM_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_PEM_CL_anode',network=pn, fluid=Nafion,geometry=geom_PEM_CL_anode,loglevel=50)
phys_Nafion_PEM_CL_anode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for Nafion in PEM_CL_cathode interface
###----------------------------------------------------------------------------
phys_Nafion_PEM_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_PEM_CL_cathode',network=pn, fluid=Nafion,geometry=geom_PEM_CL_cathode,loglevel=50)
phys_Nafion_PEM_CL_cathode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for hydrogen in GDL_anode
###----------------------------------------------------------------------------
phys_hydrogen_GDL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_GDL_anode',network=pn, fluid=hydrogen,geometry=geom_GDL_anode,loglevel=50)
phys_hydrogen_GDL_anode.add_property(prop='diffusive_conductance', model='bulk_diffusion')

#### Physics for hydrogen in CL_anode
###----------------------------------------------------------------------------
phys_hydrogen_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_CL_anode',network=pn, fluid=hydrogen,geometry=geom_CL_anode,loglevel=50)
phys_hydrogen_CL_anode.add_property(prop='diffusive_conductance', model='bulk_diffusion')

#### Physics for hydrogen in GDL_CL_anode interface
###----------------------------------------------------------------------------
phys_hydrogen_GDL_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_GDL_CL_anode',network=pn, fluid=hydrogen,geometry=geom_GDL_CL_anode,loglevel=50)
phys_hydrogen_GDL_CL_anode.add_property(prop='diffusive_conductance', model='bulk_diffusion')

#### Physics for hydrogen in GDL_CL_anode interface
###----------------------------------------------------------------------------
phys_hydrogen_PEM_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_PEM_CL_anode',network=pn, fluid=hydrogen,geometry=geom_PEM_CL_anode,loglevel=50)
phys_hydrogen_PEM_CL_anode.add_property(prop='diffusive_conductance', model='constant',value=1e-100)

#### Physics for air in GDL_cathode
###----------------------------------------------------------------------------
phys_air_GDL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_GDL_cathode',network=pn, fluid=air,geometry=geom_GDL_cathode,loglevel=50)
phys_air_GDL_cathode.add_property(prop='diffusive_conductance', model='bulk_diffusion')

#### Physics for air in CL_cathode
###----------------------------------------------------------------------------
phys_air_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_CL_cathode',network=pn, fluid=air,geometry=geom_CL_cathode,loglevel=50)
phys_air_CL_cathode.add_property(prop='diffusive_conductance', model='bulk_diffusion')

#### Physics for air in GDL_CL_cathode interface
###----------------------------------------------------------------------------
phys_air_GDL_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_GDL_CL_cathode',network=pn, fluid=air,geometry=geom_GDL_CL_cathode,loglevel=50)
phys_air_GDL_CL_cathode.add_property(prop='diffusive_conductance', model='bulk_diffusion')

#### Physics for air in GDL_CL_cathode interface
###----------------------------------------------------------------------------
phys_air_PEM_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_PEM_CL_cathode',network=pn, fluid=air,geometry=geom_PEM_CL_cathode,loglevel=50)
phys_air_PEM_CL_cathode.add_property(prop='diffusive_conductance', model='constant',value=1e-100)
#### Physics for solid in GDL_anode
###----------------------------------------------------------------------------
phys_solid_GDL_anode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_anode',network=pn, fluid=solid,geometry=geom_GDL_anode,loglevel=50)
phys_solid_GDL_anode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for solid in GDL_cathode
###----------------------------------------------------------------------------
phys_solid_GDL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_cathode',network=pn, fluid=solid,geometry=geom_GDL_cathode,loglevel=50)
phys_solid_GDL_cathode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for solid in CL_anode
###----------------------------------------------------------------------------
phys_solid_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_solid_CL_anode',network=pn, fluid=solid,geometry=geom_CL_anode,loglevel=50)
phys_solid_CL_anode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for solid in CL_cathode
###----------------------------------------------------------------------------
phys_solid_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_solid_CL_cathode',network=pn, fluid=solid,geometry=geom_CL_cathode,loglevel=50)
phys_solid_CL_cathode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for solid in GDL_CL_anode interface
###----------------------------------------------------------------------------
phys_solid_GDL_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_CL_anode',network=pn, fluid=solid,geometry=geom_GDL_CL_anode,loglevel=50)
phys_solid_GDL_CL_anode.add_property(prop='electronic_conductance', model='series_resistors')

#### Physics for solid in GDL_CL_cathode interface
###----------------------------------------------------------------------------
phys_solid_GDL_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_solid_GDL_CL_cathode',network=pn, fluid=solid,geometry=geom_GDL_CL_cathode,loglevel=50)
phys_solid_GDL_CL_cathode.add_property(prop='electronic_conductance', model='series_resistors')

#Use Network's Physics regeneration method
pn.regenerate_physics()


solid.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
Nafion.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
hydrogen.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
air.set_data(prop='occupancy',throats=pn.throats('all'),data=1)



#==============================================================================
'''Begin Simulations'''
#==============================================================================

'''PEM effective resistance '''
##=============================================================================
PEM_Reff_alg = OpenPNM.Algorithms.OhmicConduction(name='PEM_Reff_alg',loglevel=10, network=pn)
PEM_Reff_alg_BC1_pores =  pn.pores('CLa')[Z[P_CL_a]==min(Z[P_CL_a])]
PEM_Reff_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue = 0.8, pores = PEM_Reff_alg_BC1_pores)
#------------------------------------------------------------------------------
PEM_Reff_alg_BC2_pores = pn.pores('CLc')[Z[P_CL_c]==max(Z[P_CL_c])] 
PEM_Reff_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue = 0.4, pores = PEM_Reff_alg_BC2_pores)
#------------------------------------------------------------------------------
Nafion.set_nan_value(element='throat',prop='electronic_conductance',value=1e-40)
#------------------------------------------------------------------------------
PEM_Reff_alg.setup(fluid=Nafion)
PEM_Reff_alg.run()
PEM_Reff_alg.update()
#------------------------------------------------------------------------------


Ps = PEM_Reff_alg.pores(labels='pore.Dirichlet')
BCs = sp.unique(PEM_Reff_alg['pore.bcval_Dirichlet'][Ps])
#Find flow through inlet face
Pn = pn.find_neighbor_pores(pores=PEM_Reff_alg_BC1_pores,excl_self=True)
Pn = Pn[-sp.in1d(Pn,pn.pores('GDLa'))]
Ts = pn.find_connecting_throat(PEM_Reff_alg_BC1_pores,Pn)
g = Nafion['throat.electronic_conductance'][Ts]
s = Nafion['throat.occupancy'][Ts]
xin = PEM_Reff_alg['pore.voltage'][PEM_Reff_alg_BC1_pores]
xout = PEM_Reff_alg['pore.voltage'][Pn]
i = g*s*(xin - xout)
Reff = sp.absolute((BCs[0]-BCs[1])/sp.sum(i))
#------------------------------------------------------------------------------

R = 8.314
F = 9.64870e4
C_ref = 40
alpha_anode = 0.5
alpha_cathode = 0.5
T = 353
I0 = 1e4
active_area = 1e-2
I = 0.01
S_H2 = I/(2*F)
S_O2 = I/(4*F)

#==============================================================================
'''Necessary transport algorithms and the fixed boundary conditions'''

##=============================================================================
'''Hydrogen Fickian Diffusion'''
##=============================================================================
Fick_hydrogen = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_hydrogen',loglevel=50, network=pn)
#------------------------------------------------------------------------------
Fick_hydrogen_BC1_pores = pn.pores('GDLa')[Z[P_GDL_a]==min(Z[P_GDL_a])]
Fick_hydrogen.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.8, pores=Fick_hydrogen_BC1_pores)
#------------------------------------------------------------------------------
#Fick_hydrogen_BC2_pores = P_CL_a
Fick_hydrogen_BC2_pores = pn.pores('CLa')[Z[P_CL_a]==max(Z[P_CL_a])]

Fick_hydrogen.set_boundary_conditions(bctype='Dirichlet', bcvalue= 0.4, pores=Fick_hydrogen_BC2_pores)

#    Fick_hydrogen.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue= -S_H2, pores=Fick_hydrogen_BC2_pores,mode='overwrite')

#------------------------------------------------------------------------------
#Fick_hydrogen_BC3_pores = ~(P_GDL_a+P_CL_a)
#Fick_hydrogen.set_boundary_conditions(bctype='Neumann_group', bcvalue = 1e-60, pores=Fick_hydrogen_BC3_pores)
#------------------------------------------------------------------------------
hydrogen.set_nan_value(element='throat',prop='diffusive_conductance',value=1e-40)
#
Fick_hydrogen.setup(fluid=hydrogen)
Fick_hydrogen.run()
Fick_hydrogen.update()
#
##=============================================================================
'''Air Fickian Diffusion'''
##=============================================================================
Fick_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_air',loglevel=50, network=pn)
#------------------------------------------------------------------------------
Fick_air_BC1_pores = pn.pores('GDLc')[Z[P_GDL_c]==max(Z[P_GDL_c])]
Fick_air.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.8, pores=Fick_air_BC1_pores)
#------------------------------------------------------------------------------
#Fick_air_BC2_pores = P_CL_c
Fick_air_BC2_pores = pn.pores('CLc')[Z[P_CL_c]==min(Z[P_CL_c])]

Fick_air.set_boundary_conditions(bctype='Dirichlet', bcvalue= 0.4, pores=Fick_air_BC2_pores)

#    Fick_air.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue= -S_H2, pores=Fick_hydrogen_BC2_pores,mode='overwrite')
#------------------------------------------------------------------------------
#Fick_air_BC3_pores = ~(P_GDL_cathode+P_CL_cathode)
#Fick_air.set_boundary_conditions(bctype='Neumann_insulated', pores=Fick_air_BC3_pores)
#------------------------------------------------------------------------------
air.set_nan_value(element='throat',prop='diffusive_conductance',value=1e-40)
Fick_air.setup(fluid=air)
Fick_air.run()
Fick_air.update()
#

##=============================================================================
'''Local mass consumption anode'''
##=============================================================================
mass_anode = Fick_hydrogen.rate(pn.pores('CLa'),mode='single')
##==============================================================================
##=============================================================================
'''Local mass consumption cathode'''
##=============================================================================
mass_cathode = Fick_air.rate(pn.pores('CLc'),mode='single')
##==============================================================================

##=============================================================================
'''Local current anode'''
##=============================================================================

I_local_anode = 2*F*mass_anode
##=============================================================================
'''Local current cathode'''
##=============================================================================

I_local_cathode = 4*F*mass_cathode

##=============================================================================
'''etta anode'''
##=============================================================================
C_H2 = hydrogen.get_data(prop='mole_fraction',pores=P_CL_a)*hydrogen.get_data(prop='molar_density',pores=P_CL_a)
etta_a = R*T/(alpha_anode)*sp.log(I_local_anode/(I0*active_area)*((C_ref/C_H2)**1.5))

##=============================================================================
'''etta cathode'''
##=============================================================================
C_O2 = air.get_data(prop='mole_fraction',pores=P_CL_c)*air.get_data(prop='molar_density',pores=P_CL_c)
etta_c = -R*T/(alpha_cathode)*sp.log(I_local_cathode/(I0*active_area)*((C_ref/C_O2)**1.5))

Vcell = 1.233-sp.average(etta_a)-sp.average(etta_c)-I*Reff

#Ohm_solid = OpenPNM.Algorithms.OhmicConduction(name='Ohmic_alg_solid',loglevel=10, network=pn)
#Ohm_solid_BC1_pores =  pn.pores('GDL_anode')[X[P_GDL_anode]==min(X[P_GDL_anode])]
#Ohm_solid.set_boundary_conditions(bctype='Dirichlet', bcvalue=1e-10, pores=Ohm_solid_BC1_pores)
##------------------------------------------------------------------------------
#Ohm_solid_BC2_pores = P_ME
#Ohm_solid.set_boundary_conditions(bctype='Neumann_insulated', pores=Ohm_solid_BC2_pores)
##------------------------------------------------------------------------------
#Ohm_solid_BC3_pores =  pn.pores('CL_anode')
##Ohm_solid_BC3_pores = pn.pores('CL_anode')[X[P_CL_anode]==max(X[P_CL_anode])]
##------------------------------------------------------------------------------
#Ohm_solid_BC4_pores = pn.pores('CL_cathode')
##Ohm_solid_BC4_pores = pn.pores('CL_cathode')[X[P_CL_cathode]==min(X[P_CL_cathode])]
#
##------------------------------------------------------------------------------
#Ohm_solid_BC5_pores =  pn.pores('GDL_cathode')[X[P_GDL_cathode]==max(X[P_GDL_cathode])]
##------------------------------------------------------------------------------
#solid.set_nan_value(element='throat',prop='electronic_conductance',value=1e-40)
#
##==============================================================================
#'''Cell Model'''
##==============================================================================
#voltage = sp.array([0.6])
#current =[]
#
#for V_cell in voltage:
#    print('==================================')
#    print('V = ',V_cell)
#    '''Step1- Guess current'''
#    guess_I = 1e-3
#    error = 10
#    counter = 0
#    Ohm_solid.set_boundary_conditions(bctype='Dirichlet', pores=Ohm_solid_BC5_pores,mode='remove')
#    Ohm_solid.set_boundary_conditions(bctype='Dirichlet', bcvalue = V_cell, pores=Ohm_solid_BC5_pores)
##    while   error>0.1:
#    #------------------------------------------------------------------------------
#    Ohm_solid.set_boundary_conditions(bctype='Neumann_rate_group', pores = Ohm_solid_BC3_pores,mode='remove')
#    Ohm_solid.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue = -guess_I, pores=Ohm_solid_BC3_pores)
#    
#    #------------------------------------------------------------------------------
#    Ohm_solid.set_boundary_conditions(bctype='Neumann_rate_group', pores = Ohm_solid_BC4_pores,mode='remove')
#    Ohm_solid.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue = guess_I, pores=Ohm_solid_BC4_pores)
#    #------------------------------------------------------------------------------
#    # Run simulation
#    Ohm_solid.run(active_fluid=solid)
#    Ohm_solid.update()
#
#    ##=============================================================================
#    '''Step2_ Obtain the mass source term based on the current'''
#    
#    S_H2 = guess_I/(2*F)
#    S_O2 = guess_I/(4*F)
#    ##=============================================================================
#    '''Step3_ Mass transfer algorithm for cathode and anode'''        
#    Fick_hydrogen.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue= -S_H2, pores=Fick_hydrogen_BC2_pores,mode='overwrite')
#    ##-----------------------------------------------------------------------------
#    ## Run Hydrogen Fickian Diffusion
#    Fick_hydrogen.run(active_fluid=hydrogen)
#    Fick_hydrogen.update()
#    #==============================================================================
#    Fick_air.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue= -S_O2, pores=Fick_air_BC2_pores,mode='overwrite')
#    ##-----------------------------------------------------------------------------
#    ## Run air Fickian Diffusion
#    Fick_air.run(active_fluid=air)
#    Fick_air.update()
#    #==============================================================================
#    Ohm_Nafion.set_boundary_conditions(bctype='Neumann_rate_group', pores = Ohm_Nafion_BC4_pores, mode='remove')        
#    Ohm_Nafion.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue = guess_I, pores = Ohm_Nafion_BC4_pores)        
#    ##-----------------------------------------------------------------------------
#    Ohm_Nafion.set_boundary_conditions(bctype='Neumann_rate_group', pores = Ohm_Nafion_BC5_pores, mode='remove')
#    Ohm_Nafion.set_boundary_conditions(bctype='Neumann_rate_group', bcvalue = -guess_I, pores = Ohm_Nafion_BC5_pores)
#    ##-----------------------------------------------------------------------------
#    # Run Simulation        
#    Ohm_Nafion.run(active_fluid=Nafion)
#    Ohm_Nafion.update()
#    ##=============================================================================
#    '''Step4_ Calculate the Tafel current based on the concentration at the catalyst layer'''
#    
#    V_e_CL = solid.get_data(prop='voltage',pores=P_CL_anode)
#    V_H_CL = Nafion.get_data(prop='voltage',pores=P_CL_anode)
#    delta_V = V_e_CL - V_H_CL 
#    C_H2 = hydrogen.get_data(prop='mole_fraction',pores=P_CL_anode)*hydrogen.get_data(prop='molar_density',pores=P_CL_anode)
#    I0_anode = 1.8e-2
#    T = hydrogen.get_data(prop='temperature',pores=P_CL_anode)
#    surface_area = 1e-2/(pn.num_pores('CL_anode'))
#    I = sp.sum(I0_anode*surface_area*(C_H2/C_ref)**1.5*sp.exp(-2*alpha_anode*F*delta_V/(R*T)))
#    counter += 1
#    print(counter,'---------------------')
#    print('guess:',guess_I)    
#    print('I:',I)
#    error = (sp.absolute(guess_I-I)/I)*100
#    if I>guess_I:
#        guess_I = (guess_I+I)/2
#    else:
#        guess_I = I + (guess_I-I)/2
##    current.append(I)
##current = sp.array(current,ndmin=1)
##import matplotlib.pylab as plt
##plt.plot(current,voltage)
#
#
##Deff = OpenPNM.Algorithms.EffectiveProperty(network=pn)
##Deff.setup(algorithm=Fickian_alg,fluid=air,conductance='diffusive_conductance',quantity='mole_fraction')
##a = Deff.run()
#
##------------------------------------------------------------------------------
#'''Export to VTK'''
##------------------------------------------------------------------------------
##OpenPNM.Visualization.Vtp.write(filename='test.vtp',fluids=[air,water],network=pn)
#vis = OpenPNM.Visualization.VTK(loglevel=10)
#vis.write(network=pn)
#
##------------------------------------------------------------------------------
#'''Perform Fickian Diffusion'''
##------------------------------------------------------------------------------
#alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=10, network=pn)
## Assign Dirichlet boundary conditions to top and bottom surface pores
#BC1_pores = pn.pores(labels=['top','boundary'],mode='intersection')
#alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=sp.log(1-0.6), pores=BC1_pores)
#
#BC2_pores = pn.pores(labels=['bottom','boundary'],mode='intersection')
#alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=sp.log(1-0.4), pores=BC2_pores)
#
## Updating data based on the result of Percolation Algorithms
#OP_1.update(Pc=1000)
## Run simulation
#alg.setup(fluid=air)
#alg.run()
#alg.update()


#Deff = OpenPNM.Algorithms.EffectiveProperty(network=pn)
#Deff.setup(algorithm=Fickian_alg,fluid=air,conductance='diffusive_conductance',quantity='mole_fraction')
#a = Deff.run()

#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
#OpenPNM.Visualization.VTK.write(filename='test.vtp',fluids=[air,water,Nafion,solid],network=pn)
vis = OpenPNM.Visualization.VTK()
vis.write(filename='test.vtp',fluids=[air,hydrogen,Nafion,solid],network=pn)






