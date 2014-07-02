import OpenPNM
import scipy as sp
import OpenPNM.Utilities.misc as misc


#==============================================================================
'''Build Topological Network'''
#==============================================================================
#pn = OpenPNM.Network.Cubic2(loglevel=10,name='net')
#
#domains = {}
#domains['GDLa'] = {'region': [[0,0,0],[0.0005,0.0005,0.0003]], 'Lc': 0.00005}
#domains['CLa']  = {'region': [[0,0,0.0003],[0.0005,0.0005,0.00035]], 'Lc': 0.00001}
#domains['PEM']  = {'region': [[0,0,0.00035],[0.0005,0.0005,0.0004]], 'Lc': 0.00001}
#domains['CLc']  = {'region': [[0,0,0.0004],[0.0005,0.0005,0.00045]], 'Lc': 0.00001}
#domains['GDLc'] = {'region': [[0,0,0.00045],[0.0005,0.0005,0.00075]], 'Lc': 0.00005}
#
#pn.generate(domains=domains)
#pn.save_object('big_net.npz')
pn = OpenPNM.Network.Cubic(loglevel=10,name='net')
pn.generate(domain_size=[0.0005,0.0005,0.00075],lattice_spacing=[0.000025],add_boundaries=False)

GDLa = 0.00030
CLa  = 0.00035
PEM  = 0.00040
CLc  = 0.00045
GDLc = 0.00075

temp = pn['pore.coords'][:,2]<=GDLa
pn['pore.GDLa'] = sp.zeros_like(pn['pore.all'],dtype=bool)
pn['pore.GDLa'][temp] = True
pn['throat.GDLa'] = sp.zeros_like(pn['throat.all'],dtype=bool)
pn['throat.GDLa'][pn.find_neighbor_throats(temp)] = True

temp = (pn['pore.coords'][:,2]>GDLa) * (pn['pore.coords'][:,2]<=CLa)
pn['pore.CLa'] = sp.zeros_like(pn['pore.all'],dtype=bool)
pn['pore.CLa'][temp] = True
pn['throat.CLa'] = sp.zeros_like(pn['throat.all'],dtype=bool)
pn['throat.CLa'][pn.find_neighbor_throats(temp)] = True

temp = (pn['pore.coords'][:,2]>CLa) * (pn['pore.coords'][:,2]<=PEM)
pn['pore.PEM'] = sp.zeros_like(pn['pore.all'],dtype=bool)
pn['pore.PEM'][temp] = True
pn['throat.PEM'] = sp.zeros_like(pn['throat.all'],dtype=bool)
pn['throat.PEM'][pn.find_neighbor_throats(temp)] = True

temp = (pn['pore.coords'][:,2]>PEM) * (pn['pore.coords'][:,2]<=CLc)
pn['pore.CLc'] = sp.zeros_like(pn['pore.all'],dtype=bool)
pn['pore.CLc'][temp] = True
pn['throat.CLc'] = sp.zeros_like(pn['throat.all'],dtype=bool)
pn['throat.CLc'][pn.find_neighbor_throats(temp)] = True

temp = (pn['pore.coords'][:,2]>CLc) * (pn['pore.coords'][:,2]<=GDLc)
pn['pore.GDLc'] = sp.zeros_like(pn['pore.all'],dtype=bool)
pn['pore.GDLc'][temp] = True
pn['throat.GDLc'] = sp.zeros_like(pn['throat.all'],dtype=bool)
pn['throat.GDLc'][pn.find_neighbor_throats(temp)] = True

R = 8.314
F = 9.64870e4
x_ref = 1
alpha_anode = 0.5
alpha_cathode = 0.5
T = 353
I0_anode = 1e3
I0_cathode = 1.8e-2
active_area = 5e5 #m2/m3
PEM_sigma = 1
porosity_CL = 0.5
C_ionomer_fraction = 0.4 
tort_CL = (porosity_CL)**-0.5
PEM_Lc = 0.00001
PEM_V = PEM_Lc**3
CL_Lc = 0.00001

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

sigma_val_CL = PEM_sigma*C_ionomer_fraction**1.5*(PEM_Lc**2/PEM_Lc)

phys_Nafion_CL_anode.add_property(prop='electronic_conductance', model='constant',value=sigma_val_CL)

#### Physics for Nafion in CL_cathode
###----------------------------------------------------------------------------
phys_Nafion_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_CL_cathode',network=pn, fluid=Nafion,geometry=geom_CL_cathode,loglevel=50)
phys_Nafion_CL_cathode.add_property(prop='electronic_conductance', model='constant',value=sigma_val_CL)

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

sigma_val_PEM = PEM_sigma*(PEM_Lc**2/PEM_Lc)

phys_Nafion_PEM.add_property(prop='electronic_conductance',  model='constant',value=sigma_val_PEM)

#### Physics for Nafion in PEM_CL_anode interface
###----------------------------------------------------------------------------
phys_Nafion_PEM_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_PEM_CL_anode',network=pn, fluid=Nafion,geometry=geom_PEM_CL_anode,loglevel=50)
phys_Nafion_PEM_CL_anode.add_property(prop='electronic_conductance', model='constant',value=sigma_val_PEM)

#### Physics for Nafion in PEM_CL_cathode interface
###----------------------------------------------------------------------------
phys_Nafion_PEM_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_Nafion_PEM_CL_cathode',network=pn, fluid=Nafion,geometry=geom_PEM_CL_cathode,loglevel=50)
phys_Nafion_PEM_CL_cathode.add_property(prop='electronic_conductance',model='constant',value=sigma_val_PEM)

#### Physics for hydrogen in GDL_anode
###----------------------------------------------------------------------------
phys_hydrogen_GDL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_GDL_anode',network=pn, fluid=hydrogen,geometry=geom_GDL_anode,loglevel=50)
phys_hydrogen_GDL_anode.add_property(prop='diffusive_conductance', model='bulk_diffusion')

#### Physics for hydrogen in CL_anode
###----------------------------------------------------------------------------
phys_hydrogen_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_CL_anode',network=pn, fluid=hydrogen,geometry=geom_CL_anode,loglevel=50)

DH_val_CL = hydrogen['pore.molar_density']*hydrogen['pore.diffusivity']*(porosity_CL/tort_CL)*(CL_Lc**2/CL_Lc)

phys_hydrogen_CL_anode.add_property(prop='diffusive_conductance', model= 'constant', value = DH_val_CL)

#### Physics for hydrogen in GDL_CL_anode interface
###----------------------------------------------------------------------------
phys_hydrogen_GDL_CL_anode = OpenPNM.Physics.GenericPhysics(name='physics_hydrogen_GDL_CL_anode',network=pn, fluid=hydrogen,geometry=geom_GDL_CL_anode,loglevel=50)
phys_hydrogen_GDL_CL_anode.add_property(prop='diffusive_conductance', model= 'constant', value = DH_val_CL)

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
DO_val_CL = air['pore.molar_density']*air['pore.diffusivity']*(porosity_CL/tort_CL)*(CL_Lc**2/CL_Lc)
phys_air_CL_cathode.add_property(prop='diffusive_conductance', model='constant', value = DO_val_CL)

#### Physics for air in GDL_CL_cathode interface
###----------------------------------------------------------------------------
phys_air_GDL_CL_cathode = OpenPNM.Physics.GenericPhysics(name='physics_air_GDL_CL_cathode',network=pn, fluid=air,geometry=geom_GDL_CL_cathode,loglevel=50)
phys_air_GDL_CL_cathode.add_property(prop='diffusive_conductance', model='constant', value = DO_val_CL)

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
#==============================================================================
PEM_Reff_alg = OpenPNM.Algorithms.OhmicConduction(name='PEM_Reff_alg',loglevel=10, network=pn)
PEM_Reff_alg_BC1_pores =  pn.pores('CLa')[Z[P_CL_a]==min(Z[P_CL_a])]
PEM_Reff_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue = 0.8, pores = PEM_Reff_alg_BC1_pores)
PEM_Reff_alg_BC2_pores = pn.pores('CLc')[Z[P_CL_c]==max(Z[P_CL_c])] 
PEM_Reff_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue = 0.4, pores = PEM_Reff_alg_BC2_pores)
Nafion.set_nan_value(element='throat',prop='electronic_conductance',value=1e-40)
PEM_Reff_alg.setup(fluid=Nafion)
PEM_Reff_alg.run()
PEM_Reff_alg.update()

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


current =[]
voltage= []
surf_area = active_area*PEM_V

hydrogen['throat.diffusive_conductance'][T_PEM_CL_a] = 1e-90
air['throat.diffusive_conductance'][T_PEM_CL_c] = 1e-90

for I in [0.0001]:
    
    current.append(I)
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
    Fick_hydrogen.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.5, pores=Fick_hydrogen_BC1_pores)
    #------------------------------------------------------------------------------
    Fick_hydrogen_BC2_pores = pn.pores('CLa')
    Fick_hydrogen.set_boundary_conditions(bctype='Neumann_group', bcvalue= S_H2, pores=Fick_hydrogen_BC2_pores)
    
    #------------------------------------------------------------------------------
    hydrogen.set_nan_value(element='throat',prop='diffusive_conductance',value=1e-40)
    
    Fick_hydrogen.setup(fluid=hydrogen)
    Fick_hydrogen.run()
    Fick_hydrogen.update()
    
    ##=============================================================================
    '''Air Fickian Diffusion'''
    ##=============================================================================
    Fick_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_air',loglevel=50, network=pn)
    #------------------------------------------------------------------------------
    Fick_air_BC1_pores = pn.pores('GDLc')[Z[P_GDL_c]==max(Z[P_GDL_c])]
    Fick_air.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.1, pores=Fick_air_BC1_pores)
    #------------------------------------------------------------------------------
    Fick_air_BC2_pores = pn.pores('CLc')
    Fick_air.set_boundary_conditions(bctype='Neumann_group', bcvalue= S_O2, pores=Fick_air_BC2_pores)
    air.set_nan_value(element='throat',prop='diffusive_conductance',value=1e-40)
    Fick_air.setup(fluid=air)
    Fick_air.run()
    Fick_air.update()
    
    ##=============================================================================
    '''Local mass consumption anode'''
    ##=============================================================================
    mass_anode = sp.absolute(Fick_hydrogen.rate(pn.pores('CLa'),mode='single'))
    ##==============================================================================
    ##=============================================================================
    '''Local mass consumption cathode'''
    ##=============================================================================
    mass_cathode = sp.absolute(Fick_air.rate(pn.pores('CLc'),mode='single'))
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
    X_H2 = hydrogen.get_data(prop='mole_fraction',pores=P_CL_a)
    etta_a = R*T/(F*alpha_anode)*sp.log(I_local_anode/(I0_anode*surf_area)*((1/X_H2)**1.5))
    
    ##=============================================================================
    '''etta cathode'''
    ##=============================================================================
    X_O2 = air.get_data(prop='mole_fraction',pores=P_CL_c)
    etta_c = -R*T/(F*alpha_cathode)*sp.log(I_local_cathode/(I0_cathode*surf_area)*((1/X_O2)**1.5))
    
    Vcell = 1.233-sp.average(etta_a)-sp.average(etta_c)-I*Reff
    voltage.append(Vcell)


#import matplotlib.pylab as plt
#plt.plot(current,voltage)


#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
#OpenPNM.Visualization.VTK.write(filename='test.vtp',fluids=[air,water,Nafion,solid],network=pn)
vis = OpenPNM.Visualization.VTK()
vis.write(filename='test.vtp',network=pn)






