import OpenPNM
import scipy as sp
import OpenPNM.Utilities.misc as misc


#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='Cell',loglevel=20)

Lc = sp.array([0.00005])
x = 30
y = 30
z = 5

pn.generate(divisions = [x, y, z], lattice_spacing = Lc,add_boundaries=True)



#==============================================================================
'''Build Geometry'''
#==============================================================================

#### GDL pores and throats for anode side
###----------------------------------------------------------------------------
geom = OpenPNM.Geometry.Toray090(name='GDL',network=pn,loglevel=10)
geom.set_locations(pores=pn.pores('internal'), throats='all')

boun = pn.add_geometry(name='boundary_geometry',subclass='Boundary')
boun.set_locations(pores=pn.pores('boundary'))

pn.regenerate_geometries()


#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(name='Air',network=pn, loglevel=50)
air.apply_conditions(temperature=350, pressure=200000)
###----------------------------------------------------------------------------
water = OpenPNM.Fluids.Water(name='Water',network=pn, loglevel=50)
water.apply_conditions(temperature=350, pressure=200000)

#Use Network's Fluid regeneration method
pn.regenerate_fluids()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================

#### Physics
###----------------------------------------------------------------------------
phys_air_GDL = OpenPNM.Physics.GenericPhysics(name='physics_air_GDL',network=pn, fluid=air,geometry=geom,loglevel=50)

phys_air_GDL.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')

#### Physics
###----------------------------------------------------------------------------
phys_water_GDL = OpenPNM.Physics.GenericPhysics(name='physics_water_GDL',network=pn, fluid=water,geometry=geom,loglevel=50)
phys_water_GDL.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')

#### Physics
###----------------------------------------------------------------------------
phys_air_boun = OpenPNM.Physics.GenericPhysics(name='physics_air_boun',network=pn, fluid=air,geometry=boun,loglevel=50)

phys_air_boun.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')

#### Physics 
###----------------------------------------------------------------------------
phys_water_boun = OpenPNM.Physics.GenericPhysics(name='physics_water_boun',network=pn, fluid=water,geometry=boun,loglevel=50)
phys_water_boun.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')

#Use Network's Physics regeneration method
pn.regenerate_physics()

water.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
air.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
#air.set_data(prop='hydraulic_conductance',throats=pn.throats('all'),data=1e-6)
##==============================================================================
#'''Begin Simulations'''
##==============================================================================
#'''PEM effective resistance '''
##==============================================================================
P_eff_alg = OpenPNM.Algorithms.StokesFlow(name='P_eff_alg',loglevel=10, network=pn)
P_eff_alg_BC1_pores =  pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
P_eff_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue = 400000, pores = P_eff_alg_BC1_pores)
P_eff_alg_BC2_pores =  pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
P_eff_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue = 200000, pores = P_eff_alg_BC2_pores)
P_eff_alg.setup(fluid=air)
P_eff_alg.run()
P_eff_alg.update()

Ps = P_eff_alg.pores(labels='pore.Dirichlet')
BCs = sp.unique(P_eff_alg['pore.bcval_Dirichlet'][Ps])
#Find flow through inlet face
Pn = pn.find_neighbor_pores(pores=P_eff_alg_BC1_pores,excl_self=True)
#Pn = Pn[-sp.in1d(Pn,pn.pores('GDL'))]
Ts = pn.find_connecting_throat(P_eff_alg_BC1_pores,Pn)
g = air['throat.hydraulic_conductance'][Ts]
s = air['throat.occupancy'][Ts]

inlets = sp.where(P_eff_alg['pore.bcval_Dirichlet']==sp.amax(BCs))[0]
outlets = sp.where(P_eff_alg['pore.bcval_Dirichlet']==sp.amin(BCs))[0]
Pin = inlets
Pout = outlets

#Fetch area and length of domain
A = pn.domain_area(face=Pin)
L = pn.domain_length(face_1=Pin,face_2=Pout)

xin = P_eff_alg['pore.pressure'][P_eff_alg_BC1_pores]
xout = P_eff_alg['pore.pressure'][Pn]
flow = g*s*(xin - xout)
D = sp.sum(flow)*L/A*air['pore.viscosity']/sp.absolute(BCs[0]-BCs[1])
#------------------------------------------------------------------------------


##------------------------------------------------------------------------------
#'''Export to VTK'''
##------------------------------------------------------------------------------
##OpenPNM.Visualization.VTK.write(filename='test.vtp',fluids=[air,water,Nafion,solid],network=pn)
#vis = OpenPNM.Visualization.VTK()
#vis.write(filename='test.vtp',network=pn)






