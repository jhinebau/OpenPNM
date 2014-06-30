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
geom = OpenPNM.Geometry.Toray090(network=pn,name='geom')
geom.set_locations(pores=pn.pores('all'),throats='all')

pn.regenerate_geometries()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, loglevel=30,name='air')
air.apply_conditions(temperature=350, pressure=200000)
air.add_property(prop='thermal_conductivity',model='constant',value=0.0262)
air.add_property(prop='electrical_conductivity',model='constant',value=1)

water = OpenPNM.Fluids.Water(network=pn,loglevel=30,name='water')
water.add_property(prop='diffusivity',prop_name='DAB',model='constant',value=5e-12)

#Use Network's Fluid regeneration method
pn.regenerate_fluids()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.BasePhysics(network=pn, fluid=water,geometry=geom,name='physwater')

phys_air = OpenPNM.Physics.BasePhysics(network=pn, fluid=air,geometry=geom,name='physair')

#Use Network's Physics regeneration method
pn.regenerate_physics()
#
##==============================================================================
#'''Begin Simulations'''
##==============================================================================
#'''Perform a Drainage Experiment (OrdinaryPercolation)'''
##------------------------------------------------------------------------------
#OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=30,network=pn)
#a = pn.pores(labels=['bottom','boundary'],mode='intersection')
#OP_1.setup(invading_fluid=water,defending_fluid=air,inlets=a,npts=20)
#OP_1.run()
##OP_1.plot_drainage_curve()
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
#OpenPNM.Visualization.Vtp.write(filename='test.vtp',fluids=[air,water],network=pn)
vis = OpenPNM.Visualization.VTK()
vis.write(network=pn)






