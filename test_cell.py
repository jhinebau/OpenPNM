

import OpenPNM
import scipy as sp
import run_function
import matplotlib.pylab as plt

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


X = pn.get_data(prop='coords',pores='all')[:,0]
Y = pn.get_data(prop='coords',pores='all')[:,1]
Z = pn.get_data(prop='coords',pores='all')[:,2]

'''Build Geometries'''
#==============================================================================
run_function.geom(pn)
#untitled2.geom(pn)
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
'''Build Fluids (air,solid,Nafion,hydrogen)'''
#==============================================================================
[air,solid,Nafion,hydrogen] = run_function.fluid(pn)
#[air,solid,Nafion,hydrogen] = untitled2.fluid(pn)
#==============================================================================
'''Build Physics'''
#==============================================================================

R = 8.314
F = 9.64870e4
x_ref = 1
alpha_anode = 0.5
alpha_cathode = 0.5
T = 353
I0_anode = 1.2246e-7 # A/cm2
I0_cathode = 1.2246e-7 #A/cm2
PEM_sigma = 50 #Siemens
porosity_CL = 0.5
C_ionomer_fraction = 0.4 
tort_CL = (porosity_CL)**1.5
PEM_Lc = 0.000025
CL_Lc = 0.000025
surf_area = 20  #cm2


run_function.physics(pn,PEM_Lc=PEM_Lc,PEM_sigma=PEM_sigma,C_ionomer_fraction=C_ionomer_fraction,CL_Lc=CL_Lc,porosity_CL=porosity_CL,tort_CL=tort_CL)
solid.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
Nafion.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
hydrogen.set_data(prop='occupancy',throats=pn.throats('all'),data=1)
air.set_data(prop='occupancy',throats=pn.throats('all'),data=1)

#==============================================================================
'''Begin Simulations'''
#==============================================================================

'''PEM effective resistance '''

T_temp = pn.find_neighbor_throats(pn.pores('CLa'))
random_val = sp.random.rand(len(pn.pores('CLa')))
max_val = max(random_val)
min_val = min(random_val)
delta = 5
Reff = sp.zeros((delta+1))
inv_percentage = sp.zeros((delta+1))

occupancy = sp.copy(Nafion['throat.occupancy'])
conductance = sp.copy(Nafion['throat.electronic_conductance'])

for j in sp.r_[:(delta+1)]:

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
    
    temp_val = min_val + (max_val - min_val)/delta*j/1.5
    pinv = pn.pores('CLa')[random_val<=temp_val]

    PEM_Reff_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue = 1e-20, pores = pinv)
    t = pn.find_neighbor_throats(pinv,mode='not_intersection')
    Nafion['throat.electronic_conductance'][t] = 1e-100
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
    i = (g*s+g*(1-s)/1e6)*(xin - xout)
    Reff[j] = sp.absolute((BCs[0]-BCs[1])/sp.sum(i))
    inv_percentage[j] = (len(pinv)+1e-50)/len(pn.pores('CLa'))
    Nafion['throat.electronic_conductance'] = conductance
    Nafion['throat.occupancy'][pn.find_neighbor_throats(pinv)] = 0

#------------------------------------------------------------------------------
plt.plot(Reff,inv_percentage)
plt.xlabel('Resistance')
plt.ylabel('Saturation at GDL-CL interface')



##------------------------------------------------------------------------------
#

current =[]
voltage= []
Power = []
IR = []

hydrogen['throat.diffusive_conductance'][T_PEM_CL_a] = 1e-90
air['throat.diffusive_conductance'][T_PEM_CL_c] = 1e-90

l1 =sp.r_[0.15:3.6:0.1]
l2 = sp.r_[3.6:4.4:0.01]
i_list = sp.concatenate((l1,l2))

for I in i_list:
    
    current.append(I)
    S_H2 = I*surf_area*1e-4/(2*F)
    S_O2 = I*surf_area*1e-4/(4*F)
    
    #==============================================================================
    '''Necessary transport algorithms and the fixed boundary conditions'''
    
    ##=============================================================================
    '''Hydrogen Fickian Diffusion'''
    ##=============================================================================
    
    Fick_hydrogen = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_hydrogen',loglevel=50, network=pn)
    #------------------------------------------------------------------------------
    Fick_hydrogen_BC1_pores = pn.pores('GDLa')[Z[P_GDL_a]==min(Z[P_GDL_a])]
    Fick_hydrogen.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.3, pores=Fick_hydrogen_BC1_pores)
    #------------------------------------------------------------------------------
    Fick_hydrogen_BC2_pores = pn.pores('CLa')
    
    Fick_hydrogen.set_boundary_conditions(bctype='Neumann_group', bcvalue= S_H2, pores=Fick_hydrogen_BC2_pores)
    

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
    Fick_air.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.3, pores=Fick_air_BC1_pores)
    #------------------------------------------------------------------------------
    
    Fick_air_BC2_pores = pn.pores('CLc')
    
    Fick_air.set_boundary_conditions(bctype='Neumann_group', bcvalue= S_O2, pores=Fick_air_BC2_pores)
    
    #------------------------------------------------------------------------------
    air.set_nan_value(element='throat',prop='diffusive_conductance',value=1e-40)
    Fick_air.setup(fluid=air)
    Fick_air.run()
    Fick_air.update()
    #
    
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
    etta_c = R*T/(F*alpha_cathode)*sp.log(I_local_cathode/(I0_cathode*surf_area)*((1/X_O2)**1.5))
    
    Vcell = 1.233-sp.average(sp.absolute(etta_a))-sp.average(sp.absolute(etta_c))-I*Reff*1e-4*surf_area
    voltage.append(Vcell)
    Power.append(Vcell*I)
    IR.append(I*Reff*1e-4*surf_area)
##=============================================================================
##=============================================================================

voltage = sp.array(voltage,ndmin=1)
current = sp.array(current,ndmin=1)
Power = sp.array(Power,ndmin=1)
IR = sp.array(IR,ndmin=1)
import matplotlib.pylab as plt
plt.figure(2)
plt.plot(current,voltage)
plt.xlabel('Current density (A/cm^2)')
plt.ylabel('Voltage (V)')
plt.ylim([0,1.5])
plt.figure(3)
plt.plot(current,Power)
plt.xlabel('Current (A/cm^2)')
plt.ylabel('Power density (W/cm^2)')
plt.ylim([0,2.5])
plt.figure(4)
plt.plot(current,IR)
plt.xlabel('Current (A/cm^2)')
plt.ylabel('IR (V)')
plt.ylim([0,1])
#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
#OpenPNM.Visualization.VTK.write(filename='test.vtp',fluids=[air,water,Nafion,solid],network=pn)
vis = OpenPNM.Visualization.VTK()
vis.write(filename='test_cell.vtp',fluids=[air,hydrogen,Nafion,solid],network=pn)

