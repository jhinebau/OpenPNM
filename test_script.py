import numpy as np

import OpenPNM
import pytest

def test_cubic_standard_call():
  pn = OpenPNM.Network.Cubic(shape=[3,4,5])
  np.testing.assert_almost_equal(pn['pore.coords'][0], [0.5,0.5,0.5])

def test_cubic_optional_call():
  image = np.random.rand(30,40)
  pn = OpenPNM.Network.Cubic.from_image(image=image)
  np.testing.assert_almost_equal(image, pn.asarray(pn['pore.values']))

def test_linear_solvers():
  pn = OpenPNM.Network.Cubic([1,40,30], spacing=0.0001)
  geom = OpenPNM.Geometry.Toray090(network=pn,pores=pn.pores(),throats=pn.throats())
  air = OpenPNM.Phases.Air(network=pn)
  phys_air = OpenPNM.Physics.Standard(network=pn,phase=air,pores=pn.pores(),throats=pn.throats())

  BC1_pores = pn.pores(labels=['left'])
  BC2_pores = pn.pores(labels=['right'])

  alg_1 = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
  alg_1.set_boundary_conditions(bctype='Dirichlet', bcvalue=1, pores=BC1_pores)
  alg_1.set_boundary_conditions(bctype='Dirichlet', bcvalue=0, pores=BC2_pores)
  alg_1.run(phase=air)

  alg_2 = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
  alg_2.set_boundary_conditions(bctype='Neumann', bcvalue = -1e-11, pores=BC1_pores)
  alg_2.set_boundary_conditions(bctype='Dirichlet', bcvalue=0, pores=BC2_pores)
  alg_2.run(phase=air)

  alg_3 = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
  alg_3.set_boundary_conditions(bctype='Neumann_group', bcvalue = -3e-10, pores=BC1_pores)
  alg_3.set_boundary_conditions(bctype='Dirichlet', bcvalue=0, pores=BC2_pores)
  alg_3.run(phase=air)

  print( alg_1['pore.mole_fraction'][BC1_pores] )
  print( alg_2['pore.mole_fraction'][BC1_pores] )
  print( alg_3['pore.mole_fraction'][BC1_pores] )

def test_add_boundary():
  pn = OpenPNM.Network.Cubic([5,5,5])
  pn.add_boundaries()

  keys_expected = {'pore.all', 'pore.back', 'pore.back_boundary', 'pore.bottom',
  'pore.bottom_boundary','pore.boundary', 'pore.coords', 'pore.front',
  'pore.front_boundary', 'pore.left', 'pore.left_boundary', 'pore.right',
  'pore.right_boundary', 'pore.top', 'pore.top_boundary', 'pore.values',
  'throat.all', 'throat.back', 'throat.back_boundary', 'throat.bottom',
  'throat.bottom_boundary', 'throat.boundary', 'throat.conns', 'throat.front',
  'throat.front_boundary', 'throat.left', 'throat.left_boundary', 'throat.right',
  'throat.right_boundary', 'throat.top', 'throat.top_boundary'}

  keys_found = set(pn.keys())

  symmetric_diff = keys_found ^ keys_expected
  assert not symmetric_diff

def test_rectilinear_integrity():
  R = np.random.rand(3,4,5)
  # prune the network way
  network = OpenPNM.Network.Cubic.from_image(R)
  network.trim( R.ravel() <= R.mean() )
  O = network.asarray(network['pore.values'])
  # what it would look like normally
  M = np.where(R > R.mean(), R, 0)
  np.testing.assert_almost_equal(M, O)
  
def test_open_air_diffusivity():
    pn = OpenPNM.Network.Cubic([5,5,5], spacing=1, name='net', loglevel=30)
    pn.add_boundaries()
    Ps = pn.pores('boundary',mode='difference')
    Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
    geom = OpenPNM.Geometry.Cube_and_Cuboid(network=pn,pores=Ps,throats=Ts)
    geom['pore.diameter'] = 0.999999
    geom['throat.diameter'] = 0.999999
    geom.regenerate(['pore.diameter','throat.diameter'],mode='exclude')
    Ps = pn.pores('boundary')
    Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
    boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts,shape='cubes')
    air = OpenPNM.Phases.Air(network=pn,name='air')
    Ps = pn.pores()
    Ts = pn.throats()
    phys_air = OpenPNM.Physics.Standard(network=pn,phase=air,pores=Ps,throats=Ts,name='phys_air')
    BC1_pores = pn.pores(labels=['top_boundary'])
    BC2_pores = pn.pores(labels=['bottom_boundary'])
    print('length =',round(pn.domain_length(BC1_pores,BC2_pores)))
    print('area = ',round(pn.domain_area(BC1_pores),2),round(pn.domain_area(BC2_pores),2))
    Diff = OpenPNM.Algorithms.FickianDiffusion(loglevel=30, network=pn)
    # Assign Dirichlet boundary conditions to top and bottom surface pores
    Diff.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    Diff.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
    Diff.run(phase=air)
    Diff.update()
    Diff_deff = Diff.calc_eff_diffusivity()/np.mean(air['pore.diffusivity'])
    assert np.round(Diff_deff,3) == 1

if __name__ == '__main__':
  pytest.main()
