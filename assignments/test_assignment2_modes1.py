# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:59:03 2025

@author: GREY
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA
import scipy.optimize as opt
import pydynsm as PDM

Assembler = PDM.Assembler
s1 = Assembler('bridge')

# functions for root-finder
def BVP(ww):
    return s1.GlobalConstrainedStiffness(ww)

def det_func(ww):
  dets = np.linalg.det(BVP(ww)/1e8)
  # dets = np.heaviside(dets.real, 1) * np.log(abs(dets.real)) - np.heaviside(-dets.real, 0) * np.log(abs(dets.real))
  return dets

def find_eigen_frequencies(omega):
    Det_M = np.array([det_func(ww) for ww in omega])
    omega_initial = omega[np.where(np.isclose(abs(np.diff(np.angle(Det_M)))/np.pi,1, atol=.1))[0]]
    # omega_initial = omega[np.where(np.isclose(Det_M,1,atol=.1))[0]]
    omega_m = []
    for ww in omega_initial:
        omega_m.append(opt.newton(det_func,ww).real)
    return np.unique(omega_m)

def remove_close_roots(roots, tol=1e-6):
    '''
    The function removes close roots from the array of root-finders
    '''
    roots = np.sort(roots)
    filtered_roots = [roots[0]]

    for i in range(1, len(roots)):
        if not np.isclose(roots[i], filtered_roots[-1], atol=tol):
            filtered_roots.append(roots[i])

    return np.array(filtered_roots)
# %% inputs
l1 = 4
l2 = 11.75
l3 = 5.33

# material props
E = 210e9  # Elasticity modulus steel
rho =  7800
# for beam
I_b = 0.0395
A_b = 0.5

# for cable
r_c = 9.5*0.005
A_c = np.pi*r_c**2
I_c = 0.25*np.pi*r_c**4

# for truss 
A_t = 0.1*A_b # truss erea
I_t = I_c


k_d = 1e6
# c_d = 1e6
c_d = 0
# %% sturctural plot
# %%% nodes
node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(l2,0)
node3 = s1.CreateNode(2*l2,0)
node4 = s1.CreateNode(3*l2,0)
node5 = s1.CreateNode(4*l2,0)
node6 = s1.CreateNode(4*l2+l1,0)
node7 = s1.CreateNode(4*l2+l1,-l1,dof_config = [['x'],['z']])
node8 = s1.CreateNode(4*l2+l1+l3,0)
node9 = s1.CreateNode(4*l2+l1+2*l3,0,)
node10 = s1.CreateNode(4*l2+l1+3*l3,0)
node11 = s1.CreateNode(4*l2+l1+4*l3,0)
node12 = s1.CreateNode(4*l2+l1,l1,dof_config = [['x'],['z']])
node13 = s1.CreateNode(4*l2,l1,dof_config = [['x'],['z']])
node14 = s1.CreateNode(4*l2+l1,2*l1,dof_config = [['x'],['z']])
node15 = s1.CreateNode(4*l2,2*l1,dof_config = [['x'],['z']])
node16 = s1.CreateNode(4*l2+l1,3*l1,dof_config = [['x'],['z']])
node17 = s1.CreateNode(4*l2,3*l1,dof_config = [['x'],['z']])
node18 = s1.CreateNode(4*l2,4*l1,dof_config = [['x'],['z']])
s1.PlotStructure()
# %%% elements
# beam
beams = []

# Define node pairs for beam elements
beam_node_pairs = [
    (node1, node2),
    (node2, node3),
    (node3, node4),
    (node4, node5),
    (node5, node6),
    (node6, node8),
    (node8, node9),
    (node10, node11)
]

# Loop to create and append beam elements
for pair in beam_node_pairs:
    beam = s1.CreateElement(list(pair))
    beams.append(beam)  # Append beam element to list

# Truss
trusses = []

# Define node pairs for truss elements
truss_node_pairs = [
    (node5, node7),(node6, node7),
    (node5, node12), (node6, node12),
    (node12, node13), (node12, node14),
    (node12, node15), (node13, node15),
    (node13, node5), (node14, node15),
    (node14, node16), (node15, node16),
    (node15, node17), (node16, node17),
    (node16, node18), (node17, node18)
]

# Loop to create and append truss elements
for pair in truss_node_pairs:
    truss = s1.CreateElement(list(pair))
    trusses.append(truss)

# front cables
front_cables = []
front_cable_node_pairs = [(node18, node2), (node3, node17),
(node4, node15)]

for pair in front_cable_node_pairs:
    cable = s1.CreateElement(list(pair))
    front_cables.append(cable)
# back cables    
back_cables = []
back_cable_node_pairs = [(node9, node14),(node8, node12),(node10, node16)]
for pair in back_cable_node_pairs:
    cable = s1.CreateElement(list(pair))
    back_cables.append(cable)

# %% constraint nodes
node1.fix_node('z')
node7.fix_node('x','z')

# %% assign sections
for beam in beams:
    beam.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0, 'Wb':I_b})
    beam.SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0})

beams[-1].SetSection('EulerBernoulli Beam with foundation', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0, 'kd':k_d,'cd':c_d})
beams[-1].SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0})
for cable in front_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':0})
    # cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':0})
for cable in back_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':0})
    # cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':0})
for truss in trusses:
    truss.SetSection('Rod', {'E': E, 'A':A_t, 'rho':rho,'ksi':0})
    # truss.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_t, 'rho':rho, 'Ib':I_t,'ksi':0})
s1.run_connectivity()
p = 
nodes[1].add_load(z=p)