# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:59:03 2025

@author: GREY
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# %% add package to sys and import 

import pydynsm as PDM

Assembler = PDM.Assembler

s1 = Assembler('bridge',analysis_type='new')
# %% inputs
l1 = 8
l2 = 22
l3 = 6

# material props
Es = 210e9  # Elasticity modulus steel
Ec = 3.1e10 # Elasticity modulus concrete
rhoc = 2500 # density concrete
rhos = 7800 # density steel
Ab = 2 # beam cross sectional area 
Ib = 0.00395 # beam 2nd mmt of intertial
rc = 9.5*0.001
Ac = np.pi*rc**2
Ic = 0.25*np.pi*rc**4
omega_f = 0.001


# %% sturctural plot
# %%% nodes
node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(l2,0)
node3 = s1.CreateNode(2*l2,0)
node4 = s1.CreateNode(3*l2,0)
node5 = s1.CreateNode(4*l2,0)
node6 = s1.CreateNode(4*l2+l1,0)
node7 = s1.CreateNode(4*l2+l1,-l1)
node8 = s1.CreateNode(4*l2+l1+l3,0)
node9 = s1.CreateNode(4*l2+l1+2*l3,0)
node10 = s1.CreateNode(4*l2+l1+3*l3,0)
node11 = s1.CreateNode(4*l2+l1+4*l3,0)
node12 = s1.CreateNode(4*l2+l1,l1)
node13 = s1.CreateNode(4*l2,l1)
node14 = s1.CreateNode(4*l2+l1,2*l1)
node15 = s1.CreateNode(4*l2,2*l1)
node16 = s1.CreateNode(4*l2+l1,3*l1)
node17 = s1.CreateNode(4*l2,3*l1)
node18 = s1.CreateNode(4*l2,4*l1)
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
    (node6, node8),
    (node8, node9),
    (node9, node10),
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
    (node5, node6),
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
    beam.SetSection('EulerBernoulli Beam', {'E': Ec, 'A':Ab, 'rho':rhoc, 'Ib':Ib, 'Wb': Ib})
    beam.SetSection('Rod', {'E': Ec, 'A':Ac, 'rho':rhoc})
for cable in front_cables:
    cable.SetSection('Rod', {'E': Es, 'A':Ac, 'rho':rhos})
    cable.SetSection('EulerBernoulli Beam', {'E': Es, 'A':Ac, 'rho':rhos, 'Ib':Ic, 'Wb': Ic})
for cable in back_cables:
    cable.SetSection('Rod', {'E': Es, 'A':Ac, 'rho':rhos})
    cable.SetSection('EulerBernoulli Beam', {'E': Es, 'A':Ac, 'rho':rhos, 'Ib':Ic, 'Wb': Ic})
for truss in trusses:
    truss.SetSection('Rod', {'E': Es, 'A':Ac, 'rho':rhos})
    truss.SetSection('EulerBernoulli Beam', {'E': Es, 'A':Ac, 'rho':rhos, 'Ib':Ic, 'Wb': Ic})

# s_cable.SetSection('Rod', {'E': E, 'A':A, 'rho':rho})

s1.PlotStructure(plot_elements=True)

# %% loading definition
p_x = lambda omega: 2.9e3 if omega == omega_f else 0

# wind load

w_x1 = lambda omega: 20e3 if omega == omega_f else 0
w_x2 = lambda omega: 20e3 if omega == omega_f else 0
w_x3 = lambda omega: 20e3 if omega == omega_f else 0
w_x4 = lambda omega: 20e3 if omega == omega_f else 0

for beam in beams:
    beam.AddDistributedLoad(z=p_x)

s1.run_connectivity()

# %% Get the global stiffness and force matrices

K_global = s1.GlobalStiffness(omega_f)
F_global = s1.GlobalForce(omega_f)
Kc_global = s1.GlobalConstrainedStiffness(omega_f)
Fc_global = s1.GlobalConstrainedForce(omega_f)

# %% solve for free node displacements
u_free = s1.SolveUfree(Kc_global, Fc_global)
f_supp = s1.SupportReactions(s1.GlobalStiffness(omega_f), u_free, s1.GlobalForce(omega_f))

# %% post-processing
u_elem = s1.FullDisplacement(u_free)
disp = s1.ElementDisplacements(u_elem, omega_f,num_points=20)
force = s1.ElementForces(u_elem, omega_f,num_points=100)
s1.PlotMoments(force,scale = 1e-5)
s1.PlotElementDisplacements(disp,scale=1)
s1.PlotAxialforces(force,scale = 1e-5)