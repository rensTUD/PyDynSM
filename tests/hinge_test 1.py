# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 18:24:14 2025

@author: thana
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# %% add package to sys and import

import pydynsm as PDM

Assembler = PDM.Assembler

s1 = Assembler('test')
# %% inputs
l1 = 8
l2 = 22
l3 = 6

# material props
Es = 210e9  # Elasticity modulus steel
Ec = 3.1e10  # Elasticity modulus concrete
rhoc = 2500  # density concrete
rhos = 7800  # density steel
Ab = 2  # beam cross sectional area
Ib = 0.00395  # beam 2nd mmt of intertial
rc = 9.5 * 0.001
Ac = np.pi * rc**2
Ic = 0.25 * np.pi * rc**4
omega_f = 100

# %%% nodes
node1 = s1.CreateNode(0, 0, dof_config=['z', 'phi_y'])
node2 = s1.CreateNode(l2, 0, dof_config=['z', 'phi_y'])
# node3 = s1.CreateNode(1.8 * l2, 0)
node3 = s1.CreateNode(2*l2, 0, dof_config=['z', 'phi_y'])

# s1.PlotStructure()

# %%% elements
# beam
beams = []

# Define node pairs for beam elements
beam_node_pairs = [
    (node1, node2),
    (node2, node3)]

# Loop to create and append beam elements
for pair in beam_node_pairs:
    beam = s1.CreateElement(list(pair))
    beams.append(beam)  # Append beam element to list

# %% constraint nodes
node1.fix_node('z')
node2.fix_node('z')
node3.fix_node('z')

# beams[0].free_dof(node2, 'phi_y')
beams[0].decouple_dof(node2, 'phi_y')
beams[1].decouple_dof(node2, 'phi_y')


# %% assign sections
for beam in beams:
    beam.SetSection('EulerBernoulli Beam', {'E': Ec, 'A': Ab, 'rho': rhoc, 'Ib': Ib, 'Wb': Ib})
    # beam.SetSection('Rod', {'E': Ec, 'A': Ac, 'rho': rhoc})

# s1.PlotStructure(plot_elements=True)

# %% loading definition

for beam in beams[:]:
    beam.AddDistributedLoad(z=1e03)

s1.run_connectivity()

# %% Get the global stiffness and force matrices

# K_global = s1.GlobalStiffness(omega_f)
# F_global = s1.GlobalForce(omega_f)
Kc_global = s1.GlobalConstrainedStiffness(omega_f)
Fc_global = s1.GlobalConstrainedForce(omega_f)


# %% solve for free node displacements

u_free = s1.SolveUfree(Kc_global, Fc_global)
# f_supp = s1.SupportReactions(s1.GlobalStiffness(omega_f), u_free, s1.GlobalForce(omega_f))

# %% post-processing
u_elem = s1.FullDisplacement(u_free)
disp = s1.ElementDisplacements(u_elem, omega_f, num_points=100)
force = s1.ElementForces(u_elem, omega_f, num_points=200)
s1.PlotMoments(force, scale=1e-3)
# s1.PlotShearforces(force, scale=1e-3)
s1.PlotElementDisplacements(disp, scale=1e05)
# s1.PlotAxialforces(force, scale=1e-5)
