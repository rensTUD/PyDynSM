# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:23:06 2024

@author: rensv
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# %% add package to sys and import 

import pydynsm as PDM


Assembler = PDM.Assembler

# %% Example

# %%% Initialise an assembler with your project name

s2 = Assembler('beam')

# %%% Parameters
E = 210e9
EA = 7e6
A = EA/E
EI = 1.5 * 7e5
I = EI/E
rhoA = 1e03 
rho = rhoA/A
L  = 1
ksi = 0
omega_f = 1



# %%% Create nodes from the assembler

node1 = s2.CreateNode(0,0)
node2 = s2.CreateNode(L,0)
node3 = s2.CreateNode(2*L,0)
# node4 = s2.CreateNode(0,-L)

elem1 = s2.CreateElement([node1,node2])
elem2 = s2.CreateElement([node2,node3])
# elem3 = s2.CreateElement([node2,node4])


elem1.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'ksi':0, 'Wb': I})
elem2.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'ksi':0, 'Wb': I})
# elem3.SetSection('Rod', {'E': E, 'A':A, 'rho':rho,'ksi':0})
# elem3.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'ksi':0})

# %%% Now set the constraints directly onto the nodes

node1.fix_node('x','z')
node2.fix_node('x')
node3.fix_node('x','z')


# %%%% test adding a load to the node

q_z = lambda Omega: 1e6 if Omega == omega_f else 0
elem1.AddDistributedLoad(z=q_z)
elem2.AddDistributedLoad(z=q_z)
s2.run_connectivity()

# %%% Get the global stiffness and force matrices

K_global = s2.GlobalStiffness(omega_f)
F_global = s2.GlobalForce(omega_f)

# %%% get the global constrained stiffness and force

Kc_global = s2.GlobalConstrainedStiffness(omega_f)
Fc_global = s2.GlobalConstrainedForce(omega_f)

# %%% Solve for the free DOFs

u_free = s2.SolveUfree(Kc_global, Fc_global)

# %%% Solve for the support reactions

f_supp = s2.SupportReactions(s2.GlobalStiffness(omega_f), u_free, s2.GlobalForce(omega_f))


# %%%%
print(f'Solution of u_free = \n{u_free}\n')

# %%%%
print(f'Global constrained stiffness matrix = \n{Kc_global}\n')

# %%%%

print(f'Global constrained force vector = \n{Fc_global}\n')

# %%%%

print(f'Global support reactions = \n{f_supp}\n')

# %%%
u_elem = s2.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')

# %%% get element displacements

disp = s2.ElementDisplacements(u_elem, omega_f)
# forces = s2.ElementForces(u_elem, omega)
s2.PlotElementDisplacements(disp,scale=1)
# disp_mid_beam = disp[elem1.id][1][-1]
# disp_mid_rod = disp[elem3.id][1][0]
# s2.PlotMoments(forces,scale=1)