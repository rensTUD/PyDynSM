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

s1 = Assembler('beam',analysis_type='new')

# %%% Parameters
EA = 7e6
EI = 1.5 * 7e06 
rhoA = 1e03 
q_r = 0*1e02 + 0j 
q_b = 1*1e06 + 0j 
L  = 1
omega = 100  
ksi = 0.01 



# %%% Create nodes from the assembler


node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,L,dof_config = [['x'],['phi_y']])
node3 = s1.CreateNode(2*L,L)
node4 = s1.CreateNode(L,2*L)

# %%% Now set the constraints directly onto the nodes
node1.fix_node('x', 'z')
node2.fix_node('x')
node3.fix_node('x','z', 'phi_y')
node4.fix_node('x','z')
# %%%% test adding a load to the node

# define a lambda function running over omega for p_x
# TODO - see how to add a time domain load instead of freq domain load
omega_f = 100
p_x = lambda omega: 1 if omega == omega_f else 0

# add a load directly to node2

node1.AddLoad([p_x, 0, 0])

# %%%% Plot nodes

# s1.PlotStructure()

# %%% Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])
elem2 = s1.CreateElement([node2, node3])
elem3 = s1.CreateElement([node2, node4])
# %%% try to set a section that is not implemented yet
elem.SetSection('Rod2', {'EA': EA, 'rhoA':rhoA})
 
# %%% now set the sections (or element types for that matter)
elem.SetSection('Rod', {'EA': EA, 'rhoA':rhoA})
elem.SetSection('EB Beam', {'EI': EI, 'rhoA':rhoA})

elem2.SetSection('Rod', {'EA': EA, 'rhoA':rhoA})
elem2.SetSection('EB Beam', {'EI': EI, 'rhoA':rhoA})

elem3.SetSection('Rod', {'EA': EA, 'rhoA':rhoA})
elem3.SetSection('EB Beam', {'EI': EI, 'rhoA':rhoA})
# %%% test adding dof to an element that is not supported by the local dofs of that element

elem.fix_dof(node2, 'z')


# %%%%
elem.free_dof(node2, 'z')

# %%%% 

elem3.free_dof(node2,'x')
elem3.fix_dof(node2,'z')

# %% Testing for connectivity, B, and L matrices

# %%% get global dofs of the elements
dofs_indices = s1.get_dofs_elements()

# %%% run the connectivity
B = s1.get_B_matrix()

# %%%%

s1.run_connectivity()

# %%%% plot elements too

# s1.PlotStructure(plot_elements=True)

# %% get global stiffness matrix

K_gobal = s1.GlobalStiffness(omega)

# %%% Add distributed load
# TODO - should check whether this is correct, for now it just adds the load as a function
q_r = lambda omega: 1 if omega == omega_f else 0
q_b = lambda omega: 1 if omega == omega_f else 0
elem.AddDistributedLoad(x=q_r, z=q_b)

# %% get global force vector

F_global = s1.GlobalForce(omega)


# %%% Get the global stiffness and force matrices

global_K = s1.GlobalStiffness(omega)
global_F = s1.GlobalForce(omega)

# %%% get the global constrained stiffness and force

global_kc = s1.GlobalConstrainedStiffness(omega)
global_fc = s1.GlobalConstrainedForce(omega)

# %%% solve as well for u_free
u_free = s1.SolveUfree(global_kc, global_fc)

# %%% and for the support reactions

f_supp = s1.SupportReactions(s1.GlobalStiffness(omega), u_free, s1.GlobalForce(omega))

# %%%%
# TODO: IMPORTANT: RESULTS ARE SLIGHTLY OFF FROM THE NOTEBOOK 3.3 (SHOULD INVESTIGATE)

# %%%%
print(f'Solution of u_free = \n{u_free}\n')

# %%%%
print(f'Global constrained stiffness matrix = \n{global_kc}\n')

# %%%%

print(f'Global constrained force vector = \n{global_fc}\n')

# %%%%

print(f'Global support reactions = \n{f_supp}\n')

# %%% check global dofs 

elem.GlobalDofs()

# %%%

u_elem = s1.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')

