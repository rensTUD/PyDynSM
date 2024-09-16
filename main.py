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

Node = PDM.Node
Element = PDM.Element
Assembler = PDM.Assembler

# %% Example

# %%% Initialise an assembler with your project name

s1 = Assembler('beam')

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
node2 = s1.CreateNode(L,0)

# %%% Now set the constraints directly onto the nodes
node1.FixDofs('x', 'z')
node2.FixDofs('z')

# %%%% test adding a load to the node

# define a lambda function running over omega for p_x
# TODO - see how to add a time domain load instead of freq domain load
omega_f = 100
p_x = lambda omega: 1 if omega == omega_f else 0

# add a load directly to node2

node1.AddLoad([0, 0, 0])

# %%%% Plot nodes

# s1.PlotStructure()

# %%% Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])

# %%% try to set a section that is not implemented yet
elem.SetSection('Rod2', {'EA': EA, 'rhoA':rhoA})
 
# %%% now set the sections (or element types for that matter)
elem.SetSection('Rod', {'EA': EA, 'rhoA':rhoA})
elem.SetSection('EB Beam', {'EI': EI, 'rhoA':rhoA})

# %%%% plot elements too

s1.PlotStructure(plot_elements=True)


# %%% Add distributed load
# TODO - should check whether this is correct, for now it just adds the load as a function
q_r = lambda omega: 0 if omega == omega_f else 0
q_b = lambda omega: 1e06 if omega == omega_f else 0
elem.AddDistributedLoad([q_r,q_b,0])

# %%% Get the global stiffness and force matrices

# global_K = s1.GlobalStiffness(omega)
# global_F = s1.GlobalForce(omega)

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

