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
Constrainer = PDM.Constrainer
Element = PDM.Element
Assembler = PDM.Assembler

# %% Example

# %%% Initialise an assembler with your project name

s1 = Assembler('beam')

# %%% Parameters
EA = 7e6
EI = 1.5 * 7e06 # MODIFIED
rhoA = 1e03  # MODIFIED
q_r = 0*1e02 + 0j # MODIFIED
q_b = 1*1e06 + 0j # MODIFIED
L  = 1
omega = 100  # MODIFIED
ksi = 0.01 # MODIFIED


Node.clear()
Element.clear()

# %%% Create nodes

node1 = Node(0,0,s1)
node2 = Node(L,0,s1)

# %%% Now set it directly on the nodes
node1.fix_dofs('x', 'z')
node2.fix_dofs('x', 'z')


# %%%% Plot nodes

# s1.PlotStructure()

# %%% Create element

# initialise element
elem = Element ( [node1, node2] , s1)
 
# %%% trial setting elements
elem.SetSection('Rod', {'EA': EA, 'rhoA':rhoA})
elem.SetSection('EB Beam', {'EI': EI, 'rhoA':rhoA})

# %%%% plot elements too

# s1.PlotStructure(plot_elements=True)

# %%% set boundary conditions



# %%% Add distributed load
elem.AddDistributedLoad([q_r,q_b,0], omega)


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

