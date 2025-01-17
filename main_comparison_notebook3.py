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
E = 210e9
EA = 7e6
A = EA/E
EI = 1.5 * 7e06 
I = EI/E
rhoA = 1e03 
rho = rhoA/A
q_r = 1*1e02 + 0j 
q_b = 1*1e06 + 0j 
L  = 1
omega = 100  
ksi = 0.01 
omega_f = omega


# %%% Create nodes from the assembler


node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,0)

# %%% Now set the constraints directly onto the nodes
node1.fix_node('x', 'z')
node2.fix_node('z')


# %%% Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])


 
# %%% now set the sections (or element types for that matter)
# elem.SetSection('Rod', {'EA': EA, 'rhoA':rhoA})
# elem.SetSection('EB Beam', {'EI': EI, 'rhoA':rhoA})

elem.SetSection('Rod', {'E': E, 'A':A, 'rho':rho})
elem.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I})


# %%%% run connectivity

s1.run_connectivity()

# %%%% plot elements too

# s1.PlotStructure(plot_elements=True)


# %%% Add distributed load per DOF in global coord system

q_r = lambda omega: 1e2 if omega == omega_f else 0
q_b = lambda omega: 1e6 if omega == omega_f else 0
elem.AddDistributedLoad(x=q_r, z=q_b)



# %%% Get the global stiffness and force matrices

K_lgobal = s1.GlobalStiffness(omega)
F_global = s1.GlobalForce(omega)

# %%% get the global constrained stiffness and force

Kc_global = s1.GlobalConstrainedStiffness(omega)
Fc_global = s1.GlobalConstrainedForce(omega)

# %%% solve as well for u_free
u_free = s1.SolveUfree(Kc_global, Fc_global)

# %%% and for the support reactions

f_supp = s1.SupportReactions(s1.GlobalStiffness(omega), u_free, s1.GlobalForce(omega))

# %%%%
# TODO: IMPORTANT: RESULTS ARE SLIGHTLY OFF FROM THE NOTEBOOK 3.3 (SHOULD INVESTIGATE)

# %%%%
print(f'Solution of u_free = \n{u_free}\n')

# %%%%
print(f'Global constrained stiffness matrix = \n{Kc_global}\n')

# %%%%

print(f'Global constrained force vector = \n{Fc_global}\n')

# %%%%

print(f'Global support reactions = \n{f_supp}\n')


# %%%

u_elem = s1.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')

# %%% get element displacements

disp = s1.ElementDisplacements(u_elem, omega)

# %%% Plot displacements

s1.PlotElementDisplacements(disp,scale=100.0)

