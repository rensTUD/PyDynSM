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
q_r = 0*1e02 + 0j 
q_b = 1*1e06 + 0j 
L  = 1
omega = 100  
ksi = 0.01 
omega_f = omega



# %%% Create nodes from the assembler

node1 = s1.CreateNode(0,0)
# node2 will have no 'z' displacement - will be handled without setting stiffness to infinity (results are still to be verified)
node2 = s1.CreateNode(L,L,dof_config = [['x'],['phi_y']]) 
node3 = s1.CreateNode(0,L)

# %%% Now set the constraints directly onto the nodes

node1.fix_node('x', 'z')
node3.fix_node('x','z', 'phi_y')

# %%%% test adding a load to the node

# define a lambda function running over omega for p_x - if constant is given that value will be evaluated for all frequencies
omega_f = 100
p_x = lambda omega: 10 if omega == omega_f else 0

# add a load directly to node2
node2.add_load(x=p_x)

# %%%% Plot nodes

s1.PlotStructure()

# %%% Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])
elem1 = s1.CreateElement([node2, node3])

# plot elements too
s1.PlotStructure(plot_elements=True)
 
# %%% now set the sections onto the elements

elem.SetSection('Rod', {'E': E, 'A':A, 'rho':rho})
elem.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I})

elem1.SetSection('Rod', {'E': E, 'A':A, 'rho':rho})
elem1.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I})

# %%% Run connectivity

s1.run_connectivity()

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

# %%% Solve for the free DOFs

u_free = s1.SolveUfree(Kc_global, Fc_global)

# %%% Solve for the support reactions

f_supp = s1.SupportReactions(s1.GlobalStiffness(omega), u_free, s1.GlobalForce(omega))


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

