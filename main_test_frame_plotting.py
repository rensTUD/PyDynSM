# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:33:47 2025

@author: GREY
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# %% add package to sys and import 

import pydynsm as PDM


Assembler = PDM.Assembler

# %% test portal frame

# %%% Initialise an assembler with your project name

s1 = Assembler('beam',analysis_type='new')

# %%% Parameters
E = 210e9
EA = 7e7
A = EA/E
EI = 1.5 * 7e07
I = EI/E
W = I
rhoA = 1e03 
rho = rhoA/A
q_r = 0*1e02 + 0j 
q_b = 1*1e06 + 0j 
L  = 1
omega = 0.001
ksi = 0.01
ksi = 0
omega_f = omega

# %%% Create nodes from the assembler
node1 = s1.CreateNode(0,0)
# node2 will have no 'z' displacement - will be handled without setting stiffness to infinity (results are still to be verified)
node2 = s1.CreateNode(0,L) 
node3 = s1.CreateNode(L,L)
node4 = s1.CreateNode(L,0)

node1.fix_node('x','z','phi_y')
node4.fix_node('x','z','phi_y')

elem = s1.CreateElement([node1, node2])
elem1 = s1.CreateElement([node2, node3])
elem2= s1.CreateElement([node3, node4])

# %%% plot structural elements
s1.PlotStructure(plot_elements=True)

# elem.SetSection('Rod', {'E': E, 'A':A, 'rho':rho})
elem.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': W})

elem1.SetSection('Rod', {'E': 1e10*E, 'A':A, 'rho':rho})
elem1.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':1e15*I, 'Wb': W})

elem2.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': W})


p_x = lambda omega: 1e6 if omega == omega_f else 0
node3.add_load(x=p_x)
s1.run_connectivity()

# %%% Get the global stiffness and force matrices

K_global = s1.GlobalStiffness(omega)
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
force = s1.ElementForces(u_elem, omega,num_points=201)
stress = s1.ElementStresses(u_elem, omega,num_points=201)

# %%% Plot displacements

s1.PlotElementDisplacements(disp,scale=100)
s1.PlotMoments(force,scale = 1e-12)
s1.PlotAxialforces(force,scale = 1e-7)
s1.PlotShearforces(force,scale = 1e-8)
