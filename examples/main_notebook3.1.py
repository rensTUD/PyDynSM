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

s1 = Assembler('beam')

# %%% Parameters
E = 210e9
EA = 7e6
A = EA/E
EI = 1.5 * 7e06
I = EI/E
W = I
rhoA = 1e03 
rho = rhoA/A
q_r = 0*1e02 + 0j 
q_b = 1*1e06 + 0j 
L  = 1
ksi = 0.01
omega_f = 100
F = 1e06

# %%% Create nodes from the assembler
node1 = s1.CreateNode(0,0)
# node2 will have no 'z' displacement - will be handled without setting stiffness to infinity (results are still to be verified)
node2 = s1.CreateNode(L,0) 
node1.fix_node('x','z','phi_y')
node2.fix_node('z','phi_y')
elem = s1.CreateElement([node1, node2])

# %%% plot structural elements
s1.PlotStructure(plot_elements=True)

elem.SetSection('Rod', {'E': E, 'A':A, 'rho':rho, 'ksi':ksi})
F_0 = lambda omega: F if omega == omega_f else 0
node2.add_load(x=F_0)
s1.run_connectivity()

# %%% Get the global stiffness and force matrices

K_global = s1.GlobalStiffness(omega_f)
F_global = s1.GlobalForce(omega_f)

# %%% get the global constrained stiffness and force

Kc_global = s1.GlobalConstrainedStiffness(omega_f)
Fc_global = s1.GlobalConstrainedForce(omega_f)

# %%% Solve for the free DOFs

u_free = s1.SolveUfree(Kc_global, Fc_global)

# %%% Solve for the support reactions

f_supp = s1.SupportReactions(s1.GlobalStiffness(omega_f), u_free, s1.GlobalForce(omega_f))


# %%%%
print(f'Solution of u_free = \n{u_free}\n')
