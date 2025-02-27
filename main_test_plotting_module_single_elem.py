# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:23:06 2024

@author: rensv
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# import sys
# import os

# # Import the pydynsm module

# ## find the module in the parent directory (DELETE THESE LINES AFTER WHEN THE NOTEBOOK IS READY, DIRECTLY IMPORT THE MODULE)

# ### Get the parent directory
# parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
# ### Add the parent directory to the system path
# sys.path.append(parent_dir)

## Now you can import the entire pydynsm module
import pydynsm as PDM

Assembler = PDM.Assembler

# %% test 1: 1d rod with UDL axial loading

# %%% Initialise an assembler with your project name

Assembler = PDM.Assembler  #Import the Assembler class from the PDM module

s1 = Assembler('1D Rod',analysis_type= 'new')

# %%% Parameters
E = 210e9
A = 1/300
rho = 3e6
F_0 = 1e06 + 0j 
L  = 1
Omega = 1
ksi = 0.01
EI = 1.5 * 7e07
I = EI/E

# %%% Create nodes from the assembler

node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,0)
# node3 = s1.CreateNode(L,0)

# %%% print the coords of 2 nodes here
print(f'The coordinate of node 1 is {node1.get_coords()}')
print(f'The coordinate of node 2 is {node2.get_coords()}')

# %%% create elements
elem = s1.CreateElement([node1, node2])
# elem1 = s1.CreateElement([node2,node3])

# %%% plot structural elements
s1.PlotStructure(plot_elements=True)
node1.fix_node('x','z')
node2.fix_node('x','z')
# elem.SetSection('Rod', {'E': E, 'rho':rho, 'A':A})
elem.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': I})
q_r = lambda omega: 1e02 if omega == Omega else 0
q_b = lambda omega: -1e6 if omega == Omega else 0
elem.AddDistributedLoad(z=q_b)
#%% test 1: free end disp @ omega = 100 rad/s need being revisited in maple
#%%% add UDL
# elem.AddDistributedLoad(x=q_r,z=q_b)
# p_x = lambda omega: -1e2 if omega == Omega else 0
# p_x2 = lambda omega: 1e2 if omega == Omega else 0
# node2.add_load(x=p_x)

#%%% num_dof
dof_indices, num_dof = s1.get_dofs_elements()

# %%% run the connectivity
dof_indices, B = s1.get_B_matrix()

# %%%%

s1.run_connectivity()

# %%% Get the global stiffness and force matrices

K_global = s1.GlobalStiffness(Omega)
F_global = s1.GlobalForce(Omega)

# %%% Constrain matrix
Kc_global = s1.GlobalConstrainedStiffness(Omega)
Fc_global = s1.GlobalConstrainedForce(Omega)

# %%% solve for free end disp
u_free = s1.SolveUfree(Kc_global,Fc_global)
u_elem =  s1.FullDisplacement(u_free)
disp = s1.ElementDisplacements(u_elem, Omega)
s1.PlotElementDisplacements(disp,scale=1000)
# s1.PlotAxialforces(force,scale=1e5)
# s1.PlotMoments(force,scale=1e-6)
# s1.PlotShearforces(force,scale = 1e-7)
# s1.PlotAxialstresses(force,scale = 100)

# # %% test 2: 1d eb beam with UDL
# s2 = Assembler('1D EB',analysis_type= 'new')
# node3 = s2.CreateNode(0,0)
# node4 = s2.CreateNode(L,0)
# elem1 = s2.CreateElement([node3, node4])
# node3.fix_node('x','z')
# node4.fix_node('x','z')
# s2.PlotStructure(plot_elements=True)
# elem1.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': I})
# s2.run_connectivity()
# q_b = lambda omega: 1e6 if omega == Omega else 0
# elem1.AddDistributedLoad(z=q_b)
# K_global = s2.GlobalStiffness(Omega)
# F_global = s2.GlobalForce(Omega)
# Kc_global = s2.GlobalConstrainedStiffness(Omega)
# Fc_global = s2.GlobalConstrainedForce(Omega)
# u_free = s2.SolveUfree(Kc_global, Fc_global)
# f_supp = s2.SupportReactions(s2.GlobalStiffness(Omega), u_free, s2.GlobalForce(Omega))
# u_elem = s2.FullDisplacement(u_free)
# disp = s2.ElementDisplacements(u_elem, Omega)
# force = s2.ElementForces(u_elem, Omega)
# s2.PlotElementDisplacements(disp,scale=1000)
# s2.PlotShearforces(force,scale=1e-6)
# s2.PlotMoments(force,scale=1e-6)

# # %% test 3: 1d eb cantilever with UDL
# s3 = Assembler('1D EB cantilever',analysis_type= 'new')
# node5 = s3.CreateNode(0,0)
# node6 = s3.CreateNode(L,0)
# elem2 = s3.CreateElement([node5, node6])
# node5.fix_node('x','z','phi_y')
# node6.fix_node('x')
# s3.PlotStructure(plot_elements=True)
# elem2.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': I})
# s3.run_connectivity()
# q_b = lambda omega: 1e6 if omega == Omega else 0
# elem2.AddDistributedLoad(z=q_b)
# K_global = s3.GlobalStiffness(Omega)
# F_global = s3.GlobalForce(Omega)
# Kc_global = s3.GlobalConstrainedStiffness(Omega)
# Fc_global = s3.GlobalConstrainedForce(Omega)
# u_free = s3.SolveUfree(Kc_global, Fc_global)
# f_supp = s3.SupportReactions(s3.GlobalStiffness(Omega), u_free, s3.GlobalForce(Omega))
# u_elem = s3.FullDisplacement(u_free)
# disp = s3.ElementDisplacements(u_elem, Omega)
# force = s3.ElementForces(u_elem, Omega)
# s3.PlotElementDisplacements(disp,scale=100)
# s3.PlotShearforces(force,scale=1e-7)
# s3.PlotMoments(force,scale=1e-7)

# # %% test 3: 1d eb cantilever with single point load
# s4 = Assembler('1D EB cantilever',analysis_type= 'new')
# node7 = s4.CreateNode(0,0)
# node8 = s4.CreateNode(L,0)
# elem3 = s4.CreateElement([node7, node8])
# node7.fix_node('x','z','phi_y')
# node8.fix_node('x')
# s4.PlotStructure(plot_elements=True)
# elem3.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': I})
# s4.run_connectivity()
# q_b = lambda omega: 1e6 if omega == Omega else 0
# p_z = lambda omega: -1e5 if omega == Omega else 0
# node8.add_load(z=p_z)
# K_global = s4.GlobalStiffness(Omega)
# F_global = s4.GlobalForce(Omega)
# Kc_global = s4.GlobalConstrainedStiffness(Omega)
# Fc_global = s4.GlobalConstrainedForce(Omega)
# u_free = s4.SolveUfree(Kc_global, Fc_global)
# f_supp = s4.SupportReactions(s4.GlobalStiffness(Omega), u_free, s4.GlobalForce(Omega))
# u_elem = s4.FullDisplacement(u_free)
# disp = s4.ElementDisplacements(u_elem, Omega)
# force = s4.ElementForces(u_elem, Omega)
# s4.PlotElementDisplacements(disp,scale=1000)
# s4.PlotShearforces(force,scale=1e-6)
# s4.PlotMoments(force,scale=1e-6)