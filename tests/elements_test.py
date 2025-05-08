# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:23:06 2024

@author: rensv
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# import package
import pydynsm as PDM


# %% start assembler
Assembler = PDM.Assembler

# %% Set analysis frequency

omega = 100 

# %% EB Beam

# Initialise an assembler with your project name

s1 = Assembler('EB beam')

# list element types
s1.ListElementTypes()

# Parameters

# material parameters
rho = 1
E = 1

# section parameters
A = 1
I = 1

# Length
L  = 1

# damping
ksi = 0.0 

# Nodes

node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,0)

# Set the constraints of the nodes
node1.fix_node('x','z')
node2.fix_node('x','z')

# Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])

# plot the nodes and elements
# s1.PlotStructure(plot_elements = True)

# set the element type
elem.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': 0})

# add distributed unit load 
q_z = 1
elem.AddDistributedLoad(z=q_z)

# Solve 

# run connectivity
s1.run_connectivity()


# get the global constrained stiffness and force
Kc_global = s1.GlobalConstrainedStiffness(omega)
print(f'Global constrained stiffness matrix = \n{Kc_global}\n')

Fc_global = s1.GlobalConstrainedForce(omega)
print(f'Global constrained force vector = \n{Fc_global}\n')

# solve for u_free
u_free = s1.SolveUfree(Kc_global, Fc_global)
print(f'Solution of u_free = \n{u_free}\n')

# and for the support reactions
f_supp = s1.SupportReactions(s1.GlobalStiffness(omega), u_free, s1.GlobalForce(omega))
print(f'Global support reactions = \n{f_supp}\n')

# Get full displacements
u_elem = s1.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')

# get element displacements
disp = s1.ElementDisplacements(u_elem, omega, num_points=100)

# Plot displacements
s1.PlotElementDisplacements(disp,scale=100.0)

# %% EB Beam with foundation 

# Initialise an assembler with your project name

s1 = Assembler('EB beam foundation')

# Parameters

# foundation stiffness
kd = 1000


# Nodes

node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,0)

# Set the constraints of the nodes
node1.fix_node('x','z')
node2.fix_node('x','z')

# Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])

# plot the nodes and elements
# s1.PlotStructure(plot_elements = True)

# set the element type
elem.SetSection(
    'EulerBernoulli Beam with foundation', 
    {'E': E, 
     'A':A, 
     'rho':rho, 
     'Ib':I, 
     'Wb': 0, 
     'kd': kd,
     }
    )

# add distributed unit load 
q_z = 1
elem.AddDistributedLoad(z=q_z)

# Solve 

# run connectivity
s1.run_connectivity()


# get the global constrained stiffness and force
Kc_global = s1.GlobalConstrainedStiffness(omega)
print(f'Global constrained stiffness matrix = \n{Kc_global}\n')

Fc_global = s1.GlobalConstrainedForce(omega)
print(f'Global constrained force vector = \n{Fc_global}\n')

# solve for u_free
u_free = s1.SolveUfree(Kc_global, Fc_global)
print(f'Solution of u_free = \n{u_free}\n')

# and for the support reactions
f_supp = s1.SupportReactions(s1.GlobalStiffness(omega), u_free, s1.GlobalForce(omega))
print(f'Global support reactions = \n{f_supp}\n')

# Get full displacements
u_elem = s1.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')

# get element displacements
disp = s1.ElementDisplacements(u_elem, omega, num_points=100)

# Plot displacements
s1.PlotElementDisplacements(disp,scale=100.0)

# %% EB Beam with foundation and no end attachments

# Initialise an assembler with your project name

s1 = Assembler('EB beam foundation end attachments')

# Parameters

# foundation stiffness
kd = 1000


# Nodes

node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,0)

# Set the constraints of the nodes
node1.fix_node('x','z')
node2.fix_node('x','z')

# Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])

# plot the nodes and elements
# s1.PlotStructure(plot_elements = True)

# set the element type
elem.SetSection(
    'EulerBernoulli Beam foundation attachments', 
    {'E': E, 
     'A':A, 
     'rho':rho, 
     'Ib':I, 
     'Wb': 0, 
     'kd': kd,
     }
    )

# add distributed unit load 
q_z = 1
elem.AddDistributedLoad(z=q_z)

# Solve 

# run connectivity
s1.run_connectivity()


# get the global constrained stiffness and force
Kc_global = s1.GlobalConstrainedStiffness(omega)
print(f'Global constrained stiffness matrix = \n{Kc_global}\n')

Fc_global = s1.GlobalConstrainedForce(omega)
print(f'Global constrained force vector = \n{Fc_global}\n')

# solve for u_free
u_free = s1.SolveUfree(Kc_global, Fc_global)
print(f'Solution of u_free = \n{u_free}\n')

# and for the support reactions
f_supp = s1.SupportReactions(s1.GlobalStiffness(omega), u_free, s1.GlobalForce(omega))
print(f'Global support reactions = \n{f_supp}\n')

# Get full displacements
u_elem = s1.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')

# get element displacements
disp = s1.ElementDisplacements(u_elem, omega, num_points=100)

# Plot displacements
s1.PlotElementDisplacements(disp,scale=100.0)

# %% EB Beam with foundation and end attachments

# Initialise an assembler with your project name

s1 = Assembler('EB beam foundation end attachments')

# Parameters

# foundation stiffness
kd = 0
Pm2 = 1
J2 = 1

# Nodes

node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,0)

# Set the constraints of the nodes
node1.fix_node('x','z')
node2.fix_node('x')

# Create element

# initialise element by setting its nodes and calling it from the assembler
elem = s1.CreateElement([node1, node2])

# plot the nodes and elements
# s1.PlotStructure(plot_elements = True)

# set the element type
elem.SetSection(
    'EulerBernoulli Beam foundation attachments', 
    {'E': E, 
     'A':A, 
     'rho':rho, 
     'Ib':I, 
     'Wb': 0, 
     'kd': kd,
     'Pm2': Pm2,
     'J2': J2
     }
    )

# add distributed unit load 
q_z = 1
elem.AddDistributedLoad(z=q_z)

# Solve 

# run connectivity
s1.run_connectivity()


# get the global constrained stiffness and force
Kc_global = s1.GlobalConstrainedStiffness(omega)
print(f'Global constrained stiffness matrix = \n{Kc_global}\n')

Fc_global = s1.GlobalConstrainedForce(omega)
print(f'Global constrained force vector = \n{Fc_global}\n')

# solve for u_free
u_free = s1.SolveUfree(Kc_global, Fc_global)
print(f'Solution of u_free = \n{u_free}\n')

# and for the support reactions
f_supp = s1.SupportReactions(s1.GlobalStiffness(omega), u_free, s1.GlobalForce(omega))
print(f'Global support reactions = \n{f_supp}\n')

# Get full displacements
u_elem = s1.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')

# get element displacements
disp = s1.ElementDisplacements(u_elem, omega, num_points=100)

# Plot displacements
s1.PlotElementDisplacements(disp,scale=100.0)

