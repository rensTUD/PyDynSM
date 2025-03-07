# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 18:24:14 2025

@author: thana
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# %% add package to sys and import

import pydynsm as PDM

Assembler = PDM.Assembler

s1 = Assembler('test')
# %% inputs
l1 = 8
h1 = 4

# material props
Es = 210e9  # Elasticity modulus steel
Ec = 3.1e10  # Elasticity modulus concrete
rhoc = 2500  # density concrete
rhos = 7800  # density steel
Ab = 2  # beam cross sectional area
Ib = 0.00395  # beam 2nd mmt of intertial
rc = 9.5 * 0.001
Ac = np.pi * rc**2
Ic = 0.25 * np.pi * rc**4
omega_f = 100

# %% Create Nodes
'''
To create nodes we call the CreateNode method. The input is the x,z location 
of the node. 

Default, the assigned dofs will be 'x', 'z', and 'phi_y'. However, you can choose
which dofs to assign. Note that similar behaviour is obtained when the dof you leave out is 
fixed. The difference is that in the latter case, support forces are generated.

In this example, we limit node 3 to having vertical (z) and rotational (phi_y) motion,
as we will only assign an euler bernoulli beam to this node, which has no stiffness in x.
Node 4 will be assigned a rod, which will be placed diagonally, and thus only having stiffness
in horizontal (x) and vertical (z) motion.
'''
node0 = s1.CreateNode(0, 0)
node1 = s1.CreateNode(l1, 0)
node2 = s1.CreateNode(2*l1, 0, dof_config=['z', 'phi_y'])
node3 = s1.CreateNode(0, h1, dof_config=['x', 'z'])

# %% Plotting nodes
'''
Once the nodes are created, you can perform a visual check by calling the PlotStructure method.
'''
s1.PlotStructure()

# %% Elements
'''
To create elements, we use the CreateElement method, which takes two nodes as input.

In this example, we will 3 elements.
'''
# beam
beams = []

# Define node pairs for beam elements
beam_node_pairs = [
    (node0, node1),
    (node1, node2),
    (node3, node1)]

# Loop to create and append beam elements
for pair in beam_node_pairs:
    beam = s1.CreateElement(list(pair))
    beams.append(beam)  # Append beam element to list
    
# %%
'''
After creating the elements, we can inspect our structure again. However we 
use the kwarg 'plot_elements=True' to enforce this.
'''

s1.PlotStructure(plot_elements=True)

# %% Constraint nodes, decouple dofs
'''
Regarding the dofs, we have multiple options. Once a node is initialised, all assigned 
dofs will be set to 'free'. To set it to fixed, you use the method 'node.fix_node(dof1, dof2, ...)'

In this example, we apply the following constraints:
    node0: 'x', 'z'
    node1: 'phi_y'
    node2: 'z'
    node3: 'x', 'z'
    
Important to note is that nodal dofs are saved separately from element dofs. However, 
once an element has been assigned, its dofs will be 'monolithically' connected to the node.
You can relax this condition, and by doing so create 'hinges'. 

In this example, we will create two hinges:
    element0: at node2
    element1: at node2

once decoupled, you can also change the specific dof of the elements (which we won't do in this example).
To do so, you can the methods: 
    fix dofs: element.fix_dof(dof1, dof2, ...)                                                                      
    free dofs: element.free_dof(dof1, dof2, ...) 

NOTE: As will be shown later, beams[2] will be the only element coupled to node1
whilst still having the dof 'phi_y'. In order to avoid a singular matrix, we have
to constrain this dof as well. Think why!                                                                  

'''
node0.fix_node('x','z')
node1.fix_node('phi_y')
node2.fix_node('z')
node3.fix_node('x','z')


beams[0].decouple_dof(node1, 'phi_y')
beams[1].decouple_dof(node1, 'phi_y')


# %% assign sections
'''
Currently we have only created elements, to give them 'stiffness', we have to add sections.
There are various sections available. In this example we will use:
    Euler bernoulli beam
    Euler bernoulli beam on foundation
    Rod

To assign sections, we use the element.SetSection(element_type, parameters), where:
    element_type: string with the element type:
                                        - 'EulerBernoulli Beam'
                                        - 'EulerBernoulli Beam with foundation'
                                        - 'Rod'

'''
# set first element to EB - rod
beams[0].SetSection('EulerBernoulli Beam', {'E': Ec, 'A': Ab, 'rho': rhoc, 'Ib': Ib, 'Wb': Ib})
beams[0].SetSection('Rod', {'E': Ec, 'A': Ac, 'rho': rhoc})

# set second element to EB
beams[1].SetSection('EulerBernoulli Beam', {'E': Ec, 'A': Ab, 'rho': rhoc, 'Ib': Ib, 'Wb': Ib})

# set third element to rod
beams[2].SetSection('Rod', {'E': Ec, 'A': Ac, 'rho': rhoc})


# %% loading definition
'''
Currently we can add nodal loads and distributed loads. 

In case of nodal loads, these are in the GLOBAL coordinates. To add a load, we
use the method:
    node.add_load(dof1=f1(omega), dof2=f2(omega), ...)
    
Here, f1(omega) can be any function depending on omega. Furthermore, if you define it 
as a CONSTANT, it means that for every frequency that load will be present.

In case of the distributed load, we add loads in the LOCAL coordinates, To add a load, we
use the method:
    element.AddDistributedLoad(dof1=f1(omega), dof2=f2(omega), ...)

Where the exact same logic applies regarding the forcing functions. 
'''

# define a lambda function running over omega for p_z, which return 1e6 if omega = omega_p
omega_p = 50
p_z = lambda omega: 1e6 if omega == omega_p else 0

# add a load directly to node1
node1.add_load(z=p_z)

# add distributed loads on the horizontal beams
beams[0].AddDistributedLoad(z=1e03)
beams[1].AddDistributedLoad(z=1e03)


# %% connectivity
'''
Once you have finalised your system, you have to run a command that 'glues' the
structure together. This is done by the method:
    - run_connectivity()
    
You have to do this only once!
'''

s1.run_connectivity()

# %% Get the global stiffness and force matrices
'''
To obtain the stiffness matrix and force vector, we run the following two methods:
    K = GlobalConstrainedStiffness(omega)
    F = GlobalConstrainedForce(omega)
    
These will return the constrained versions, which means any constrains (e.g. zero displacements)
are applied. 

NOTE: for every frequency you have to run these commands!

If you wish to obtain the unconstrained versions, you can use the following methods:
    k = GlobalStiffness(omega)
    f = GlobalForce(omega)
'''
K_global = s1.GlobalStiffness(omega_f)
F_global = s1.GlobalForce(omega_f)
Kc_global = s1.GlobalConstrainedStiffness(omega_f)
Fc_global = s1.GlobalConstrainedForce(omega_f)


# %% Solving
'''
Now we are ready to solve for the free displacements. To calculate these the following method is applied:
    SolveUfree(K,F)
    
Where:
    K = global constrained stiffness matrix
    F = global constrained Force vector
    
To obtain the support forces, the following method is used:
    SupportReactions(k,u_free,f)
    
Where:
    K = global stiffness matrix
    F = global Force vector
    
'''
# solve free displacements
u_free = s1.SolveUfree(Kc_global, Fc_global)

# solve support forces
f_supp = s1.SupportReactions(s1.GlobalStiffness(omega_f), u_free, s1.GlobalForce(omega_f))

# %% post-processing displacements
'''
In terms of post processing we have two choices:
    - Plot GLOBAL
    - Plot LOCAL

When we plot GLOBAL, the values for the whole structure will be plotted. 
In case of LOCAL plots, the values of the chosen element will be plotted

As an example, let us plot the global displacements of the whole structure, to 
do so, we run the following commands:
    u_elem = s1.FullDisplacement(u_free)
    diplacements = s1.ElementDisplacements(u_elem, omega, num_points=100)
    s1.PlotElementDisplacements(displacements, scale=1e5)
    
where:
    u_elem = the 'full' displacement vector, including the constrained ones
    num_points = the amount of plotting points.
    scale = the value with which to scale the displacements to make them more visible 
    
NOTE: plotted displacements may not be shown realistically in case an element does NOT
      contain a stiffness in one of the linear dofs.
      
When plotting forces, we first have to calculate the forces:
    forces = s1.ElementForces(u_elem,omega,num_points)

After which we can choose for:
    s1.PlotMoments(forces, scale)
    s1.PlotShearforces(forces, scale)
    s1.PlotAxialforces(forces, scale)
    
To plot LOCAL forces we call either of the following methods:
    s1.StructurePlotter.plot_element_moment(element, forces)
    s1.StructurePlotter.plot_element_axial_force(element, forces)
    s1.StructurePlotter.plot_element_shear_force(element, forces)
    
We can do the same for stresses:
    stresses = s1.ElementStresses(u_elem, omega, num_points)
    s1.StructurePlotter.plot_element_bending_stress(element, stresses)
    s1.StructurePlotter.plot_element_axial_stress(element, stresses)
    s1.StructurePlotter.plot_element_shear_stress(element, stresses)
    
'''
# plot global displacements
u_elem = s1.FullDisplacement(u_free)
displacements = s1.ElementDisplacements(u_elem, omega_f, num_points=100)
s1.PlotElementDisplacements(displacements, scale=1e05)

# plot GLOBAL forces
forces = s1.ElementForces(u_elem, omega_f, num_points=200)
s1.PlotMoments(forces, scale=1e-3)
s1.PlotShearforces(forces, scale=1e-3)
s1.PlotAxialforces(forces, scale=1)

# plot LOCAL forces
s1.StructurePlotter.plot_element_moment(beams[0], forces)
s1.StructurePlotter.plot_element_axial_force(beams[0], forces)
s1.StructurePlotter.plot_element_shear_force(beams[0], forces)

# plot LOCAL stresses
stresses = s1.ElementStresses(u_elem, omega_f, num_points=200)
s1.StructurePlotter.plot_element_bending_stress(beams[0], stresses)
s1.StructurePlotter.plot_element_axial_stress(beams[0], stresses)
s1.StructurePlotter.plot_element_shear_stress(beams[0], stresses)

