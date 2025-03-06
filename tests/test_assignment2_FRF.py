# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:59:03 2025

@author: GREY
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# %% add package to sys and import 

import pydynsm as PDM
import scipy.integrate as spi
Assembler = PDM.Assembler

s1 = Assembler('bridge',analysis_type='new')
# %% inputs
l1 = 4
l2 = 10
l3 = 6

# material props
E = 210e9  # Elasticity modulus steel
rho =  7800

# for beam
EI_b = 400e6
I_b = 0.0395
EA_b = 1.3e9
A_b = 2.2

# for cable
r_c = 9.5*0.01
A_c = np.pi*r_c**2
I_c = 0.25*np.pi*r_c**4

# for truss 
A_t = 0.1*A_b # truss erea
I_t = I_c

omega_w = 4
omega_c = 0.001

k_d = 1e2
# c_d = 1e6
c_d = 0
ksi = 0

# %% sturctural plot
# %%% nodes
node1 = s1.CreateNode(0,0)
node19 = s1.CreateNode(l2/5,0)
node20 = s1.CreateNode(l2/5+3,0)
node2 = s1.CreateNode(l2,0)
node3 = s1.CreateNode(2*l2,0)
node4 = s1.CreateNode(3*l2,0)
node5 = s1.CreateNode(4*l2,0)
node6 = s1.CreateNode(4*l2+l1,0)
node7 = s1.CreateNode(4*l2+l1,-l1)
node8 = s1.CreateNode(4*l2+l1+l3,0)
node9 = s1.CreateNode(4*l2+l1+2*l3,0)
node21 = s1.CreateNode(4*l2+l1+2*l3+l3/5,0)
node22 = s1.CreateNode(4*l2+l1+2*l3+l3/5+3,0)
node10 = s1.CreateNode(4*l2+l1+3*l3,0)
node11 = s1.CreateNode(4*l2+l1+4*l3,0)
node12 = s1.CreateNode(4*l2+l1,l1)
node13 = s1.CreateNode(4*l2,l1)
node14 = s1.CreateNode(4*l2+l1,2*l1)
node15 = s1.CreateNode(4*l2,2*l1)
node16 = s1.CreateNode(4*l2+l1,3*l1)
node17 = s1.CreateNode(4*l2,3*l1)
node18 = s1.CreateNode(4*l2,4*l1)
s1.PlotStructure()
# %%% elements
# beam
beams = []

# Define node pairs for beam elements
beam_node_pairs = [
    (node1, node19),
    (node19,node20),
    (node20,node2),
    (node2, node3),
    (node3, node4),
    (node4, node5),
    (node5, node6),
    (node6, node8),
    (node8, node9),
    (node9, node21),
    (node21, node22),
    (node22,node10),
    (node10, node11)
]

# Loop to create and append beam elements
for pair in beam_node_pairs:
    beam = s1.CreateElement(list(pair))
    beams.append(beam)  # Append beam element to list

# Truss
trusses = []

# Define node pairs for truss elements
truss_node_pairs = [
    (node5, node7),(node6, node7),
    (node5, node12), (node6, node12),
    (node12, node13), (node12, node14),
    (node12, node15), (node13, node15),
    (node13, node5), (node14, node15),
    (node14, node16), (node15, node16),
    (node15, node17), (node16, node17),
    (node16, node18), (node17, node18)
]

# Loop to create and append truss elements
for pair in truss_node_pairs:
    truss = s1.CreateElement(list(pair))
    trusses.append(truss)

# front cables
front_cables = []
front_cable_node_pairs = [(node18, node2), (node3, node17),
(node4, node15)]

for pair in front_cable_node_pairs:
    cable = s1.CreateElement(list(pair))
    front_cables.append(cable)
# back cables    
back_cables = []
back_cable_node_pairs = [(node9, node14),(node8, node12),(node10, node16)]
for pair in back_cable_node_pairs:
    cable = s1.CreateElement(list(pair))
    back_cables.append(cable)

# %% constraint nodes
node1.fix_node('z')
# node11.fix_node('z')
node7.fix_node('x','z')

# %% assign sections
for beam in beams[:-1]:
    beam.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':ksi})
    beam.SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':ksi})

beams[-1].SetSection('EulerBernoulli Beam with foundation', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':ksi, 'kd':k_d,'cd':c_d})
beams[-1].SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':ksi})
for cable in front_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':ksi})
    cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':ksi})
for cable in back_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':ksi})
    cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':ksi})
for truss in trusses:
    truss.SetSection('Rod', {'E': E, 'A':A_t, 'rho':rho,'ksi':ksi})
    truss.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_t, 'rho':rho, 'Ib':I_t,'ksi':ksi})
s1.PlotStructure(plot_elements=True)

# %% loading definition
# wind load. change the value later
def trap(a,b):

    '''
    a,b are start and end locations of the wind load, defined in the local direction.
    '''
    
    alpha = 0.4
    qw = 15  
    result, error = spi.quad(lambda x: np.exp(alpha * x), a, b)
    return result*qw

load1 = trap(0,0.5*l1)
print(load1)
load2 = trap(0.5*l1,1.5*l1)
print(load2)
load3 = trap(1.5*l1,2.5*l1)
print(load3)
load4 = trap(2.5*l1,3.5*l1)
print(load4)
load5 = trap(3.5*l1,4*l1)
print(load5)


x_coor = np.array([0,l1,2*l1,3*l1,4*l1])
exp_load = np.array([load1, load2, load3, load4, load5])
plt.figure()
plt.plot(exp_load/1e3,x_coor)
plt.ylabel('Height [m]')
plt.xlabel('Equivalent Wind load [kN]')
plt.grid()
plt.show()

w_x1 = lambda omega: load1 if omega == omega_w else 0
w_x2 = lambda omega: load2 if omega == omega_w else 0
w_x3 = lambda omega: load3 if omega == omega_w else 0
w_x4 = lambda omega: load4 if omega == omega_w else 0
w_x5 = lambda omega: load5 if omega == omega_w else 0

node18.add_load(x=load1)
node17.add_load(x=load2)
node15.add_load(x=load3)
node13.add_load(x=load4)
node5.add_load(x=load4)

# for beam in beams:
#     beam.AddDistributedLoad(z=1e2) # test for FRF


# define the traffic loading
q_v1 = lambda omega: -l2*(2+0.54)*1e3/2 if omega == omega_c else 0
q_v2 = lambda omega: -l3*(2+0.54)*1e3/2 if omega == omega_c else 0
s1.run_connectivity()


omega_range = np.linspace(1,200*2*np.pi,5000)
u_18 = np.zeros(len(omega_range),complex)
u_18_test = np.zeros(len(omega_range),complex)
for ii, omega in enumerate(omega_range):
# %% Get the global stiffness and force matrices
    Kc_global = s1.GlobalConstrainedStiffness(omega)
    Fc_global = s1.GlobalConstrainedForce(omega)
    
    # %% solve for free node displacements
    u_free = s1.SolveUfree(Kc_global, Fc_global)
    
    # %% post-processing
    u_18[ii] = u_free[19]
    

omega_test = np.array([4.58937,5.28712,9.8458,11.454,12.7329,18.4059,24.8094])
plt.figure()
plt.plot(omega_range,abs(u_18))
plt.plot(omega_range,u_18.real)
plt.plot(omega_range,u_18.imag)
plt.scatter(omega_test,np.zeros_like(omega_test))
plt.xlim([0,50]);


