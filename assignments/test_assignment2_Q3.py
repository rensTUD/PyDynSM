# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:59:03 2025

@author: GREY
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA
import scipy.optimize as opt
import pydynsm as PDM
import scipy.integrate as spi

Assembler = PDM.Assembler
s1 = Assembler('bridge')
omega_m = np.loadtxt("eigenvalues.csv", delimiter=",")
#%% function for calculating wind load
def trap(a,b):

    '''
    a,b are start and end locations of the wind load, defined in the local direction.
    '''
    
    alpha = 0.4
    qw = 20  
    result, error = spi.quad(lambda x: np.exp(alpha * x), a, b)
    return result*qw

# %% inputs
l1 = 4
l2 = 11.75
l3 = 5.33

# material props
E = 210e9  # Elasticity modulus steel
rho =  7800
# for beam
I_b = 0.040 + 0.003*9
A_b = 0.5


# for cable
r_c = 0.0475+0.001*9
A_c = np.pi*r_c**2

# for truss 
A_t = 0.05 # truss erea

k_d = 1e6
# c_d = 1e6
c_d = 0
ksi = 0.01

EA_b = E*A_b
EI_b = E*I_b
rhoA_b = rho*A_b
EA_c = E*A_c
rhoA_c = rho*A_c
EA_t = E*A_t
rhoA_t = rho*A_t
# %% sturctural plot
# %%% nodes
node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(l2,0)
node3 = s1.CreateNode(2*l2,0)
node4 = s1.CreateNode(3*l2,0)
node5 = s1.CreateNode(4*l2,0)
node6 = s1.CreateNode(4*l2+l1,0)
node7 = s1.CreateNode(4*l2+l1,-l1,dof_config = ['x','z'])
node8 = s1.CreateNode(4*l2+l1+l3,0)
node9 = s1.CreateNode(4*l2+l1+2*l3,0,)
node10 = s1.CreateNode(4*l2+l1+3*l3,0)
node11 = s1.CreateNode(4*l2+l1+4*l3,0)
node12 = s1.CreateNode(4*l2+l1,l1,dof_config = ['x','z'])
node13 = s1.CreateNode(4*l2,l1,dof_config = ['x','z'])
node14 = s1.CreateNode(4*l2+l1,2*l1,dof_config = ['x','z'])
node15 = s1.CreateNode(4*l2,2*l1,dof_config = ['x','z'])
node16 = s1.CreateNode(4*l2+l1,3*l1,dof_config = ['x','z'])
node17 = s1.CreateNode(4*l2,3*l1,dof_config = ['x','z'])
node18 = s1.CreateNode(4*l2,4*l1,dof_config = ['x','z'])
s1.PlotStructure()
# %%% elements
# beam
beams = []

# Define node pairs for beam elements
beam_node_pairs = [
    (node1, node2),
    (node2, node3),
    (node3, node4),
    (node4, node5),
    (node5, node6),
    (node6, node8),
    (node8, node9),
    (node9, node10),
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

cables = front_cables + back_cables

# %% constraint nodes
node1.fix_node('z')
node7.fix_node('x','z')

# %% assign sections
for beam in beams[:-1]:
    beam.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':ksi, 'Wb':I_b})
    beam.SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':ksi})

beams[-1].SetSection('EulerBernoulli Beam with foundation', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':ksi, 'kd':k_d,'cd':c_d})
beams[-1].SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':ksi})

for cable in cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':ksi})
for truss in trusses:
    truss.SetSection('Rod', {'E': E, 'A':A_t, 'rho':rho,'ksi':ksi})
s1.run_connectivity()
#%% define the real wind load here
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

# %% plot wind load along height
plt.figure()
plt.plot(np.array([load1, load2, load3, load4, load5])/load1,np.array([0,l1,2*l1,3*l1,4*l1]))
plt.ylabel('Height [m]')
plt.xlabel('Equivalent Wind load [kN]')
plt.grid()
plt.show()

node18.add_load(x=load5/load1)
node17.add_load(x=load4/load1)
node15.add_load(x=load3/load1)
node13.add_load(x=load2/load1)
node5.add_load(x=load1/load1)

omega = np.linspace(1, 200*2*np.pi,500)
u_18 = np.zeros(len(omega),complex)
u_5 = np.zeros(len(omega),complex)
for ii, ww in enumerate(omega):
    Kc_global = s1.GlobalConstrainedStiffness(ww)
    Fc_global = s1.GlobalConstrainedForce(ww)
    u_free = s1.SolveUfree(Kc_global, Fc_global)
    u_elem = s1.FullDisplacement(u_free)
    u_18[ii] = s1.ElementDisplacements(u_elem,ww)[trusses[0].id][0][0]
    u_5[ii] = s1.ElementDisplacements(u_elem,ww)[beams[3].id][0][-1]

# %% plot FRF
fig, axs = plt.subplots(2, sharex=True,figsize=(10,6))
axs[0].plot(omega/2/np.pi,u_18.real,label='Real')
axs[0].plot(omega/2/np.pi,u_18.imag,label='Imag')
axs[0].plot(omega/2/np.pi,abs(u_18),label='Abs')
axs[1].plot(omega/2/np.pi,u_5.real,label='Real')
axs[1].plot(omega/2/np.pi,u_5.imag,label='Imag')
axs[1].plot(omega/2/np.pi,abs(u_5),label='Abs')
axs[0].scatter(omega_m/2/np.pi,np.ones_like(omega_m)*np.min(abs(u_18)),marker='x', c='r',label='Eigenvalues')
axs[1].scatter(omega_m/2/np.pi,np.ones_like(omega_m)*np.min(abs(u_5)),marker='x', c='r',label='Eigenvalues')
axs[1].set_xlabel('$f$[Hz]')
axs[0].set_ylabel(r'FRF $\tilde{U}_{18}$ [m/N]')
axs[1].set_ylabel(r'FRF $\tilde{U}_5$ [m/N]')
for ax in axs.flat:
    ax.grid()
    ax.legend()

