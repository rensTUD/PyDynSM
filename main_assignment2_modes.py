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

# %% add package to sys and import 

import pydynsm as PDM

Assembler = PDM.Assembler
s1 = Assembler('bridge',analysis_type='new')

# functions for root-finder

def BVP(ww):
    return s1.GlobalConstrainedStiffness(ww)

def det_func(ww):
  dets = np.linalg.det(BVP(ww)/1e5)
  # dets = np.heaviside(dets.real, 1) * np.log(abs(dets.real)) - np.heaviside(-dets.real, 0) * np.log(abs(dets.real))
  return dets

def find_eigen_frequencies(omega):
    """ Finds natural frequencies by detecting determinant phase changes and refining with root finding. """
    Det_M = np.array([det_func(ww) for ww in omega])
    # Identify initial guesses for natural frequencies
    index = np.isclose(abs(np.diff(np.angle(Det_M)))/np.pi,1, atol=.1)
    omega_initial = omega[np.where(np.isclose(abs(np.diff(np.angle(Det_M)))/np.pi,1, atol=.1))[0]]
    omega_m = []
    for ww in omega_initial:
        omega_m.append(opt.newton(det_func, ww))
    return omega_m
# %% inputs
l1 = 3
l2 = 11.75
l3 = 6

# material props
Es = 210e9 # Elasticity modulus steel
Ec = 3.1e10 # Elasticity modulus concrete
rhoc = 2500 # density concrete
rhos = 7800 # density steel
Ab = 4.8   # beam cross sectional area 
Ib = 0.0395 # beam 2nd mmt of intertial
rc = 9.5*0.001
Ac = np.pi*rc**2
Ic = 0.25*np.pi*rc**4
ksi = 0


# %% sturctural plot
# %%% nodes
node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(l2,0)
node3 = s1.CreateNode(2*l2,0)
node4 = s1.CreateNode(3*l2,0)
node5 = s1.CreateNode(4*l2,0)
node6 = s1.CreateNode(4*l2+l1,0)
node7 = s1.CreateNode(4*l2+l1,-l1)
node8 = s1.CreateNode(4*l2+l1+l3,0)
node9 = s1.CreateNode(4*l2+l1+2*l3,0)
node10 = s1.CreateNode(4*l2+l1+3*l3,0)
node11 = s1.CreateNode(4*l2+l1+4*l3,0)
node12 = s1.CreateNode(4*l2+l1,l1)
node13 = s1.CreateNode(4*l2,l1)
node14 = s1.CreateNode(4*l2+l1,2*l1)
node15 = s1.CreateNode(4*l2,2*l1)
node16 = s1.CreateNode(4*l2+l1,3*l1)
node17 = s1.CreateNode(4*l2,3*l1)
node18 = s1.CreateNode(4*l2,4*l1)
# s1.PlotStructure()
# %%% elements
# beam
beams = []

# Define node pairs for beam elements
beam_node_pairs = [
    (node1, node2),
    (node2, node3),
    (node3, node4),
    (node4, node5),
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
    (node5, node6),
    (node5, node7),   (node6, node7),
    (node5, node12),  (node6, node12),
    (node12, node13), (node12, node14),
    (node12, node15), (node13, node15),
    (node13, node5),  (node14, node15),
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
node7.fix_node('x','z')

# %% assign sections
for beam in beams:
    beam.SetSection('EulerBernoulli Beam', {'E': Ec, 'A':Ab, 'rho':rhoc, 'Ib':Ib, 'Wb': Ib,'ksi': ksi})
    beam.SetSection('Rod', {'E': Ec, 'A':Ab, 'rho':rhoc})
for cable in front_cables:
    # cable.SetSection('Rod', {'E': Es, 'A':Ac, 'rho':rhos})
    # cable.SetSection('EulerBernoulli Beam', {'E': Ec, 'A':Ab, 'rho':rhoc, 'Ib':Ib, 'Wb': Ib})
    cable.SetSection('EulerBernoulli Beam', {'E': Es, 'A':Ac, 'rho':rhos, 'Ib':Ic, 'Wb': Ic,'ksi': ksi})
    cable.SetSection('Rod', {'E': Es, 'A':Ac, 'rho':rhos,'ksi': ksi})
for cable in back_cables:
    cable.SetSection('EulerBernoulli Beam', {'E': Es, 'A':Ac, 'rho':rhos, 'Ib':Ic, 'Wb': Ic,'ksi': ksi})
    cable.SetSection('Rod', {'E': Es, 'A':Ac, 'rho':rhos,'ksi': ksi})
for truss in trusses:
    truss.SetSection('EulerBernoulli Beam', {'E': Ec, 'A':Ab, 'rho':rhoc, 'Ib':Ib, 'Wb': Ib,'ksi': ksi})
    truss.SetSection('Rod', {'E': Ec, 'A':Ab, 'rho':rhoc,'ksi': ksi})

s1.PlotStructure(plot_elements=True)
s1.run_connectivity()

omega_range = np.linspace(1e-4,100,3000)  # Adjust range as needed
omega_m = np.array(find_eigen_frequencies(omega_range))

dets = np.array([np.linalg.det(BVP(ww)/1e5) for ww in omega_range])
arg_dets = abs(np.diff(np.angle(dets)))/np.pi
log_dets = np.heaviside(dets.real, 1) * np.log(abs(dets.real)) - np.heaviside(-dets.real, 0) * np.log(abs(dets.real))
# # log_dets_imag = np.heaviside(dets.imag, 1) * np.log(abs(dets.imag)) - np.heaviside(-dets.imag, 0) * np.log(abs(dets.imag))
plt.figure()
plt.plot(omega_range, dets.real)
plt.scatter(omega_m, np.zeros_like(omega_m),marker='x', c='r')

plt.figure()
plt.plot(omega_range[1:], arg_dets)

omega_initial = omega_range[np.where(np.isclose(dets.real,0, atol=.1))[0]]
omega_initial1 = omega_range[np.where(np.isclose(abs(np.diff(np.angle(log_dets)))/np.pi,1, atol=.1))[0]]
# # # natural_frequencies = find_eigen_frequencies(omega_range)

# fig, axs = plt.subplots(3, sharex=True,figsize=(10,6))
# axs[0].plot(omega_range / 2 / np.pi,np.abs(dets), label='abs')
# axs[0].set_yscale('log')
# # axs[1].plot(omega_range / 2 / np.pi,dets.imag/omega_range**6, label='Imag')
# axs[1].plot(omega_range / 2 / np.pi,log_dets, label='Real')
# axs[2].plot(omega_range / 2/ np.pi, arg_dets, label='argument')
# axs[2].plot(omega_range[:-1] / 2/ np.pi, abs(np.diff(arg_dets))/np.pi, label='abs phase jump / $\pi$')
# axs[2].set_ylim([-np.pi, np.pi])

# axs[0].scatter(omega_m / 2  / np.pi, np.ones_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
# axs[1].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
# axs[2].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
# for ax in axs.flat:
#     ax.grid()
#     ax.legend()

# axs[0].set_ylabel('Determinant')
# axs[1].set_ylabel('Determinant / $\omega^6$')
# axs[2].set_ylabel('phase ($-\pi$ to $\pi$)')
# #plt.rcParams['figure.figsize'] = [5,10]
# plt.tight_layout()

# # %% constraint nodes
omega_test = omega_m[0]
u_1 = lambda omega: 1 if omega == omega_test else 0
F_1 = lambda omega: 0.001 if omega == omega_test else 0
node2.prescribe_node(z=u_1)
node3.add_load(z=F_1)


# K_global = s1.GlobalStiffness(omega_test)
# det_test =  np.linalg.det(BVP(omega_test)/1e10)
Kc_global = s1.GlobalConstrainedStiffness(omega_test)
Fc_global = s1.GlobalConstrainedForce(omega_test)

u_free = s1.SolveUfree(Kc_global, Fc_global)
u_elem = s1.FullDisplacement(u_free)
disp = s1.ElementDisplacements(u_elem, omega_test)
s1.PlotElementDisplacements(disp,scale=1e7)




#%% normalize modes

# export coordinates of structrues
#%% orthogonality check


