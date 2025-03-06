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

Assembler = PDM.Assembler
s1 = Assembler('bridge')

# functions for root-finder
def BVP(ww):
    return s1.GlobalConstrainedStiffness(ww)

def det_func(ww):
  dets = np.linalg.det(BVP(ww/1e10))
  # dets = np.heaviside(dets.real, 1) * np.log(abs(dets)) - np.heaviside(-dets.real, 0) * np.log(abs(dets))
  return dets

def find_eigen_frequencies(omega):
    Det_M = np.array([det_func(ww) for ww in omega])
    omega_initial = omega[np.where(np.isclose(abs(np.diff(np.angle(Det_M)))/np.pi,1, atol=.1))[0]]
    # omega_initial = omega[np.where(np.isclose(Det_M,1,atol=.1))[0]]
    omega_m = []
    for ww in omega_initial:
        omega_m.append(opt.newton(det_func,ww).real)
    return np.unique(omega_m)

def remove_close_roots(roots, tol=1e-6):
    '''
    The function removes close roots from the array of root-finders
    '''
    roots = np.sort(roots)
    filtered_roots = [roots[0]]

    for i in range(1, len(roots)):
        if not np.isclose(roots[i], filtered_roots[-1], atol=tol):
            filtered_roots.append(roots[i])

    return np.array(filtered_roots)
# %% inputs
l1 = 4
l2 = 11.75
l3 = 5.33

# material props
E = 210e9  # Elasticity modulus steel
rho =  7800
# for beam
I_b = 0.0395
A_b = 0.5

# for cable
r_c = 9.5*0.005
A_c = np.pi*r_c**2
I_c = 0.25*np.pi*r_c**4

# for truss 
A_t = 0.1*A_b # truss erea
I_t = I_c


k_d = 1e6
# c_d = 1e6
c_d = 0
# %% sturctural plot
# %%% nodes
node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(l2,0)
node3 = s1.CreateNode(2*l2,0)
node4 = s1.CreateNode(3*l2,0)
node5 = s1.CreateNode(4*l2,0)
node6 = s1.CreateNode(4*l2+l1,0)
node7 = s1.CreateNode(4*l2+l1,-l1,dof_config = [['x'],['z']])
node8 = s1.CreateNode(4*l2+l1+l3,0)
node9 = s1.CreateNode(4*l2+l1+2*l3,0,)
node10 = s1.CreateNode(4*l2+l1+3*l3,0)
node11 = s1.CreateNode(4*l2+l1+4*l3,0)
node12 = s1.CreateNode(4*l2+l1,l1,dof_config = [['x'],['z']])
node13 = s1.CreateNode(4*l2,l1,dof_config = [['x'],['z']])
node14 = s1.CreateNode(4*l2+l1,2*l1,dof_config = [['x'],['z']])
node15 = s1.CreateNode(4*l2,2*l1,dof_config = [['x'],['z']])
node16 = s1.CreateNode(4*l2+l1,3*l1,dof_config = [['x'],['z']])
node17 = s1.CreateNode(4*l2,3*l1,dof_config = [['x'],['z']])
node18 = s1.CreateNode(4*l2,4*l1,dof_config = [['x'],['z']])
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
node7.fix_node('x','z')

# %% assign sections
for beam in beams:
    beam.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0, 'Wb':I_b})
    beam.SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0})

# beams[-1].SetSection('EulerBernoulli Beam with foundation', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0, 'kd':k_d,'cd':c_d})
# beams[-1].SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0})
for cable in front_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':0})
    # cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':0})
for cable in back_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':0})
    # cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':0})
for truss in trusses:
    truss.SetSection('Rod', {'E': E, 'A':A_t, 'rho':rho,'ksi':0})
    # truss.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_t, 'rho':rho, 'Ib':I_t,'ksi':0})
s1.run_connectivity()

omega_range = np.linspace(1,100,5000)
# plot the det
Det_M = np.array([det_func(omega) for omega in omega_range])
dets = np.heaviside(Det_M.real, 1) * np.log(abs(Det_M.real)) - np.heaviside(-Det_M.real, 0) * np.log(abs(Det_M.real))
arg_dets = np.angle(dets)
plt.figure()
plt.plot(omega_range,dets);
omega_m = find_eigen_frequencies(omega_range)
plt.scatter(omega_m,np.ones_like(omega_m),marker='x', c='r');

#%% visualize the eigenvalues
omega_m = omega_m[:30]
omega = omega_range

fig, axs = plt.subplots(3, sharex=True,figsize=(10,6))
axs[0].plot(omega / 2 / np.pi,abs(Det_M), label='abs')
axs[0].set_yscale('log')
axs[1].plot(omega / 2 / np.pi,Det_M.real/omega**6, label='Real')
axs[1].plot(omega / 2 / np.pi,Det_M.imag/omega**6, label='Imag')
axs[2].plot(omega / 2/ np.pi, arg_dets, label='argument')
axs[2].plot(omega[:-1] / 2/ np.pi, abs(np.diff(arg_dets))/np.pi, label='abs phase jump / $\pi$')
axs[2].set_ylim([-np.pi, np.pi])

axs[0].scatter(omega_m / 2  / np.pi, np.ones_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[1].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[2].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
for ax in axs.flat:
    ax.grid()
    ax.legend()
    ax.set_xlim([0,8])

axs[0].set_ylim(bottom=1e-1)
plt.tight_layout()

# #%% compute eigenmodes

# for omega in omega_m[:6]:
#     Kc_global = s1.GlobalConstrainedStiffness(omega)
#     # solve for the eigenvector by setting one component to amplitude 1
#     u_free = s1.SolveUfree(Kc_global[:-1,:-1], -Kc_global[:-1,-1])
#     u_free = np.concatenate((u_free, np.array([1])))
#     u_elem = s1.FullDisplacement(u_free)
#     disp = s1.ElementDisplacements(u_elem,omega)
#     s1.PlotElementDisplacements(disp)

# def get_local_coordinates(elements,num):
#     coor_local = {}
#     for element in elements:
#         coor_local[element.id] = element.L*np.linspace(0,1,num)
#     return coor_local
# # %% disp dict
# disp_dict = {}
# for omega in omega_m[:6]:
#     Kc_global = s1.GlobalConstrainedStiffness(omega)
#     u_free = s1.SolveUfree(Kc_global[:-1,:-1], -Kc_global[:-1,-1])
#     u_free = np.concatenate((u_free, np.array([1])))
#     u_elem = s1.FullDisplacement(u_free)
#     disp = s1.ElementDisplacements(u_elem, omega)
#     disp_dict[omega] = disp
# local_coordinates = get_local_coordinates(beams, 20)
# #%% normalize modes
# Gamma_beam = np.zeros(len(beams),complex)
# Gamma_fc = np.zeros(len(front_cables),complex)
# Gamma_bc = np.zeros(len(back_cables),complex)
# Gamma_truss = np.zeros(len(trusses),complex)
# def Normalize_modes(Omega,displacement,local_coor_b,local_coor_fc,local_coor_bc,local_coor_t):
#     disp_norm = {}
#     for omega in Omega:
#         for ii, beam in enumerate(beams):
#             Gamma_beam[ii] = rho*A_b * (np.trapz(displacement[omega][beam.id][0]**2,local_coor_b[beam.id] +
#                                         np.trapz(displacement[omega][beam.id][1]**2,local_coor_b[beam.id]))          
#     )
#         for jj, fcable in enumerate(front_cables):
#                 Gamma_fc[ii] = rho*A_c * (np.trapz(displacement[omega][fcable.id][0]**2,local_coor_fc[fcable.id] +
#                                           np.trapz(displacement[omega][fcable.id][1]**2,local_coor_fc[fcable.id]))          
#         )
#         for kk, bcable in enumerate(back_cables):
#                 Gamma_bc[ii] = rho*A_c * (np.trapz(displacement[omega][bcable.id][0]**2,local_coor_bc[bcable.id] +
#                                           np.trapz(displacement[omega][bcable.id][1]**2,local_coor_bc[bcable.id]))          
#         )
#         for ll,truss in enumerate(trusses):
#                 Gamma_truss[ii] = rho*A_b * (np.trapz(displacement[omega][truss.id][0]**2,local_coor_t[truss.id] +
#                                             np.trapz(displacement[omega][truss.id][1]**2,local_coor_[truss.id]))          
#         )
        
#         # disp_norm[omega] = {elem: disp / np.sqrt(np.sum(Gamma_beam)) for elem, disp in displacement[omega].items()}
#     return np.sum()

# np.sum(Normalize_modes(omega_m[:6],disp_dict,local_coordinates))





