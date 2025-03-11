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
  dets = np.linalg.det(BVP(ww)/1e10)
  # dets = np.heaviside(dets.real, 1) * np.log(abs(dets.real)) - np.heaviside(-dets.real, 0) * np.log(abs(dets.real))
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
    beam.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0, 'Wb':I_b})
    beam.SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0})

beams[-1].SetSection('EulerBernoulli Beam with foundation', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0, 'kd':k_d,'cd':c_d})
beams[-1].SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0})

for cable in cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':0})
for truss in trusses:
    truss.SetSection('Rod', {'E': E, 'A':A_t, 'rho':rho,'ksi':0})
s1.run_connectivity()

omega = np.linspace(1,200*2*np.pi,10000)
omega_m = remove_close_roots(find_eigen_frequencies(omega))

#%% visualize the plot of roots
Det_M = np.array([det_func(ww) for ww in omega])
arg_dets = np.angle(Det_M)
fig, axs = plt.subplots(3, sharex=True,figsize=(10,6))
axs[0].plot(omega / 2 / np.pi,abs(Det_M), label='abs')
axs[0].set_yscale('log')
axs[1].plot(omega / 2 / np.pi,Det_M.real/omega**10, label='Real')
axs[1].plot(omega / 2 / np.pi,Det_M.imag/omega**10, label='Imag')
axs[2].plot(omega / 2/ np.pi, arg_dets, label='argument')
axs[2].plot(omega[:-1] / 2/ np.pi, abs(np.diff(arg_dets))/np.pi, label='abs phase jump / $\pi$')
axs[2].set_ylim([-np.pi, np.pi])

axs[0].scatter(omega_m / 2  / np.pi, np.ones_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[1].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[2].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
for ax in axs.flat:
    ax.grid()
    ax.legend()

#%% normalize modes
num_points = 2000
local_coor_beam = np.array([np.linspace(0,1,num_points)*beam.L for beam in beams])
local_coor_cable = np.array([np.linspace(0,1,num_points)*cable.L for cable in cables])
local_coor_truss = np.array([np.linspace(0,1,num_points)*truss.L for truss in trusses])

#%%% firstly create a dictionary for storing the displacements of each mode
disp_dict = {}
for omega in omega_m[:30]:
    Kc_global = s1.GlobalConstrainedStiffness(omega)
    u_free = s1.SolveUfree(Kc_global[:-1,:-1], -Kc_global[:-1,-1])
    u_free = np.concatenate((u_free, np.array([1])))
    u_elem = s1.FullDisplacement(u_free)
    disp = s1.ElementDisplacements(u_elem, omega,num_points)
    disp_dict[omega] = disp
#%% work out the gamma factors
Gamma_beam = np.zeros((len(omega_m[:30]),len(beams)),complex)
Gamma_truss = np.zeros((len(omega_m[:30]),len(trusses)),complex)
Gamma_cable = np.zeros((len(omega_m[:30]),len(cables)),complex)
Gamma_beam_winkler = np.zeros((len(omega_m[:30]),1),complex)
for ii, ww in enumerate(omega_m[:30]):
    for jj, bb in enumerate(beams):
        Gamma_beam[ii,jj] = (rho*A_b*np.trapz(disp_dict[ww][bb.id][0]**2,local_coor_beam[jj])+
                                rho*A_b*np.trapz(disp_dict[ww][bb.id][1]**2,local_coor_beam[jj]))
    for kk, tt in enumerate(trusses):
        Gamma_truss[ii,kk] = rho*A_t*np.trapz(disp_dict[ww][tt.id][0]**2,local_coor_truss[kk])
    for ll, cc in enumerate(cables):
        Gamma_cable[ii,ll] = rho*A_c*np.trapz(disp_dict[ww][cc.id][0]**2,local_coor_cable[ll])
Gamma = (np.sum(Gamma_beam,axis=1,keepdims=True)+
         np.sum(Gamma_truss,axis=1,keepdims=True)+
         np.sum(Gamma_cable,axis=1,keepdims=True))
#%% normalize displacements
disp_norm = {}
for ii, omega in enumerate(omega_m[:30]):
    disp_norm[omega] = {elem: disp / np.sqrt(Gamma[ii]) for elem, disp in disp_dict[omega].items()}
#%% plot normalized displacements
for omega in omega_m[:5]:
    s1.PlotElementDisplacements(disp_norm[omega],scale=1e3)
#%% check orthogonality
N = 30
Ortho_beam = np.zeros((N,N,len(beams)), dtype=complex)
Ortho_truss = np.zeros((N,N,len(trusses)), dtype=complex)
Ortho_cable = np.zeros((N,N,len(cables)), dtype=complex)
for ii, ww1 in enumerate(omega_m[:N]):
        for jj, ww2 in enumerate(omega_m[:N]):
            for kk, bb in enumerate(beams):
             Ortho_beam[ii,jj,kk] = (rho*A_b*np.trapz(disp_norm[ww1][bb.id][0]*disp_norm[ww2][bb.id][0],local_coor_beam[kk])+
                                     rho*A_b*np.trapz(disp_norm[ww1][bb.id][1]*disp_norm[ww2][bb.id][1],local_coor_beam[kk]))
            for ll, cc in enumerate(cables):
             Ortho_cable[ii,jj,ll] = (rho*A_c*np.trapz(disp_norm[ww1][cc.id][0]*disp_norm[ww2][cc.id][0],local_coor_cable[ll]))
            for mm, tt in enumerate(trusses):
              Ortho_truss[ii,jj,mm] = (rho*A_t*np.trapz(disp_norm[ww1][tt.id][0]*disp_norm[ww2][tt.id][0],local_coor_truss[mm]))
    
Ortho = np.sum(Ortho_beam,axis=2)+ np.sum(Ortho_cable,axis=2)+ np.sum(Ortho_truss,axis=2)
plt.figure()
plt.imshow(abs(Ortho))
plt.colorbar();


#%% write first 30 omega_m as .csv.file
np.savetxt("eigenvalues.csv", omega_m, delimiter=",")
