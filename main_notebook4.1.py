# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:33:47 2025

@author: GREY
"""

# import python standard dependencies

import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import numpy.fft as fft
import scipy.optimize as opt

# %% add package to sys and import 

import pydynsm as PDM

# all parameters
h = 6
b1 = 6
b2 = 8
E = 210e9
rho = 7850

EI_G = 2e5
I_G = EI_G/E

rhoA_G = 60
A_G = rhoA_G/rho


EI_C = 8e4
I_C = EI_C/E
rhoA_C = 40
A_C = rhoA_C/rho

f = 0.001
Omega = 2 * np.pi * f

P_0 = 70
Q_1 = 5
Assembler = PDM.Assembler
s1 = Assembler('Frame',analysis_type='new')

def BVP(ww):
    return s1.GlobalConstrainedStiffness(ww)

def det_func(ww):
  dets = np.linalg.det(BVP(ww)/1e5)
  # dets = np.heaviside(dets.real, 1) * np.log(abs(dets.real)) - np.heaviside(-dets.real, 0) * np.log(abs(dets.real))
  return dets

def find_eigen_frequencies(omega):
    """ Finds natural frequencies by detecting determinant phase changes and refining with root finding. """
    Det_M = np.array([det_func(ww) for ww in omega])
    omega_initial = omega[np.where(np.isclose(abs(np.diff(np.angle(Det_M)))/np.pi,1, atol=.1))[0]]
    omega_m = []
    for ww in omega_initial:
        omega_m.append(opt.newton(det_func,ww).real)
    return np.unique(omega_m)

def remove_close_roots(roots, tol=1e-6):
    """ Remove repetitive roots found by the eigen-solver """
    roots = np.sort(roots)  # 先排序，方便處理
    filtered_roots = [roots[0]]  # 保留第一個 root

    for i in range(1, len(roots)):
        if not np.isclose(roots[i], filtered_roots[-1], atol=tol):
            filtered_roots.append(roots[i])

    return np.array(filtered_roots)
# %%% Initialise an assembler with your project name

nodes = []

nodes.append(s1.CreateNode(0,0))
nodes.append(s1.CreateNode(0,h))
nodes.append(s1.CreateNode(b2,h))
nodes.append(s1.CreateNode(b2,0))


s1.PlotStructure()

elements = []

elements.append(s1.CreateElement([nodes[0],nodes[1]]))
elements.append(s1.CreateElement([nodes[1],nodes[2]]))
elements.append(s1.CreateElement([nodes[2],nodes[3]]))
s1.PlotStructure(plot_elements=True)


nodes[0].fix_node('x','z','phi_y')
nodes[3].fix_node('x','z','phi_y')


# assign parameters
column = {}
column['E'] = E
column['A'] = A_C
column['Ib'] = I_C
column['Wb'] = I_C
column['ksi'] = 0
column['rho'] = rho

beam = {}
beam['E'] = E
beam['A'] = A_G
beam['Ib'] = I_G
beam['Wb'] = I_G
beam['ksi'] = 0
beam['rho'] = rho

rodb = {}
rodb['E'] = 1*E
rodb['A'] = A_G
rodb['ksi'] = 0
rodb['rho'] = rho

rodc = {}
rodc['E'] = 1*E
rodc['A'] = A_C
rodc['ksi'] = 0
rodc['rho'] = rho


elements[0].SetSection('EulerBernoulli Beam', column)
elements[2].SetSection('EulerBernoulli Beam', column)
elements[1].SetSection('EulerBernoulli Beam', beam)



elements[0].SetSection('Rod', rodc)
elements[2].SetSection('Rod', rodc)
elements[1].SetSection('Rod', rodb)

s1.run_connectivity()

f = np.linspace(1e-4, 200, 10001)
omega_range = 2* np.pi * f
omega_m = np.array(find_eigen_frequencies(omega_range))
Det_M = np.array([det_func(omega) for omega in omega_range])
arg_dets = np.angle(Det_M)

fig, axs = plt.subplots(3, sharex=True,figsize=(10,6))
axs[0].plot(omega_range / 2 / np.pi,np.abs(Det_M), label='abs')
axs[0].set_yscale('log')
# axs[1].plot(omega_range / 2 / np.pi,dets.imag/omega_range**6, label='Imag')
axs[1].plot(omega_range / 2 / np.pi,Det_M.real/omega_range**6/1e43, label='real')
axs[1].plot(omega_range / 2 / np.pi,Det_M.imag/omega_range**6/1e43, label='imag')
axs[1].set_yscale('log')
axs[2].plot(omega_range / 2/ np.pi, arg_dets, label='argument')
axs[2].plot(omega_range[:-1] / 2/ np.pi, abs(np.diff(arg_dets))/np.pi, label='abs phase jump / $\pi$')
axs[2].set_ylim([-np.pi, np.pi])

axs[0].scatter(omega_m / 2  / np.pi, np.ones_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[1].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[2].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
for ax in axs.flat:
    ax.grid()
    ax.legend()
    ax.set_xlim([0,50])

axs[0].set_ylim(bottom=1e-1)
axs[0].set_ylabel('Determinant')
axs[1].set_ylabel('Determinant / $\omega^6$')
axs[2].set_ylabel('phase ($-\pi$ to $\pi$)')
#plt.rcParams['figure.figsize'] = [5,10]
plt.tight_layout()
def plot_modes(omega_m,n=5):
    """
    n is the maximum number of modes 
    """
    omega_m = remove_close_roots(omega_m, tol=1e-6)[:n]
    N_modes = len(omega_m)
    for omega in omega_m:
        F_1 = lambda omega: 0.001 if omega == omega else 0
        nodes[1].add_load(z=F_1)
        Kc_global = s1.GlobalConstrainedStiffness(omega)
        Fc_global = s1.GlobalConstrainedForce(omega)
        u_free = s1.SolveUfree(Kc_global, Fc_global)
        u_elem = s1.FullDisplacement(u_free)
        disp = s1.ElementDisplacements(u_elem, omega,500)
        s1.PlotElementDisplacements(disp,scale=1/np.max(np.abs(u_elem)))
    return 

plot_modes(omega_m,6)
n = 15 #number of modes you are interested in
omega_used = remove_close_roots(omega_m, tol=1e-6)[:n]

disp_dict = {}
num = 8000
for omega in omega_used:
    F_1 = lambda omega: 0.001 if omega == omega else 0
    nodes[1].add_load(z=F_1)
    Kc_global = s1.GlobalConstrainedStiffness(omega)
    Fc_global = s1.GlobalConstrainedForce(omega)
    u_free = s1.SolveUfree(Kc_global, Fc_global)
    u_elem = s1.FullDisplacement(u_free)
    disp = s1.ElementDisplacements(u_elem, omega,num)
    disp_dict[omega] = disp


def get_local_coordinates(elements):
    coor_local = {}
    for element in elements:
        coor_local[element.id] = element.L*np.linspace(0,1,num)
    return coor_local

coor_local = get_local_coordinates(elements)
    
def Normalize_modes(omega_used,disp):
    Gamma = {}
    disp_bar = {}
    for omega in omega_used:
        Gamma[omega] =  rho*A_C * np.trapz(disp_dict[omega][elements[0].id][0]**2,coor_local[elements[0].id]) +\
                        rho*A_G * np.trapz(disp_dict[omega][elements[1].id][1]**2,coor_local[elements[1].id]) +\
                        rho*A_C * np.trapz(disp_dict[omega][elements[2].id][0]**2,coor_local[elements[2].id]) +\
                        rho*A_G * (b2)*disp_dict[omega][elements[0].id][0][-1]**2
        
        disp_bar[omega] = {elem: disp / np.sqrt(Gamma[omega]) for elem, disp in disp_dict[omega].items()}
    return disp_bar

disp_bar = Normalize_modes(omega_used,disp_dict)

s1.PlotElementDisplacements(disp_bar[omega_used[0]],scale=10)
s1.PlotElementDisplacements(disp_bar[omega_used[1]],scale=10)
s1.PlotElementDisplacements(disp_bar[omega_used[2]],scale=10)
s1.PlotElementDisplacements(disp_bar[omega_used[3]],scale=10)
s1.PlotElementDisplacements(disp_bar[omega_used[4]],scale=10)
s1.PlotElementDisplacements(disp_bar[omega_used[5]],scale=10)


def Check_Orthogonality(omega_used,disp_bar,coor_local):
    # Get number of omega values
    num_omega = len(omega_used)
    
    # Initialize the Ortho matrix
    Ortho_matrix = np.zeros((num_omega, num_omega))
    
    # Convert omega_used to a list (if needed for indexing)
    omega_list = list(omega_used)
    
    # Iterate over indices instead of dictionary keys
    for i, omega_1 in enumerate(omega_used):
        for j, omega_2 in enumerate(omega_used):
            
            # Compute each term separately for clarity
            term1 = rho * A_C * np.trapz(
                disp_bar[omega_1][elements[0].id][0] * disp_bar[omega_2][elements[0].id][0],
                coor_local[elements[0].id]
            )
    
            term2 = rho * A_G * np.trapz(
                disp_bar[omega_1][elements[1].id][1] * disp_bar[omega_2][elements[1].id][1],
                coor_local[elements[1].id]
            )
    
            term3 = rho * A_C * np.trapz(
                disp_bar[omega_1][elements[2].id][0] * disp_bar[omega_2][elements[2].id][0],
                coor_local[elements[2].id]
            )
    
            term4 = rho * A_G * b2 * disp_bar[omega_1][elements[0].id][0][-1] * disp_bar[omega_2][elements[0].id][0][-1]
    
            # Assign to matrix
            Ortho_matrix[i, j] = term1 + term2 + term3 + term4
    return Ortho_matrix

Ortho_coeff = Check_Orthogonality(omega_used,disp_bar,coor_local)
plt.imshow(abs(Ortho_coeff))
plt.colorbar()

    







                   


    

    
        



    
    



