# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:05:12 2025

@author: GREY
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:23:06 2024

@author: rensv
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA
import scipy.optimize as opt

# %% add package to sys and import 

import pydynsm as PDM


Assembler = PDM.Assembler

# %% Example

# %%% Initialise an assembler with your project name

s1 = Assembler('SSbeam',analysis_type='new')

# functions for root-finder

def BVP(ww):
    return s1.GlobalConstrainedStiffness(ww)

def det_func(ww):
    return np.linalg.det(BVP(ww))

def find_eigen_frequencies(f):
    omega = 2 * np.pi * f
    Det_M = np.array([np.linalg.det(BVP(ww)) for ww in omega])
    omega_initial = omega[np.where(np.isclose(abs(np.diff(np.angle(Det_M)))/np.pi,1, atol=.1))[0]]
    omega_m = []
    for ww in omega_initial:
        omega_m.append(opt.newton(det_func, ww).real)
    return omega_m



# %%% Parameters
E = 210e5
EA = 7e6
A = EA/E
EI = 1.5e7
I = EI/E
W = I
rhoA = 1e03 
rho = rhoA/A
L  = 1
omega = 1
ksi = 0.01
ksi = 0
omega_f = omega


node1 = s1.CreateNode(0,0)
# node2 will have no 'z' displacement - will be handled without setting stiffness to infinity (results are still to be verified)
node2 = s1.CreateNode(L/4,0)
node3 = s1.CreateNode(L,0) 

node1.fix_node('x', 'z')
node2.fix_node('x')
node3.fix_node('x','z')


elem = s1.CreateElement([node1, node2])
elem1 = s1.CreateElement([node2, node3])

elem.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': W, 'ksi': ksi})
elem1.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'rho':rho, 'Ib':I, 'Wb': W, 'ksi': ksi})


s1.run_connectivity()
f = np.linspace(1,2000,10001)

omega = 2 * np.pi * f
omega_m = np.array(find_eigen_frequencies(f))
N_modes = len(omega_m)
print(N_modes)

Det_M = np.array([np.linalg.det(BVP(ww)) for ww in omega])
arg_det = np.angle(Det_M)

fig, axs = plt.subplots(3, sharex=True,figsize=(10,6))
axs[0].plot(omega / 2 / np.pi,abs(Det_M), label='abs')
axs[0].set_yscale('log')
axs[1].plot(omega / 2 / np.pi,Det_M.real, label='Real')
axs[1].plot(omega / 2 / np.pi,Det_M.imag, label='Imag')
axs[2].plot(omega / 2/ np.pi, arg_det, label='argument')
axs[2].plot(omega[:-1] / 2/ np.pi, abs(np.diff(arg_det))/np.pi, label='abs phase jump / $\pi$')
axs[2].set_ylim([-np.pi, np.pi])

axs[0].scatter(omega_m / 2  / np.pi, np.ones_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[1].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
axs[2].scatter(omega_m / 2  / np.pi, np.zeros_like(omega_m),marker='x', c='r' , label='Eigenfrequencies')
for ax in axs.flat:
    ax.grid()
    ax.legend()
    ax.set_xlim([0,2000])

axs[0].set_ylim(bottom=1e-2)

axs[0].set_ylabel('Determinant')
axs[1].set_ylabel('Determinant / $\omega^6$')
axs[2].set_ylabel('phase ($-\pi$ to $\pi$)')
#plt.rcParams['figure.figsize'] = [5,10]
plt.tight_layout()

omega_1 = 4*np.pi**2/L**2*np.sqrt(E*I/rho/A)
print(f'1st natural frequencies analytical should be \n{omega_1}\n')

omega_test = omega_m[3]

F_1 = lambda omega: 0.001 if omega == omega_test else 0
node2.add_load(z=F_1)




# K_global = s1.GlobalStiffness(omega_test)
# det_test =  np.linalg.det(BVP(omega_test)/1e10)
Kc_global = s1.GlobalConstrainedStiffness(omega_test)
Fc_global = s1.GlobalConstrainedForce(omega_test)

u_free = s1.SolveUfree(Kc_global, Fc_global)
u_elem = s1.FullDisplacement(u_free)
disp = s1.ElementDisplacements(u_elem, omega_test)
s1.PlotElementDisplacements(disp,scale=1e-2)

