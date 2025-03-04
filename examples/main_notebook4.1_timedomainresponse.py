# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 23:56:31 2025

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

f = 20
Omega = 2 * np.pi * f

P_0 = 70
Q_1 = 5
Assembler = PDM.Assembler
s1 = Assembler('Frame',analysis_type='new')

# %%% Initialise an assembler with your project name

nodes = []

nodes.append(s1.CreateNode(0,0))
nodes.append(s1.CreateNode(0,h))
nodes.append(s1.CreateNode(b1,h))
nodes.append(s1.CreateNode(b2,h))
nodes.append(s1.CreateNode(b2,0))

s1.PlotStructure()

elements = []

elements.append(s1.CreateElement([nodes[0],nodes[1]]))
elements.append(s1.CreateElement([nodes[1],nodes[2]]))
elements.append(s1.CreateElement([nodes[2],nodes[3]]))
elements.append(s1.CreateElement([nodes[3],nodes[4]]))
s1.PlotStructure(plot_elements=True)

nodes[0].fix_node('x','z','phi_y')
nodes[4].fix_node('x','z','phi_y')


# assign parameters
column = {}
column['E'] = E
column['A'] = A_C
column['Ib'] = I_C
column['Wb'] = I_C
column['ksi'] = 0.05
column['rho'] = rho

beam = {}
beam['E'] = E
beam['A'] = A_G
beam['Ib'] = I_G
beam['Wb'] = I_G
beam['ksi'] = 0.05
beam['rho'] = rho

rodb = {}
rodb['E'] = 1*E
rodb['A'] = A_G
rodb['ksi'] = 0.05
rodb['rho'] = rho

rodc = {}
rodc['E'] = 1*E
rodc['A'] = A_C
rodc['ksi'] = 0.05
rodc['rho'] = rho


elements[0].SetSection('EulerBernoulli Beam', column)
elements[2].SetSection('EulerBernoulli Beam', beam)
elements[1].SetSection('EulerBernoulli Beam', beam)
elements[3].SetSection('EulerBernoulli Beam', column)



elements[0].SetSection('Rod', rodc)
elements[2].SetSection('Rod', rodb)
elements[1].SetSection('Rod', rodb)
elements[3].SetSection('Rod', rodc)

# s1.run_connectivity()

# P = lambda omega: -P_0  if omega == Omega else 0
# q = lambda omega: Q_1  if omega == Omega else 0

# elements[0].AddDistributedLoad(x=q)
# nodes[2].add_load(z=P)

# s1.run_connectivity()

# # %% Get the global stiffness and force matrices

# K_global = s1.GlobalStiffness(Omega)
# F_global = s1.GlobalForce(Omega)
# Kc_global = s1.GlobalConstrainedStiffness(Omega)
# Fc_global = s1.GlobalConstrainedForce(Omega)
# # %% solve for free node displacements
# u_free = s1.SolveUfree(Kc_global, Fc_global)
# # f_supp = s1.SupportReactions(s1.GlobalStiffness(omega_f), u_free, s1.GlobalForce(omega_f))

# # %% post-processing
# u_elem = s1.FullDisplacement(u_free)
# disp = s1.ElementDisplacements(u_elem, Omega,num_points=20)
# s1.PlotElementDisplacements(disp,scale=1e4)

def force_time(Amp,f0,t0,t): 
    return Amp*(1-np.cos(2*np.pi*f0*(t-t0)))*np.sin(2*np.pi*f0*(t-t0))*(t>t0)*((t-t0)<1/f0)

t=np.linspace(0,20,10000,endpoint=False)
N = len(t)
q_load_time = force_time(80,10,0.15,t) 
P_load_time = force_time(600,5,0.2,t) 
Q_omega = fft.rfft(q_load_time)   * 2  / N
P_omega = fft.rfft(P_load_time)   * 2  / N

freq = fft.rfftfreq(N,t[1])
fig,axs = plt.subplots(2,2,figsize=(15,10))
axs[0,0].plot(t,q_load_time)
axs[0,0].set_xlim([0,1])
axs[1,0].plot(freq,abs(Q_omega))
axs[1,0].plot(freq,np.real(Q_omega))
axs[1,0].plot(freq,np.imag(Q_omega))
axs[0,1].plot(t,P_load_time)
axs[0,1].set_xlim([0,1])
axs[1,1].plot(freq,abs(P_omega))
axs[1,1].plot(freq,np.real(P_omega))
axs[1,1].plot(freq,np.imag(P_omega))
for ax in axs.flat:
    ax.grid()
axs[0,0].set_xlabel('t')
axs[0,0].set_ylabel('q')
axs[1,0].set_xlabel('f')
axs[1,0].set_ylabel('$Q_0$')
axs[1,0].legend(['abs','Re', 'Im'])
omega = 2 * np.pi * freq
omega[0] = 1e-3 # Approach 0 frequency that gives a numerical error, because of defiding by 0.
freq = omega/2/np.pi


# Define a function that behaves like lambdify
def Q_omega_callable(omega):
    """
    A callable function that mimics sympy.lambdify for FFT-based loading functions.
    
    Parameters:
    - omega (float): The frequency at which to evaluate the function.

    Returns:
    - The FFT result at the given frequency if an exact match is found.
    
    Raises:
    - ValueError if the requested frequency is not available in the FFT bins.
    """
    frequency = omega/2/np.pi
    if frequency in freq:
        idx = np.where(freq == frequency)[0][0]  # Get index of the exact match
        return Q_omega[idx]
    else:
        raise ValueError(f"Requested frequency {omega} not found in FFT bins!")
def P_omega_callable(omega):
    """
    A callable function that mimics sympy.lambdify for FFT-based loading functions.
    
    Parameters:
    - omega (float): The frequency at which to evaluate the function.

    Returns:
    - The FFT result at the given frequency if an exact match is found.
    
    Raises:
    - ValueError if the requested frequency is not available in the FFT bins.
    """
    frequency = omega/2/np.pi
    if frequency in freq:
        idx = np.where(freq == frequency)[0][0]  # Get index of the exact match
        return P_omega[idx]
    else:
        raise ValueError(f"Requested frequency {omega} not found in FFT bins!")

s1.run_connectivity()
elements[0].AddDistributedLoad(x=Q_omega_callable)
nodes[2].add_load(z=P_omega_callable)
FRF_midpoint_column = []
FRF_midpoint_beam = []

# need to further check better ways of dealing with broadband excitations:
for omega_value, p, q in zip(omega, P_omega, Q_omega):
    P = lambda omega: p  if omega == omega_value else 0
    Q = lambda omega: q  if omega == omega_value else 0
    elements[0].AddDistributedLoad(x=Q)
    nodes[2].add_load(z=P)
    Kc_global = s1.GlobalConstrainedStiffness(omega_value)
    Fc_global = s1.GlobalConstrainedForce(omega_value)
    u_free = s1.SolveUfree(Kc_global, Fc_global)
    u_elem = s1.FullDisplacement(u_free)
    disp = s1.ElementDisplacements(u_elem, omega_value,num_points=4001)
    ucol = disp[elements[0].id][0]
    wbeam = disp[elements[1].id][1]
    FRF_midpoint_column.append(ucol[2000])
    FRF_midpoint_beam.append(wbeam[2000])
    
FRF_midpoint_column = np.array(FRF_midpoint_column)
FRF_midpoint_beam = np.array(FRF_midpoint_beam)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(14,10))
axs[0].plot(freq,abs(FRF_midpoint_column),label='Abs')
axs[0].plot(freq,FRF_midpoint_column.real,label='Re')
axs[0].plot(freq,FRF_midpoint_column.imag,label='Im')
axs[0].set_ylabel('$U$ center')
axs[1].plot(freq,abs(FRF_midpoint_beam),label='Abs')
axs[1].plot(freq,FRF_midpoint_beam.real,label='Re')
axs[1].plot(freq,FRF_midpoint_beam.imag,label='Im')
axs[1].set_ylabel('$W$ center')
axs[1].set_xlabel('f')

for ax in axs.flat:
    ax.grid()
    ax.legend(loc='best')
    ax.set_xlim([0,50])

# Step 2: Create symmetric frequency range
omega_positive = omega
omega_negative = np.linspace(-omega[-1], -omega[0], len(omega))  # Negative frequency range: -5001 to -1
omega_values = np.concatenate([omega_negative, omega_positive])  # Full symmetric frequency range

w_zero_mean = [0.0]

w1_center_freq = np.concatenate([
    w_zero_mean,
    FRF_midpoint_beam,  # Positive frequencies
    np.conj(FRF_midpoint_beam[::-1])  # Negative frequencies
])
u1_center_freq = np.concatenate([
    w_zero_mean,
    FRF_midpoint_column,  # Positive frequencies
    np.conj(FRF_midpoint_column[::-1])  # Negative frequencies
])
w1_center_time = np.fft.ifft(w1_center_freq)*0.5*N
u1_center_time = np.fft.ifft(u1_center_freq)*0.5*N
# u1_center = fft.irfft(FRF_midpoint_column) * N /2
# Frequency-domain parameters
domega = omega_values[1] - omega_values[0]  # Angular frequency step
df = domega / (2 * np.pi)                   # Frequency step
N = len(w1_center_freq)          # Length of frequency-domain data
dt = 1 / (df * N)                           # Time step

# Create time vector matching the length of w_free_frequency
t = np.linspace(0, N * dt, N)

fig, axs = plt.subplots(2,figsize=(14,10))
axs[0].plot(t, w1_center_time)
axs[1].plot(t, u1_center_time)

for ax in axs.flat:
    ax.grid()
    ax.set_xlabel('t')

axs[0].set_ylabel('$w$ center')
axs[1].set_ylabel('$u_1$ center')


    
    
