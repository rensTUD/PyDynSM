# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:23:06 2024

@author: rensv
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# import sys
# import os

# # Import the pydynsm module

# ## find the module in the parent directory (DELETE THESE LINES AFTER WHEN THE NOTEBOOK IS READY, DIRECTLY IMPORT THE MODULE)

# ### Get the parent directory
# parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
# ### Add the parent directory to the system path
# sys.path.append(parent_dir)

## Now you can import the entire pydynsm module
import pydynsm as PDM

Assembler = PDM.Assembler

# %% Example

# %%% Initialise an assembler with your project name

Assembler = PDM.Assembler  #Import the Assembler class from the PDM module

s1 = Assembler('1D Rod',analysis_type= 'new')

# %%% Parameters
E = 210e9
A = 1/3000
rho = 3e6
F_0 = 1e06 + 0j 
L  = 1
Omega = 100  
ksi = 0.01 

# %%% Create nodes from the assembler

node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,0)

# %%% print the coords of 2 nodes here
print(f'The coordinate of node 1 is {node1.get_coords()}')
print(f'The coordinate of node 2 is {node2.get_coords()}')

# %%% create elements
elem = s1.CreateElement([node1, node2])

# %%% plot structural elements
s1.PlotStructure(plot_elements=True)
node1.fix_node('x','z','phi_y')
node2.fix_node('z','phi_y')
elem.SetSection('Rod', {'E': E, 'rho':rho, 'A':A})
# F_omega = lambda omega: F_0 if omega == Omega else 0
q_r = lambda omega: F_0 if omega == Omega else 0
# Generate the range of omega values
# omega_values = np.linspace(0, 200, 41)
# F_values_lambda = [F_omega(omega) for omega in omega_values]

#%%% Plot the result using the lambda function
# plt.figure(figsize=(8, 4))
# plt.stem(omega_values, F_values_lambda, basefmt=" ")
# plt.xlabel('Frequency ($\omega$)')
# plt.ylabel('F($\omega$)')
# plt.title(r'Plot of $\tilde{F}(\omega)$')
# plt.grid(True)
# plt.show()

#%% test 1: free end disp @ omega = 100 rad/s need being revisited in maple
#%%% add nodal load
elem.AddDistributedLoad(x=q_r)

#%%% num_dof
dof_indices, num_dof = s1.get_dofs_elements()

# %%% run the connectivity
dof_indices, B = s1.get_B_matrix()

# %%%%

s1.run_connectivity()

# %%% Get the global stiffness and force matrices

K_global = s1.GlobalStiffness(Omega)
F_global = s1.GlobalForce(Omega)

# %%% Constrain matrix
Kc_global = s1.GlobalConstrainedStiffness(Omega)
Fc_global = s1.GlobalConstrainedForce(Omega)

# %%% solve for free end disp
u_free = s1.SolveUfree(Kc_global, Fc_global)
print(u_free.real)

#%% end of test 1
#%% test 2: free end disp FRF
# %%% create FRF inputs
omega_max = 2000
omega_values_check_test = np.linspace(1, omega_max, omega_max+1)

# %%% initialize the FRF list
FRF_U2 = []
FRF_H = []

# %%% assembler
s2 = Assembler('1D Rod', analysis_type='new')
# %%% create nodes and elements
node3 = s2.CreateNode(0,0)
node4 = s2.CreateNode(L,0)
elem2 = s2.CreateElement([node3, node4])
# %%% add constraints
node3.fix_node('x','z', 'phi_y')
node4.fix_node('z', 'phi_y')
elem2.SetSection('Rod', {'E': E, 'A': A, 'rho': rho, 'ksi': ksi})
s2.run_connectivity()
# %%% add transverse distributed load with all frequencies, need to be checked later, seems unphysical...
elem2.AddDistributedLoad(x=F_0)
# %%% main loop for generating FRF
for omega_test in omega_values_check_test:

    K_global_test = s2.GlobalStiffness(omega_test)
    F_global_test = s2.GlobalForce(omega_test)

    Kc_global_test = s2.GlobalConstrainedStiffness(omega_test)
    Fc_global_test = s2.GlobalConstrainedForce(omega_test)

    u_free = s2.SolveUfree(Kc_global_test, Fc_global_test) # Free end solution
    H_fix =  s2.SupportReactions(K_global_test,u_free,F_global_test) # Reaction forces at the fixed end

    FRF_U2.append(u_free)
    FRF_H.append(H_fix[0])

#%%% Calculate natural frequencies

n_index = np.arange(1, 5)
natural_frequencies = []

for n in n_index:
    E_complex = E * (1 + 2 * 1j * ksi)
    c_b = np.sqrt(E_complex / rho)
    omega_n = (2 * n - 1) * np.pi / (2 * L) * c_b
    natural_frequencies.append(omega_n)

#%%% Calculate limiting case under static loading
u_static = F_0 * L**2 / (2*E * A)

#%%% Plot displacement FRF
plt.figure(figsize=(8, 5))
plt.plot(omega_values_check_test, [u.real for u in FRF_U2], 'k', label='$\mathrm{Re}(u)$')
plt.plot(omega_values_check_test, [u.imag for u in FRF_U2], 'r--', label='$\mathrm{Im}(u)$')

# Add natural frequency lines with only one legend entry
for i, omega_n in enumerate(natural_frequencies):
    if i == 0:
        plt.axvline(x=omega_n.real, color='g', linestyle='--', linewidth=1, label='Natural Frequencies')
    else:
        plt.axvline(x=omega_n.real, color='g', linestyle='--', linewidth=1)

#%%% Static displacement line
plt.axhline(y=u_static.real, color='b', linestyle='--', linewidth=1, label='Static Displacement')
plt.xlabel('Frequency ($\omega_\mathrm{test}$)')
plt.ylabel('Free End Displacement $\mathrm{Re}(u)$ and $\mathrm{Im}(u)$')
plt.title('Free End Displacement vs Loading Frequency')
plt.xlim(0, 2000)
plt.legend()
plt.grid(True)
plt.show()
#%%% Plot reaction FRF

#%%%% Plot the reaction results
plt.figure(figsize=(8, 4))
plt.plot(omega_values_check_test, [H.real for H in FRF_H], 'k',label='$\mathrm{Re}(H)$')
plt.plot(omega_values_check_test, [H.imag for H in FRF_H], 'r--',label='$\mathrm{Im}(H)$')
plt.xlabel('Frequency ($\omega_\mathrm{test}$)')
plt.ylabel('Support Reaction $\mathrm{Re}(H)$ and $\mathrm{Im}(H)$')
plt.title('Support Reaction vs Loading Frequency')


#%%%% Add natural frequency lines with only one legend entry
for i, omega_n in enumerate(natural_frequencies):
    if i == 0:
        plt.axvline(x=omega_n.real, color='g', linestyle='--', linewidth=1, label='Natural Frequencies')
    else:
        plt.axvline(x=omega_n.real, color='g', linestyle='--', linewidth=1)

#%%%% Static reaction line
plt.axhline(y=(-F_0*L).real, color='b', linestyle='--', linewidth=1, label='Static Reaction')
plt.xlabel('Frequency ($\omega_\mathrm{test}$)')
plt.xlim(0, 2000)
plt.legend()
plt.grid(True)
plt.show()

#%% end of test 2
