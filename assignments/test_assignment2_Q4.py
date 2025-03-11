# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:59:03 2025

@author: GREY
"""

#%% import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA
import scipy.optimize as opt
import pydynsm as PDM
import scipy.integrate as spi

#%% create assembler
Assembler = PDM.Assembler
s1 = Assembler('bridge')
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

EA_b = E*A_b
EI_b = E*I_b
rhoA_b = rho*A_b
EA_c = E*A_c
rhoA_c = rho*A_c
EA_t = E*A_t
rhoA_t = rho*A_t
# %% nodes
node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(l2,0)
node3 = s1.CreateNode(2*l2,0)
node19 = s1.CreateNode(1*l2+l2/5,0)
node20 = s1.CreateNode(1*l2+l2/5+3,0)
node4 = s1.CreateNode(3*l2,0)
node5 = s1.CreateNode(4*l2,0)
node6 = s1.CreateNode(4*l2+l1,0)
node7 = s1.CreateNode(4*l2+l1,-l1,dof_config = ['x','z'])
node8 = s1.CreateNode(4*l2+l1+l3,0)
node21 = s1.CreateNode(4*l2+l1+l3+l3/5,0)
node22 = s1.CreateNode(4*l2+l1+l3+l3/5+3,0)
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
    (node2, node19),
    (node19, node20),
    (node20, node3),
    (node3, node4),
    (node4, node5),
    (node5, node6),
    (node6, node8),
    (node8, node21),
    (node21, node22),
    (node22, node9),
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

# %% constraint nodes
node1.fix_node('z')
node7.fix_node('x','z')

# %% assign sections
for beam in beams[:-1]:
    beam.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0.01, 'Wb':I_b})
    beam.SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0.01})

beams[-1].SetSection('EulerBernoulli Beam with foundation', {'E': E, 'A':A_b, 'rho':rho, 'Ib':I_b, 'ksi':0.01, 'kd':k_d,'cd':c_d})
beams[-1].SetSection('Rod', {'E': E, 'A':A_b, 'rho':rho,'ksi':0.01})
for cable in front_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':0.01})
    # cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':0})
for cable in back_cables:
    cable.SetSection('Rod', {'E': E, 'A':A_c, 'rho':rho,'ksi':0.01})
    # cable.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_c, 'rho':rho, 'Ib':I_c,'ksi':0})
for truss in trusses:
    truss.SetSection('Rod', {'E': E, 'A':A_t, 'rho':rho,'ksi':0.01})
    # truss.SetSection('EulerBernoulli Beam', {'E': E, 'A':A_t, 'rho':rho, 'Ib':I_t,'ksi':0})
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

#%% plot wind load along height
plt.figure()
plt.plot(np.array([load1, load2, load3, load4, load5])/1e3,np.array([0,l1,2*l1,3*l1,4*l1]))
plt.ylabel('Height [m]')
plt.xlabel('Equivalent Wind load [kN]')
plt.grid()
plt.show()

#%% apply wind load

omega_w = 12
w_x1 = lambda omega: load1 if omega == omega_w else 0
w_x2 = lambda omega: load2 if omega == omega_w else 0
w_x3 = lambda omega: load3 if omega == omega_w else 0
w_x4 = lambda omega: load4 if omega == omega_w else 0
w_x5 = lambda omega: load5 if omega == omega_w else 0


node18.add_load(x=w_x5)
node17.add_load(x=w_x4)
node15.add_load(x=w_x3)
node13.add_load(x=w_x2)
node5.add_load(x=w_x1)

#%% vehicle loading
omega_c = 34
w_c = (-0.6-0.5*9)*1e4
p_c = lambda omega: w_c if omega == omega_c else 0


node19.add_load(z=p_c)
node20.add_load(z=p_c)
node21.add_load(z=p_c)
node22.add_load(z=p_c)

Kc_global = s1.GlobalConstrainedStiffness(omega_w)
Fc_global = s1.GlobalConstrainedForce(omega_w)

Kc_global2 = s1.GlobalConstrainedStiffness(omega_c)
Fc_global2 = s1.GlobalConstrainedForce(omega_c)

u_free = s1.SolveUfree(Kc_global, Fc_global)
u_free2 = s1.SolveUfree(Kc_global2, Fc_global2)

u_elem = s1.FullDisplacement(u_free)
u_elem2 = s1.FullDisplacement(u_free2)
disp_wind = s1.ElementDisplacements(u_elem,omega_w,400)
disp_vehicle = s1.ElementDisplacements(u_elem2,omega_c,400)
force_wind = s1.ElementForces(u_elem,omega_w,400)
force_vehicle = s1.ElementForces(u_elem2,omega_c,400)

num_points = 400
local_coords_girder = np.zeros((len(beams),num_points))
local_w_girder = np.zeros((len(beams),num_points),complex)
local_u_girder = np.zeros((len(beams),num_points),complex)
local_phi_girder = np.zeros((len(beams),num_points),complex)
local_N_girder = np.zeros((len(beams),num_points),complex)
local_M_girder = np.zeros((len(beams),num_points),complex)
local_V_girder = np.zeros((len(beams),num_points),complex)

for ii, beam in enumerate(beams):
    local_coords_girder[ii,:] = beams[ii].L*np.linspace(0,1,num_points)+local_coords_girder[ii-1,-1]  
    local_u_girder[ii,:] = disp_wind[beam.id][0]
    local_w_girder[ii,:] = disp_wind[beam.id][1]
    local_phi_girder[ii,:] = disp_wind[beam.id][-1]
    local_N_girder[ii,:] = force_wind[beam.id][0]
    local_V_girder[ii,:] = force_wind[beam.id][1]
    local_M_girder[ii,:] = force_wind[beam.id][-1]
    
local_coords_girder = local_coords_girder.ravel()

u_girder_wind = local_u_girder.ravel()
w_girder_wind = local_w_girder.ravel()
N_girder_wind = local_N_girder.ravel()
phi_girder_wind = local_phi_girder.ravel()
V_girder_wind = local_V_girder.ravel()
M_girder_wind = local_M_girder.ravel()

local_w_girder = np.zeros((len(beams),num_points),complex)
local_u_girder = np.zeros((len(beams),num_points),complex)
local_phi_girder = np.zeros((len(beams),num_points),complex)
local_N_girder = np.zeros((len(beams),num_points),complex)
local_M_girder = np.zeros((len(beams),num_points),complex)
local_V_girder = np.zeros((len(beams),num_points),complex)

for ii, beam in enumerate(beams):
    local_u_girder[ii,:] = disp_vehicle[beam.id][0]
    local_w_girder[ii,:] = disp_vehicle[beam.id][1]
    local_phi_girder[ii,:] = disp_vehicle[beam.id][-1]
    local_N_girder[ii,:] = force_vehicle[beam.id][0]
    local_V_girder[ii,:] = force_vehicle[beam.id][1]
    local_M_girder[ii,:] = force_vehicle[beam.id][-1]


u_girder_car = local_u_girder.ravel()
w_girder_car = local_w_girder.ravel()
phi_girder_car = local_phi_girder.ravel()
N_girder_car = local_N_girder.ravel()
V_girder_car = local_V_girder.ravel()
M_girder_car = local_M_girder.ravel()



post_results = [
                u_girder_car + u_girder_wind,
                w_girder_car + w_girder_wind,
                phi_girder_car + phi_girder_wind,
                N_girder_car + N_girder_wind,
                V_girder_car + V_girder_wind,
                M_girder_car + M_girder_wind
                ]

axis_name = ['u [m]','w [m]','phi [rad]','N [N]','V [N]','M [Nm]']

fig, axs = plt.subplots(len(axis_name), 1, sharex=True, figsize=(14,10))

for ii, result in enumerate(post_results):
    axs[ii].plot(local_coords_girder,result.real,'k',label='Re')
    axs[ii].plot(local_coords_girder,result.imag,'r',label='Im')
    axs[ii].set_ylabel(axis_name[ii])
    
axs[-1].set_xlabel('Length[m]')

for ax in axs:
    ax.legend()
    ax.grid()
plt.show()
