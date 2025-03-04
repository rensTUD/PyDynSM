
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:33:47 2025

@author: GREY
"""

# Import Python standard dependencies
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import numpy.fft as fft
import scipy.optimize as opt

# %% Add package to sys and import 
import pydynsm as PDM

# All parameters
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
s1 = Assembler('Frame', analysis_type='new')

# %%% Initialise an assembler with your project name
nodes = []
nodes.append(s1.CreateNode(0, 0))
nodes.append(s1.CreateNode(0, h))
nodes.append(s1.CreateNode(b2, h))
nodes.append(s1.CreateNode(b2, 0))

s1.PlotStructure()

elements = []
elements.append(s1.CreateElement([nodes[0], nodes[1]]))
elements.append(s1.CreateElement([nodes[1], nodes[2]]))
elements.append(s1.CreateElement([nodes[2], nodes[3]]))
# %%%
s1.PlotStructure(plot_elements=True)

# %%
nodes[0].fix_node('x', 'z', 'phi_y')
nodes[1].fix_node('z')
nodes[2].fix_node('z')
nodes[3].fix_node('x', 'z', 'phi_y')

# %% 
# Assign parameters
column = {
    'E': E,
    'A': A_C,
    'Ib': I_C,
    'ksi': 0,
    'rho': rho
}

beam = {
    'E': E,
    'A': A_G,
    'Ib': I_G,
    'ksi': 0,
    'rho': rho
}

rodb = {
    'E': E,
    'A': A_G,
    'ksi': 0,
    'rho': rho
}

rodc = {
    'E': E,
    'A': A_C,
    'ksi': 0,
    'rho': rho
}

elements[0].SetSection('EulerBernoulli Beam', column)
elements[2].SetSection('EulerBernoulli Beam', column)
elements[1].SetSection('EulerBernoulli Beam', beam)

s1.run_connectivity()

p_x = lambda omega: 2.9e3 if omega == Omega else 0

Kc_global = s1.GlobalConstrainedStiffness(Omega)
Fc_global = s1.GlobalConstrainedForce(Omega)
