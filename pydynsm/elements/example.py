# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 21:10:58 2024

@author: Mariana
"""
# %%
import numpy as np
import sys
# from element.structuralelement import StructuralElement
sys.path.append("C:/Github/PyDynSM/pydynsm/elements")
from rod_1D_foundation import *
# from abc import ABC, abstractmethod

# %%

'''
The first step is to define the characteristics of the element to then define
the element
'''

'Define the variables (parameters) of the problem'
'Pile--------------------------------------------'
r0 = 0.5             # [m], radius of foundation pile
L_pile = 40          # [m], lenght of pile
E_pile = 3.2e10      # [Pa], E-modulus of concrete: 32GPa for C35
rho_pile = 2400      # [kg/m^3], typical density of concrete

'Soil--------------------------------------------'
E_soil = 6.89e7    # [Pa], Young's modulus
nu_soil = 0.33     # [-], Poissoin's ratio
rho_soil = 1800    # [kg/m^3], soil density
xi_soil = 0.001    # [?],  damping coefficient of soil

'Force-------------------------------------------'
F = 1              # [N], external vertical force


'Calculation of further parameters knowing the previous'
'Pile--------------------------------------------'
A = np.pi*r0**2    # area of the pile
mu = rho_pile * A                      # linear mass of the pile
c = 0                                  # damping coefficient

'Soil--------------------------------------------'
G_soil = E_soil / (2 * (1 + nu_soil))  # shear modulus of soil
c_s = np.sqrt(G_soil / rho_soil)       # shear wave velocity


# %%

# Define element using the class
omega = np.linspace(1, 10, 10)

rod_foundation = Rod1D_foundation(rho_pile, A, E_pile, L_pile, E_soil)

alpha1, alpha2 = rod_foundation.ElementWaveNumbers(omega)