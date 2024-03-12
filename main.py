# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:23:06 2024

@author: rensv
"""

# import python standard dependencies

import numpy as np
import matplotlib.pyplot as plt

# %% add package to sys and import 

import pydynsm as PDM

Node = PDM.Node
Constrainer = PDM.Constrainer
Element = PDM.Element

# %% Example

# %%% Parameters
EA = 7e6
EI = 1.5 * 7e06 # MODIFIED
rhoA = 1e03  # MODIFIED
q_r = 0*1e02 + 0j # MODIFIED
q_b = 1*1e06 + 0j # MODIFIED
L  = 1
omega = 100  # MODIFIED
ksi = 0.01 # MODIFIED


Node.clear()
Element.clear()

# %%% Create nodes

node1 = Node(0,0)
node2 = Node(L,0)

# %%% Create element

# initialise element
elem = Element ( [node1, node2] )

# %%% trial setting elements
elem.SetSection('Rod', {'EA': EA, 'rhoA':rhoA})

# %%% set boundary conditions

con = Constrainer()

# TRY HERE DIFFERENT BOUNDARY CONDITIONS AND CHECK WITH MAPLE FILE FOR CORRECTNESS OF RESULTS:

con.fix_dof (node1,0)
con.fix_dof (node1,1)
#con.fix_dof (node1,2)
con.fix_dof (node2,0)
con.fix_dof (node2,1)
#con.fix_dof (node2,2)

elem.add_distributed_load([q_r,q_b])

# %%%

