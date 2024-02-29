# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:34:56 2024

@author: rensv
"""

# %% Import dependencies
import numpy as np

# %% class definition
class Constrainer:
    def __init__ (self):
        self.dofs = []

    def fix_dof (self, node, dof):
        self.dofs.append (node.dofs[dof])
 
    def fix_node (self, node):
        for dof in node.dofs:
            self.dofs.append (dof)       

    def constrain (self, k, f):
        kc = np.copy (k)
        fc = np.copy (f)
        
        for dof in self.dofs:
            fc[dof] = 0.0
            kc[:,dof] = kc[dof,:] = 0.0
            kc[dof,dof]           = 1.0

        return kc, fc
