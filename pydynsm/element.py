# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:34:56 2024

@author: rensv
"""

# %% Import dependencies
import numpy as np

# %% class definition
class Element:
    ne = 0

    def clear():
        Element.ne = 0
        
    def __init__ (self, nodes):
        self.nodes = nodes

        self.L = np.sqrt((nodes[1].x - nodes[0].x)**2.0 + (nodes[1].y - nodes[0].y)**2.0)

        dx = nodes[1].x - nodes[0].x
        dy = nodes[1].y - nodes[0].y

        self.cos = dx / self.L
        self.sin = dy / self.L

        R = np.zeros ((6,6))

        R[0,0] = R[1,1] = R[3,3] = R[4,4] = self.cos
        R[0,1] = R[3,4] = -self.sin
        R[1,0] = R[4,3] =  self.sin
        R[2,2] = R[5,5] = 1.0
        
        self.R  = R
        self.Rt = np.transpose(R)

        Element.ne += 1

    def set_section (self, props):
        
        if 'EA' in props:
            self.EA = props['EA']
        else:
            self.EA = 1.e20
            
        if 'ksi' in props: # MODIFIED
            self.ksi = props['ksi'] # MODIFIED
        else: # MODIFIED
            self.ksi = 0.01  # MODIFIED
            
        if 'rhoA' in props:  # MODIFIED
            self.rhoA = props['rhoA']  # MODIFIED
        else:  # MODIFIED
            self.rhoA = 1.e20  # MODIFIED
            
        if 'EI' in props:
            self.EI = props['EI']
        else:
            self.EI = 1.e20
            
        if 'omega' in props:  # MODIFIED
            self.omega = props['omega']  # MODIFIED
        else:   # MODIFIED
            self.omega = 1.e20  # MODIFIED

    def global_dofs  (self):
        return np.hstack ((self.nodes[0].dofs, self.nodes[1].dofs))

    def stiffness ( self ):
        
        k = np.zeros ((6, 6), dtype=complex) # MODIFIED
        
        ksi = self.ksi # MODIFIED
        EA = self.EA * (1 + 2j * ksi) # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        c_r = (EA/rhoA) ** 0.5  # MODIFIED
        beta_r = omega / c_r  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED

        # Extension contribution

        k[0,0] = k[3,3] = EA * beta_r * np.cos(beta_r * L) / np.sin(beta_r * L)   # MODIFIED
        k[3,0] = k[0,3] = - EA * beta_r / np.sin(beta_r * L) # MODIFIED

        # Bending contribution
        
        K_beam = np.array([[-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)]])

        k[1,1] = K_beam[0,0]
        k[1,2] = K_beam[0,1]
        k[1,4] = K_beam[0,2]
        k[1,5] = K_beam[0,3]
        k[2,1] = K_beam[1,0]
        k[2,2] = K_beam[1,1]
        k[2,4] = K_beam[1,2]
        k[2,5] = K_beam[1,3]
        k[4,1] = K_beam[2,0]
        k[4,2] = K_beam[2,1]
        k[4,4] = K_beam[2,2]
        k[4,5] = K_beam[2,3]
        k[5,1] = K_beam[3,0]
        k[5,2] = K_beam[3,1]
        k[5,4] = K_beam[3,2]
        k[5,5] = K_beam[3,3]
        
        #k[1,1] = k[4,4] =  12.0 * EI / L / L / L
        #k[1,4] = k[4,1] = -12.0 * EI / L / L / L
        #k[1,2] = k[2,1] = k[1,5] = k[5,1] = -6.0 * EI / L / L
        #k[2,4] = k[4,2] = k[4,5] = k[5,4] = 6.0 * EI / L / L
        #k[2,2] = k[5,5] = 4.0 * EI / L
        #k[2,5] = k[5,2] = 2.0 * EI / L

        return np.matmul ( np.matmul ( self.Rt, k ), self.R )
