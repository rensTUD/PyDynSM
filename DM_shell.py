#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:52:43 2024

@author: lukas
"""

# a function that returns the root based on the input parameters

import numpy as np
from .structuralelement import StructuralElement, ElementFactory


@ElementFactory.ElementType('DM Shell')
class DM_Shell(StructuralElement):

    def __init__(self, rho, E, h, R, L, n, nu = None, ksi = None):
        
        # define what dofs the eb beam contributes to and initialise
        # dofs = [1,2]
        # super().__init__(dofs)
        
        # Initialise lcoal rod element with necessary parameters
        self.ksi = ksi if ksi is not None else 0.01 # assisgn ksi if given otherwise assign a default value (think of changing this)
        self.nu = nu if nu is not None else 0.3 # assisgn poisson ratio if given otherwise assign a default value (think of changing this)
        self.rho = rho
        self.E= E * (1+2j*self.ksi)
        self.h= h
        self.R= R
        self.L = L
        self.n = n
        
        
    def LocalStiffness(self, omega, s = True):
        '''
        Determines the stiffness of the rod. As its 1D, stiffness will be 2x2 matrix as:
            [V_left, M_left, V_right, M_right] = K.[w_left, phi_left, w_right, phi_right]
        
        where:
            K = [K_V_ll, K_M_ll, K_V_lr, K_M_lr;
                 K_V_rl, K_M_rl, K_V_rr, K_M_rr]
        '''    

        k = self.ElementWaveNumbers(omega, s)
        gamma_u, gamma_v, gamma_w = self.mode_shape(omega, s)
        
        rho = self.rho
        nu = self.nu
        E = self.E
        h = self.h
        R = self.R
        L = self.L
        n = self.n
        
        
        if s == True:         # computes symmetric response
            if self.n == 0:
        
                k1 = k[0]
                k2 = k[1]
                k3 = k[2]
                k4 = k[3]
                k5 = k[4]
                k6 = k[5]
                
                gamma_u1 = gamma_u[0]
                gamma_u2 = gamma_u[1]
                gamma_u3 = gamma_u[2]
                gamma_u4 = gamma_u[3]
                gamma_u5 = gamma_u[4]
                gamma_u6 = gamma_u[5]
                
                A_mat = np.array([[gamma_u1,gamma_u2,gamma_u3,gamma_u4,gamma_u5,gamma_u6],
                           [1,1,1,1,1,1],
                           [complex(0, -1) * k1,complex(0, -1) * k2,complex(0, -1) * k3,complex(0, -1) * k4,complex(0, -1) * k5,complex(0, -1) * k6],
                           [np.exp(complex(0, -1) * k1 * L) * gamma_u1,np.exp(complex(0, -1) * k2 * L) * gamma_u2,np.exp(complex(0, -1) * k3 * L) * gamma_u3,np.exp(complex(0, -1) * k4 * L) * gamma_u4,np.exp(complex(0, -1) * k5 * L) * gamma_u5,np.exp(complex(0, -1) * k6 * L) * gamma_u6],
                           [np.exp(complex(0, -1) * k1 * L),np.exp(complex(0, -1) * k2 * L),np.exp(complex(0, -1) * k3 * L),np.exp(complex(0, -1) * k4 * L),np.exp(complex(0, -1) * k5 * L),np.exp(complex(0, -1) * k6 * L)],
                           [complex(0, -1) * np.exp(complex(0, -1) * k1 * L) * k1,complex(0, -1) * np.exp(complex(0, -1) * k2 * L) * k2,complex(0, -1) * np.exp(complex(0, -1) * k3 * L) * k3,complex(0, -1) * np.exp(complex(0, -1) * k4 * L) * k4,complex(0, -1) * np.exp(complex(0, -1) * k5 * L) * k5,complex(0, -1) * np.exp(complex(0, -1) * k6 * L) * k6]])
                
                B_mat = np.array([[(complex(0, 1) * R * k1 * gamma_u1 - nu) * h * E / (R * nu ** 2 - R),(complex(0, 1) * R * k2 * gamma_u2 - nu) * h * E / (R * nu ** 2 - R),(complex(0, 1) * R * k3 * gamma_u3 - nu) * h * E / (R * nu ** 2 - R),(complex(0, 1) * R * k4 * gamma_u4 - nu) * h * E / (R * nu ** 2 - R),(complex(0, 1) * R * k5 * gamma_u5 - nu) * h * E / (R * nu ** 2 - R),(complex(0, 1) * R * k6 * gamma_u6 - nu) * h * E / (R * nu ** 2 - R)],
                           [complex(0, -1) * k1 * E * (-R ** 2 * k1 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * h ** 3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * k2 * E * (-R ** 2 * k2 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * h ** 3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * k3 * E * (-R ** 2 * k3 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * h ** 3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * k4 * E * (-R ** 2 * k4 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * h ** 3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * k5 * E * (-R ** 2 * k5 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * h ** 3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * k6 * E * (-R ** 2 * k6 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * h ** 3 / (12 * nu ** 2 - 12) / R ** 2],
                           [-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k1 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k2 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k3 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k4 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k5 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k6 ** 2 + n ** 2 * nu) / R ** 2],
                           [-E * h * np.exp(complex(0, -1) * k1 * L) * (complex(0, 1) * R * k1 * gamma_u1 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k2 * L) * (complex(0, 1) * R * k2 * gamma_u2 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k3 * L) * (complex(0, 1) * R * k3 * gamma_u3 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k4 * L) * (complex(0, 1) * R * k4 * gamma_u4 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k5 * L) * (complex(0, 1) * R * k5 * gamma_u5 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k6 * L) * (complex(0, 1) * R * k6 * gamma_u6 - nu) / (nu ** 2 - 1) / R],
                           [complex(0, 1) * (-R ** 2 * k1 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k1 * L) * h ** 3 * E * k1 / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * (-R ** 2 * k2 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k2 * L) * h ** 3 * E * k2 / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * (-R ** 2 * k3 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k3 * L) * h ** 3 * E * k3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * (-R ** 2 * k4 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k4 * L) * h ** 3 * E * k4 / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * (-R ** 2 * k5 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k5 * L) * h ** 3 * E * k5 / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * (-R ** 2 * k6 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k6 * L) * h ** 3 * E * k6 / (12 * nu ** 2 - 12) / R ** 2],
                           [E * h ** 3 * np.exp(complex(0, -1) * k1 * L) * (R ** 2 * k1 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k2 * L) * (R ** 2 * k2 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k3 * L) * (R ** 2 * k3 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k4 * L) * (R ** 2 * k4 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k5 * L) * (R ** 2 * k5 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k6 * L) * (R ** 2 * k6 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2]])
                
                # H_n = np.array([[np.cos(n * theta),0,0,0,0,0],
                #             [0,np.cos(n * theta),0,0,0,0],
                #             [0,0,np.cos(n * theta),0,0,0],
                #             [0,0,0,np.cos(n * theta),0,0],
                #             [0,0,0,0,np.cos(n * theta),0],
                #             [0,0,0,0,0,np.cos(n * theta)]])
        
                
            else:
                k1 = k[0]
                k2 = k[1]
                k3 = k[2]
                k4 = k[3]
                k5 = k[4]
                k6 = k[5]
                k7 = k[6]
                k8 = k[7]
                
                gamma_u1 = gamma_u[0]
                gamma_u2 = gamma_u[1]
                gamma_u3 = gamma_u[2]
                gamma_u4 = gamma_u[3]
                gamma_u5 = gamma_u[4]
                gamma_u6 = gamma_u[5]
                gamma_u7 = gamma_u[6]
                gamma_u8 = gamma_u[7]
                
                gamma_v1 = gamma_v[0]
                gamma_v2 = gamma_v[1]
                gamma_v3 = gamma_v[2]
                gamma_v4 = gamma_v[3]
                gamma_v5 = gamma_v[4]
                gamma_v6 = gamma_v[5]
                gamma_v7 = gamma_v[6]
                gamma_v8 = gamma_v[7]
                
                
                A_mat = np.array([[gamma_u1,gamma_u2,gamma_u3,gamma_u4,gamma_u5,gamma_u6,gamma_u7,gamma_u8],
                           [gamma_v1,gamma_v2,gamma_v3,gamma_v4,gamma_v5,gamma_v6,gamma_v7,gamma_v8],
                           [1,1,1,1,1,1,1,1],
                           [complex(0, -1) * k1,complex(0, -1) * k2,complex(0, -1) * k3,complex(0, -1) * k4,complex(0, -1) * k5,complex(0, -1) * k6,complex(0, -1) * k7,complex(0, -1) * k8],
                           [np.exp(complex(0, -1) * k1 * L) * gamma_u1,np.exp(complex(0, -1) * k2 * L) * gamma_u2,np.exp(complex(0, -1) * k3 * L) * gamma_u3,np.exp(complex(0, -1) * k4 * L) * gamma_u4,np.exp(complex(0, -1) * k5 * L) * gamma_u5,np.exp(complex(0, -1) * k6 * L) * gamma_u6,np.exp(complex(0, -1) * k7 * L) * gamma_u7,np.exp(complex(0, -1) * k8 * L) * gamma_u8],
                           [np.exp(complex(0, -1) * k1 * L) * gamma_v1,np.exp(complex(0, -1) * k2 * L) * gamma_v2,np.exp(complex(0, -1) * k3 * L) * gamma_v3,np.exp(complex(0, -1) * k4 * L) * gamma_v4,np.exp(complex(0, -1) * k5 * L) * gamma_v5,np.exp(complex(0, -1) * k6 * L) * gamma_v6,np.exp(complex(0, -1) * k7 * L) * gamma_v7,np.exp(complex(0, -1) * k8 * L) * gamma_v8],
                           [np.exp(complex(0, -1) * k1 * L),np.exp(complex(0, -1) * k2 * L),np.exp(complex(0, -1) * k3 * L),np.exp(complex(0, -1) * k4 * L),np.exp(complex(0, -1) * k5 * L),np.exp(complex(0, -1) * k6 * L),np.exp(complex(0, -1) * k7 * L),np.exp(complex(0, -1) * k8 * L)],
                           [complex(0, -1) * np.exp(complex(0, -1) * k1 * L) * k1,complex(0, -1) * np.exp(complex(0, -1) * k2 * L) * k2,complex(0, -1) * np.exp(complex(0, -1) * k3 * L) * k3,complex(0, -1) * np.exp(complex(0, -1) * k4 * L) * k4,complex(0, -1) * np.exp(complex(0, -1) * k5 * L) * k5,complex(0, -1) * np.exp(complex(0, -1) * k6 * L) * k6,complex(0, -1) * np.exp(complex(0, -1) * k7 * L) * k7,complex(0, -1) * np.exp(complex(0, -1) * k8 * L) * k8]])
                
                B_mat = np.array([[0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * gamma_u1 * k1 * R - n * nu * gamma_v1 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k2 * gamma_u2 * R - n * nu * gamma_v2 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k3 * gamma_u3 * R - n * nu * gamma_v3 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k4 * gamma_u4 * R - n * nu * gamma_v4 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k5 * gamma_u5 * R - n * nu * gamma_v5 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k6 * gamma_u6 * R - n * nu * gamma_v6 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k7 * gamma_u7 * R - n * nu * gamma_v7 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * gamma_u8 * k8 * R - n * nu * gamma_v8 - nu) * h * E],
                          [-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k1 * gamma_v1 + complex(0, 1) * h ** 2 * k1 * n + 6 * R * n * gamma_u1) * h * E / 12,-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k2 * gamma_v2 + complex(0, 1) * h ** 2 * k2 * n + 6 * R * n * gamma_u2) * h * E / 12,-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k3 * gamma_v3 + complex(0, 1) * h ** 2 * k3 * n + 6 * R * n * gamma_u3) * h * E / 12,-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k4 * gamma_v4 + complex(0, 1) * h ** 2 * k4 * n + 6 * R * n * gamma_u4) * h * E / 12,-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k5 * gamma_v5 + complex(0, 1) * h ** 2 * k5 * n + 6 * R * n * gamma_u5) * h * E / 12,-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k6 * gamma_v6 + complex(0, 1) * h ** 2 * k6 * n + 6 * R * n * gamma_u6) * h * E / 12,-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k7 * gamma_v7 + complex(0, 1) * h ** 2 * k7 * n + 6 * R * n * gamma_u7) * h * E / 12,-0.1e1 / R ** 2 / (1 + nu) * (complex(0, 6) * R ** 2 * k8 * gamma_v8 + complex(0, 1) * h ** 2 * k8 * n + 6 * R * n * gamma_u8) * h * E / 12],
                          [complex(0, -1) * (-R ** 2 * k1 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k1 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * (-R ** 2 * k2 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k2 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * (-R ** 2 * k3 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * (-R ** 2 * k4 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k4 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * (-R ** 2 * k5 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k5 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * (-R ** 2 * k6 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k6 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * (-R ** 2 * k7 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k7 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * (-R ** 2 * k8 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * E * h ** 3 * k8 / (12 * nu ** 2 - 12) / R ** 2],
                          [-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k1 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k2 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k3 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k4 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k5 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k6 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k7 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k8 ** 2 + n ** 2 * nu) / R ** 2],
                          [-E * h * np.exp(complex(0, -1) * k1 * L) * (complex(0, 1) * gamma_u1 * k1 * R - n * nu * gamma_v1 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k2 * L) * (complex(0, 1) * k2 * gamma_u2 * R - n * nu * gamma_v2 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k3 * L) * (complex(0, 1) * k3 * gamma_u3 * R - n * nu * gamma_v3 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k4 * L) * (complex(0, 1) * k4 * gamma_u4 * R - n * nu * gamma_v4 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k5 * L) * (complex(0, 1) * k5 * gamma_u5 * R - n * nu * gamma_v5 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k6 * L) * (complex(0, 1) * k6 * gamma_u6 * R - n * nu * gamma_v6 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k7 * L) * (complex(0, 1) * k7 * gamma_u7 * R - n * nu * gamma_v7 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k8 * L) * (complex(0, 1) * gamma_u8 * k8 * R - n * nu * gamma_v8 - nu) / (nu ** 2 - 1) / R],
                          [E * h * np.exp(complex(0, -1) * k1 * L) * (complex(0, 6) * R ** 2 * k1 * gamma_v1 + complex(0, 1) * h ** 2 * k1 * n + 6 * R * n * gamma_u1) / (1 + nu) / R ** 2 / 12,E * h * np.exp(complex(0, -1) * k2 * L) * (complex(0, 6) * R ** 2 * k2 * gamma_v2 + complex(0, 1) * h ** 2 * k2 * n + 6 * R * n * gamma_u2) / (1 + nu) / R ** 2 / 12,E * h * np.exp(complex(0, -1) * k3 * L) * (complex(0, 6) * R ** 2 * k3 * gamma_v3 + complex(0, 1) * h ** 2 * k3 * n + 6 * R * n * gamma_u3) / (1 + nu) / R ** 2 / 12,E * h * np.exp(complex(0, -1) * k4 * L) * (complex(0, 6) * R ** 2 * k4 * gamma_v4 + complex(0, 1) * h ** 2 * k4 * n + 6 * R * n * gamma_u4) / (1 + nu) / R ** 2 / 12,E * h * np.exp(complex(0, -1) * k5 * L) * (complex(0, 6) * R ** 2 * k5 * gamma_v5 + complex(0, 1) * h ** 2 * k5 * n + 6 * R * n * gamma_u5) / (1 + nu) / R ** 2 / 12,E * h * np.exp(complex(0, -1) * k6 * L) * (complex(0, 6) * R ** 2 * k6 * gamma_v6 + complex(0, 1) * h ** 2 * k6 * n + 6 * R * n * gamma_u6) / (1 + nu) / R ** 2 / 12,E * h * np.exp(complex(0, -1) * k7 * L) * (complex(0, 6) * R ** 2 * k7 * gamma_v7 + complex(0, 1) * h ** 2 * k7 * n + 6 * R * n * gamma_u7) / (1 + nu) / R ** 2 / 12,E * h * np.exp(complex(0, -1) * k8 * L) * (complex(0, 6) * R ** 2 * k8 * gamma_v8 + complex(0, 1) * h ** 2 * k8 * n + 6 * R * n * gamma_u8) / (1 + nu) / R ** 2 / 12],
                          [complex(0, 1) * E * k1 * h ** 3 * (-R ** 2 * k1 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k1 * L) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * E * k2 * h ** 3 * (-R ** 2 * k2 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k2 * L) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * E * k3 * h ** 3 * (-R ** 2 * k3 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k3 * L) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * E * k4 * h ** 3 * (-R ** 2 * k4 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k4 * L) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * E * k5 * h ** 3 * (-R ** 2 * k5 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k5 * L) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * E * k6 * h ** 3 * (-R ** 2 * k6 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k6 * L) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * E * k7 * h ** 3 * (-R ** 2 * k7 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k7 * L) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * E * k8 * h ** 3 * (-R ** 2 * k8 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * np.exp(complex(0, -1) * k8 * L) / (12 * nu ** 2 - 12) / R ** 2],
                          [E * h ** 3 * np.exp(complex(0, -1) * k1 * L) * (R ** 2 * k1 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k2 * L) * (R ** 2 * k2 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k3 * L) * (R ** 2 * k3 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k4 * L) * (R ** 2 * k4 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k5 * L) * (R ** 2 * k5 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k6 * L) * (R ** 2 * k6 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k7 * L) * (R ** 2 * k7 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k8 * L) * (R ** 2 * k8 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2]])

                # H_n = np.array([[np.cos(n * theta),0,0,0,0,0,0,0],
                #            [0,np.sin(n * theta),0,0,0,0,0,0],
                #            [0,0,np.cos(n * theta),0,0,0,0,0],
                #            [0,0,0,np.cos(n * theta),0,0,0,0],
                #            [0,0,0,0,np.cos(n * theta),0,0,0],
                #            [0,0,0,0,0,np.sin(n * theta),0,0],
                #            [0,0,0,0,0,0,np.cos(n * theta),0],
                #            [0,0,0,0,0,0,0,np.cos(n * theta)]])
            K_local = np.matmul(B_mat, np.linalg.inv(A_mat))
                
        else:           # computes antisymmetric response
            if self.n == 0:
                k1 = k[0]
                k2 = k[1]
                
                K_local = np.array([[complex(0, -0.1e1 / 0.2e1) * h * E * (k1 * np.exp(complex(0, -1) * k2 * L) - k2 * np.exp(complex(0, -1) * k1 * L)) / (1 + nu) / (np.exp(complex(0, -1) * k2 * L) - np.exp(complex(0, -1) * k1 * L)),complex(0, 0.1e1 / 0.2e1) * h * E * (k1 - k2) / (1 + nu) / (np.exp(complex(0, -1) * k2 * L) - np.exp(complex(0, -1) * k1 * L))],
                           [complex(0, 0.1e1 / 0.2e1) * h * E * (k1 - k2) / (1 + nu) / (np.exp(complex(0, -1) * k2 * L) - np.exp(complex(0, -1) * k1 * L)) * np.exp(complex(0, -1) * L * (k1 + k2)),complex(0, 0.1e1 / 0.2e1) * E * h * (k2 * np.exp(complex(0, -1) * k2 * L) - np.exp(complex(0, -1) * k1 * L) * k1) / (1 + nu) / (np.exp(complex(0, -1) * k2 * L) - np.exp(complex(0, -1) * k1 * L))]])
                # H_n = np.array([[np.sin(n * theta),0],
                #                 [0,np.sin(n * theta)]])
                
            else:
                k1 = k[0]
                k2 = k[1]
                k3 = k[2]
                k4 = k[3]
                k5 = k[4]
                k6 = k[5]
                k7 = k[6]
                k8 = k[7]
                
                gamma_u1 = gamma_u[0]
                gamma_u2 = gamma_u[1]
                gamma_u3 = gamma_u[2]
                gamma_u4 = gamma_u[3]
                gamma_u5 = gamma_u[4]
                gamma_u6 = gamma_u[5]
                gamma_u7 = gamma_u[6]
                gamma_u8 = gamma_u[7]
                
                gamma_v1 = gamma_v[0]
                gamma_v2 = gamma_v[1]
                gamma_v3 = gamma_v[2]
                gamma_v4 = gamma_v[3]
                gamma_v5 = gamma_v[4]
                gamma_v6 = gamma_v[5]
                gamma_v7 = gamma_v[6]
                gamma_v8 = gamma_v[7]
                
                A_mat = np.array([[gamma_u1,gamma_u2,gamma_u3,gamma_u4,gamma_u5,gamma_u6,gamma_u7,gamma_u8],
                           [gamma_v1,gamma_v2,gamma_v3,gamma_v4,gamma_v5,gamma_v6,gamma_v7,gamma_v8],
                           [1,1,1,1,1,1,1,1],
                           [complex(0, -1) * k1,complex(0, -1) * k2,complex(0, -1) * k3,complex(0, -1) * k4,complex(0, -1) * k5,complex(0, -1) * k6,complex(0, -1) * k7,complex(0, -1) * k8],
                           [np.exp(complex(0, -1) * k1 * L) * gamma_u1,np.exp(complex(0, -1) * k2 * L) * gamma_u2,np.exp(complex(0, -1) * k3 * L) * gamma_u3,np.exp(complex(0, -1) * k4 * L) * gamma_u4,np.exp(complex(0, -1) * k5 * L) * gamma_u5,np.exp(complex(0, -1) * k6 * L) * gamma_u6,np.exp(complex(0, -1) * k7 * L) * gamma_u7,np.exp(complex(0, -1) * k8 * L) * gamma_u8],
                           [np.exp(complex(0, -1) * k1 * L) * gamma_v1,np.exp(complex(0, -1) * k2 * L) * gamma_v2,np.exp(complex(0, -1) * k3 * L) * gamma_v3,np.exp(complex(0, -1) * k4 * L) * gamma_v4,np.exp(complex(0, -1) * k5 * L) * gamma_v5,np.exp(complex(0, -1) * k6 * L) * gamma_v6,np.exp(complex(0, -1) * k7 * L) * gamma_v7,np.exp(complex(0, -1) * k8 * L) * gamma_v8],
                           [np.exp(complex(0, -1) * k1 * L),np.exp(complex(0, -1) * k2 * L),np.exp(complex(0, -1) * k3 * L),np.exp(complex(0, -1) * k4 * L),np.exp(complex(0, -1) * k5 * L),np.exp(complex(0, -1) * k6 * L),np.exp(complex(0, -1) * k7 * L),np.exp(complex(0, -1) * k8 * L)],
                           [complex(0, -1) * np.exp(complex(0, -1) * k1 * L) * k1,complex(0, -1) * np.exp(complex(0, -1) * k2 * L) * k2,complex(0, -1) * np.exp(complex(0, -1) * k3 * L) * k3,complex(0, -1) * np.exp(complex(0, -1) * k4 * L) * k4,complex(0, -1) * np.exp(complex(0, -1) * k5 * L) * k5,complex(0, -1) * np.exp(complex(0, -1) * k6 * L) * k6,complex(0, -1) * np.exp(complex(0, -1) * k7 * L) * k7,complex(0, -1) * np.exp(complex(0, -1) * k8 * L) * k8]])
                
                B_mat = np.array([[0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * gamma_u1 * k1 * R + n * nu * gamma_v1 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k2 * gamma_u2 * R + n * nu * gamma_v2 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k3 * gamma_u3 * R + n * nu * gamma_v3 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k4 * gamma_u4 * R + n * nu * gamma_v4 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k5 * gamma_u5 * R + n * nu * gamma_v5 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k6 * gamma_u6 * R + n * nu * gamma_v6 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * k7 * gamma_u7 * R + n * nu * gamma_v7 - nu) * h * E,0.1e1 / R / (nu ** 2 - 1) * (complex(0, 1) * gamma_u8 * k8 * R + n * nu * gamma_v8 - nu) * h * E],
                           [(complex(0, -6) * R ** 2 * k1 * gamma_v1 + complex(0, 1) * h ** 2 * k1 * n + 6 * R * n * gamma_u1) * h * E / R ** 2 / (1 + nu) / 12,(complex(0, -6) * R ** 2 * k2 * gamma_v2 + complex(0, 1) * h ** 2 * k2 * n + 6 * R * n * gamma_u2) * h * E / R ** 2 / (1 + nu) / 12,(complex(0, -6) * R ** 2 * k3 * gamma_v3 + complex(0, 1) * h ** 2 * k3 * n + 6 * R * n * gamma_u3) * h * E / R ** 2 / (1 + nu) / 12,(complex(0, -6) * R ** 2 * k4 * gamma_v4 + complex(0, 1) * h ** 2 * k4 * n + 6 * R * n * gamma_u4) * h * E / R ** 2 / (1 + nu) / 12,(complex(0, -6) * R ** 2 * k5 * gamma_v5 + complex(0, 1) * h ** 2 * k5 * n + 6 * R * n * gamma_u5) * h * E / R ** 2 / (1 + nu) / 12,(complex(0, -6) * R ** 2 * k6 * gamma_v6 + complex(0, 1) * h ** 2 * k6 * n + 6 * R * n * gamma_u6) * h * E / R ** 2 / (1 + nu) / 12,(complex(0, -6) * R ** 2 * k7 * gamma_v7 + complex(0, 1) * h ** 2 * k7 * n + 6 * R * n * gamma_u7) * h * E / R ** 2 / (1 + nu) / 12,(complex(0, -6) * R ** 2 * k8 * gamma_v8 + complex(0, 1) * h ** 2 * k8 * n + 6 * R * n * gamma_u8) * h * E / R ** 2 / (1 + nu) / 12],
                           [complex(0, -1) * E * h ** 3 * (-R ** 2 * k1 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k1 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * E * h ** 3 * (-R ** 2 * k2 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k2 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * E * h ** 3 * (-R ** 2 * k3 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k3 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * E * h ** 3 * (-R ** 2 * k4 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k4 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * E * h ** 3 * (-R ** 2 * k5 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k5 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * E * h ** 3 * (-R ** 2 * k6 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k6 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * E * h ** 3 * (-R ** 2 * k7 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k7 / (12 * nu ** 2 - 12) / R ** 2,complex(0, -1) * E * h ** 3 * (-R ** 2 * k8 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) * k8 / (12 * nu ** 2 - 12) / R ** 2],
                           [-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k1 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k2 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k3 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k4 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k5 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k6 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k7 ** 2 + n ** 2 * nu) / R ** 2,-E * h ** 3 / (12 * nu ** 2 - 12) * (R ** 2 * k8 ** 2 + n ** 2 * nu) / R ** 2],
                           [-E * h * np.exp(complex(0, -1) * k1 * L) * (complex(0, 1) * gamma_u1 * k1 * R + n * nu * gamma_v1 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k2 * L) * (complex(0, 1) * k2 * gamma_u2 * R + n * nu * gamma_v2 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k3 * L) * (complex(0, 1) * k3 * gamma_u3 * R + n * nu * gamma_v3 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k4 * L) * (complex(0, 1) * k4 * gamma_u4 * R + n * nu * gamma_v4 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k5 * L) * (complex(0, 1) * k5 * gamma_u5 * R + n * nu * gamma_v5 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k6 * L) * (complex(0, 1) * k6 * gamma_u6 * R + n * nu * gamma_v6 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k7 * L) * (complex(0, 1) * k7 * gamma_u7 * R + n * nu * gamma_v7 - nu) / (nu ** 2 - 1) / R,-E * h * np.exp(complex(0, -1) * k8 * L) * (complex(0, 1) * gamma_u8 * k8 * R + n * nu * gamma_v8 - nu) / (nu ** 2 - 1) / R],
                           [-E * h * np.exp(complex(0, -1) * k1 * L) * (complex(0, -6) * R ** 2 * k1 * gamma_v1 + complex(0, 1) * h ** 2 * k1 * n + 6 * R * n * gamma_u1) / (1 + nu) / R ** 2 / 12,-E * h * np.exp(complex(0, -1) * k2 * L) * (complex(0, -6) * R ** 2 * k2 * gamma_v2 + complex(0, 1) * h ** 2 * k2 * n + 6 * R * n * gamma_u2) / (1 + nu) / R ** 2 / 12,-E * h * np.exp(complex(0, -1) * k3 * L) * (complex(0, -6) * R ** 2 * k3 * gamma_v3 + complex(0, 1) * h ** 2 * k3 * n + 6 * R * n * gamma_u3) / (1 + nu) / R ** 2 / 12,-E * h * np.exp(complex(0, -1) * k4 * L) * (complex(0, -6) * R ** 2 * k4 * gamma_v4 + complex(0, 1) * h ** 2 * k4 * n + 6 * R * n * gamma_u4) / (1 + nu) / R ** 2 / 12,-E * h * np.exp(complex(0, -1) * k5 * L) * (complex(0, -6) * R ** 2 * k5 * gamma_v5 + complex(0, 1) * h ** 2 * k5 * n + 6 * R * n * gamma_u5) / (1 + nu) / R ** 2 / 12,-E * h * np.exp(complex(0, -1) * k6 * L) * (complex(0, -6) * R ** 2 * k6 * gamma_v6 + complex(0, 1) * h ** 2 * k6 * n + 6 * R * n * gamma_u6) / (1 + nu) / R ** 2 / 12,-E * h * np.exp(complex(0, -1) * k7 * L) * (complex(0, -6) * R ** 2 * k7 * gamma_v7 + complex(0, 1) * h ** 2 * k7 * n + 6 * R * n * gamma_u7) / (1 + nu) / R ** 2 / 12,-E * h * np.exp(complex(0, -1) * k8 * L) * (complex(0, -6) * R ** 2 * k8 * gamma_v8 + complex(0, 1) * h ** 2 * k8 * n + 6 * R * n * gamma_u8) / (1 + nu) / R ** 2 / 12],
                           [complex(0, 1) * k1 * h ** 3 * E * np.exp(complex(0, -1) * k1 * L) * (-R ** 2 * k1 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * k2 * h ** 3 * E * np.exp(complex(0, -1) * k2 * L) * (-R ** 2 * k2 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * k3 * h ** 3 * E * np.exp(complex(0, -1) * k3 * L) * (-R ** 2 * k3 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * k4 * h ** 3 * E * np.exp(complex(0, -1) * k4 * L) * (-R ** 2 * k4 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * k5 * h ** 3 * E * np.exp(complex(0, -1) * k5 * L) * (-R ** 2 * k5 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * k6 * h ** 3 * E * np.exp(complex(0, -1) * k6 * L) * (-R ** 2 * k6 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * k7 * h ** 3 * E * np.exp(complex(0, -1) * k7 * L) * (-R ** 2 * k7 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2,complex(0, 1) * k8 * h ** 3 * E * np.exp(complex(0, -1) * k8 * L) * (-R ** 2 * k8 ** 2 + n ** 2 * (nu - 1) * R - n ** 2) / (12 * nu ** 2 - 12) / R ** 2],
                           [E * h ** 3 * np.exp(complex(0, -1) * k1 * L) * (R ** 2 * k1 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k2 * L) * (R ** 2 * k2 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k3 * L) * (R ** 2 * k3 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k4 * L) * (R ** 2 * k4 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k5 * L) * (R ** 2 * k5 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k6 * L) * (R ** 2 * k6 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k7 * L) * (R ** 2 * k7 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2,E * h ** 3 * np.exp(complex(0, -1) * k8 * L) * (R ** 2 * k8 ** 2 + n ** 2 * nu) / (12 * nu ** 2 - 12) / R ** 2]])
                
                K_local = np.matmul(B_mat, np.linalg.inv(A_mat))
                
                # H_n = np.array([[np.sin(n * theta),0,0,0,0,0,0,0],
                #            [0,np.cos(n * theta),0,0,0,0,0,0],
                #            [0,0,np.sin(n * theta),0,0,0,0,0],
                #            [0,0,0,np.sin(n * theta),0,0,0,0],
                #            [0,0,0,0,np.sin(n * theta),0,0,0],
                #            [0,0,0,0,0,np.cos(n * theta),0,0],
                #            [0,0,0,0,0,0,np.sin(n * theta),0],
                #            [0,0,0,0,0,0,0,np.sin(n * theta)]])
                
        return K_local        
    

    def ElementWaveNumbers(self, omega, s = True):
        '''Function that returns the wavenumbers of the system
        Input:
            s = Boolean, True means symmetric response, False antisymmetric
            omega = float, angular frequency
        Output:
            roots = 1D array of size (n_roots)
        '''
        if s == True:         # computes symmetric response
            if self.n == 0:
                # 6 roots for n=0 mode
                char_poly = [self.h ** 4 * self.R ** 4 * self.E ** 2 / (self.nu ** 2 - 1) ** 2 / 12
                             ,0
                             ,self.h ** 4 * self.R ** 4 * self.E * omega ** 2 * self.rho / (12 * self.nu ** 2 - 12)
                             ,0
                             ,-self.h ** 2 * self.R ** 2 * self.E * (-self.R ** 2 * omega ** 2 * self.rho + self.E) / (self.nu ** 2 - 1)
                             ,0
                             ,self.R ** 2 * self.h ** 2 * self.rho * omega ** 2 * (self.rho * omega ** 2 * (self.nu - 1) * (self.nu + 1) * self.R ** 2 + self.E) / (self.nu ** 2 - 1)]
                roots_list = np.roots(char_poly) 
               
            else:
                # 8 roots for n>0 modes
                char_poly = [-self.E ** 3 * self.R ** 6 * self.h ** 5 / (self.self.nu + 1) ** 3 / (self.self.nu - 1) ** 2 / 24
                            ,0
                            ,self.h ** 5 * (-omega ** 2 * self.rho * (self.self.nu + 1) * (self.self.nu - 3) * self.R ** 2 + self.E * self.n ** 2 * (self.self.nu - 1) * self.R - self.E * self.n ** 2 * (self.self.nu + 3)) * self.E ** 2 * self.R ** 4 / (self.self.nu + 1) ** 3 / (self.self.nu - 1) ** 2 / 24
                            ,0
                            ,(-6 * omega ** 2 * (self.self.nu - 1) * (self.E - self.h ** 2 * self.rho * omega ** 2 * (self.self.nu + 1) / 6) * (self.self.nu + 1) * self.rho * self.R ** 4 + self.E * self.h ** 2 * self.n ** 2 * self.rho * omega ** 2 * (self.self.nu - 1) * (self.self.nu - 3) * (self.self.nu + 1) * self.R ** 3 / 2 + 6 * ((self.self.nu - 1) * self.E - self.h ** 2 * self.n ** 2 * self.rho * omega ** 2 * (self.self.nu + 2) * (self.self.nu - 3) / 12) * (self.self.nu + 1) * self.E * self.R ** 2 + self.E ** 2 * self.h ** 2 * self.n ** 4 * (self.self.nu - 1) * self.R - self.E ** 2 * self.h ** 2 * self.n ** 4 * (self.self.nu + 2)) * self.h ** 3 * self.E * self.R ** 2 / (self.self.nu + 1) ** 3 / (self.self.nu - 1) ** 2 / 12
                            ,0
                            ,self.h ** 3 * (-12 * self.rho ** 2 * omega ** 4 * (self.self.nu - 1) * (self.self.nu - 3) * (self.self.nu + 1) ** 2 * self.R ** 6 - 2 * self.h ** 2 * self.n ** 2 * self.rho ** 2 * omega ** 4 * (self.self.nu - 1) ** 2 * (self.self.nu + 1) ** 2 * self.R ** 5 - 24 * omega ** 2 * (self.self.nu - 1) * (-self.h ** 2 * self.n ** 2 * self.rho * (self.self.nu + 1) ** 2 * omega ** 2 / 12 + self.E * (self.n ** 2 + self.self.nu + 0.3e1 / 0.2e1)) * (self.self.nu + 1) * self.rho * self.R ** 4 + self.E * self.h ** 2 * self.n ** 4 * self.rho * omega ** 2 * (self.self.nu - 1) * (self.self.nu - 3) * (self.self.nu + 1) * self.R ** 3 - self.E * self.h ** 2 * self.n ** 4 * self.rho * omega ** 2 * (self.self.nu - 3) * (self.self.nu + 2) * (self.self.nu + 1) * self.R ** 2 + self.E ** 2 * self.h ** 2 * self.n ** 6 * (self.self.nu - 1) * self.R - self.E ** 2 * self.h ** 2 * self.n ** 6 * (self.self.nu + 3)) * self.E / (self.self.nu + 1) ** 3 / (self.self.nu - 1) ** 2 / 24
                            ,0
                            ,-self.h ** 3 * (-2 * omega ** 2 * self.rho * (self.self.nu + 1) * self.R ** 2 + self.E * self.n ** 2) * (12 * self.rho ** 2 * omega ** 4 * (self.self.nu - 1) ** 2 * (self.self.nu + 1) ** 2 * self.R ** 6 + 12 * self.E * self.rho * omega ** 2 * (self.self.nu - 1) * (self.self.nu + 1) * (self.n ** 2 + 1) * self.R ** 4 + self.E * self.h ** 2 * self.n ** 4 * self.rho * omega ** 2 * (self.self.nu - 1) * (self.self.nu + 1) * self.R ** 2 + self.E ** 2 * self.h ** 2 * self.n ** 6) / self.R ** 2 / (self.self.nu + 1) ** 3 / (self.self.nu - 1) ** 2 / 24]
                roots_list = np.roots(char_poly) 
                
        if s == True:         # computes antisymmetric response
            if self.n == 0:
                # 2 roots for n=0 mode
                root = 1 / self.E * np.sqrt(2) * np.sqrt(self.E * self.rho * (self.nu + 1)) * omega
                roots_list = np.array([root,-root])  
               
            else:
                # 8 roots for n>0 modes
                char_poly = [-self.E ** 3 * self.R ** 6 * self.h ** 5 / (self.nu + 1) ** 3 / (self.nu - 1) ** 2 / 24
                             ,0
                             ,self.h ** 5 * (-omega ** 2 * self.rho * (self.nu + 1) * (self.nu - 3) * self.R ** 2 + self.E * self.n ** 2 * (self.R - 1) * self.R - self.E * self.n ** 2 * (self.nu + 3)) * self.E ** 2 * self.R ** 4 / (self.nu + 1) ** 3 / (self.nu - 1) ** 2 / 24
                             ,0
                             ,(-6 * omega ** 2 * (self.nu - 1) * (self.E - self.h ** 2 * self.rho * omega ** 2 * (self.nu + 1) / 6) * (self.nu + 1) * self.rho * self.R ** 4 + self.E * self.h ** 2 * self.n ** 2 * self.rho * omega ** 2 * (self.nu - 1) * (self.nu - 3) * (self.nu + 1) * self.R ** 3 / 2 + 6 * ((self.nu - 1) * self.E - self.h ** 2 * self.n ** 2 * self.rho * omega ** 2 * (self.nu + 2) * (self.nu - 3) / 12) * (self.nu + 1) * self.E * self.R ** 2 + self.E ** 2 * self.h ** 2 * self.n ** 4 * (self.nu - 1) * self.R - self.E ** 2 * self.h ** 2 * self.n ** 4 * (self.nu + 2)) * self.h ** 3 * self.E * self.R ** 2 / (self.nu + 1) ** 3 / (self.nu - 1) ** 2 / 12
                             ,0
                             ,self.h ** 3 * (-12 * self.rho ** 2 * omega ** 4 * (self.nu - 1) * (self.nu - 3) * (self.nu + 1) ** 2 * self.R ** 6 - 2 * self.h ** 2 * self.n ** 2 * self.rho ** 2 * omega ** 4 * (self.nu - 1) ** 2 * (self.nu + 1) ** 2 * self.R ** 5 - 24 * omega ** 2 * (self.nu - 1) * (-self.h ** 2 * self.n ** 2 * self.rho * (self.nu + 1) ** 2 * omega ** 2 / 12 + self.E * (self.n ** 2 + self.nu + 0.3e1 / 0.2e1)) * (self.nu + 1) * self.rho * self.R ** 4 + self.E * self.h ** 2 * self.n ** 4 * self.rho * omega ** 2 * (self.nu - 1) * (self.nu - 3) * (self.nu + 1) * self.R ** 3 - self.E * self.h ** 2 * self.n ** 4 * self.rho * omega ** 2 * (self.nu - 3) * (self.nu + 2) * (self.nu + 1) * self.R ** 2 + self.E ** 2 * self.h ** 2 * self.n ** 6 * (self.nu - 1) * self.R - self.E ** 2 * self.h ** 2 * self.n ** 6 * (self.nu + 3)) * self.E / (self.nu + 1) ** 3 / (self.nu - 1) ** 2 / 24
                             ,0
                             ,-self.h ** 3 * (-2 * omega ** 2 * self.rho * (self.nu + 1) * self.R ** 2 + self.E * self.n ** 2) * (12 * self.rho ** 2 * omega ** 4 * (self.nu - 1) ** 2 * (self.nu + 1) ** 2 * self.R ** 6 + 12 * self.E * self.rho * omega ** 2 * (self.nu - 1) * (self.nu + 1) * (self.n ** 2 + 1) * self.R ** 4 + self.E * self.h ** 2 * self.n ** 4 * self.rho * omega ** 2 * (self.nu - 1) * (self.nu + 1) * self.R ** 2 + self.E ** 2 * self.h ** 2 * self.n ** 6) / self.R ** 2 / (self.nu + 1) ** 3 / (self.nu - 1) ** 2 / 24]
                roots_list = np.roots(char_poly) 
                
        return roots_list



    def mode_shape(self, omega, s = True):
        '''Function that computes the displacement ratios/modeshape 
        Input:
            omega = list, angular frequency
            s = Boolean, True means symmetric response, False antisymmetric
        Output:
            gamma_u = 1D array of the size (n_roots)
            gamma_v = 1D array of the size (n_roots)
            gamma_w = 1D array of the size (n_roots)
        '''
        
        k = self.ElementWaveNumbers(omega, s)
        
        if s == True:         # computes symmetric response
            if self.n == 0:
                gamma_u = complex(0, -1) / self.R * self.E * k * self.nu / (self.nu ** 2 * omega ** 2 * self.rho + self.E * k ** 2 - omega ** 2 * self.rho)
                gamma_v = np.zeros_like(gamma_u)
                gamma_w = np.ones_like(gamma_u)
                
            else:
                gamma_u = complex(0, -1) * k * self.R * (self.nu * (-2 * self.rho * omega ** 2 * (self.nu + 1) + self.E * k ** 2) * self.R ** 2 - self.E * self.n ** 2) * self.E / ((self.rho * (self.nu ** 2 - 1) * omega ** 2 + self.E * k ** 2) * self.R ** 2 + self.E * self.n ** 2) / ((-2 * self.rho * omega ** 2 * (self.nu + 1) + self.E * k ** 2) * self.R ** 2 + self.E * self.n ** 2)
                gamma_v = -self.n * ((k ** 2 * (self.nu + 2) * self.E - 2 * self.rho * omega ** 2 * (self.nu + 1)) * self.R ** 2 + self.E * self.n ** 2) * self.E / ((self.nu ** 2 * omega ** 2 * self.rho + self.E * k ** 2 - omega ** 2 * self.rho) * self.R ** 2 + self.E * self.n ** 2) / ((-2 * self.rho * omega ** 2 * (self.nu + 1) + self.E * k ** 2) * self.R ** 2 + self.E * self.n ** 2)
                gamma_w = np.ones_like(gamma_u)
                
        else:                 # computes symmetric response
            if self.n == 0:
                gamma_v = np.ones_like(k)
                gamma_u = np.zeros_like(k)
                gamma_w = np.zeros_like(k)
                
            else:
                gamma_u = complex(0, -1) * k * self.R * (self.nu * (-2 * self.rho * omega ** 2 * (self.nu + 1) + self.E * k ** 2) * self.R ** 2 - self.E * self.n ** 2) * self.E / ((self.rho * (self.nu ** 2 - 1) * omega ** 2 + self.E * k ** 2) * self.R ** 2 + self.E * self.n ** 2) / ((-2 * self.rho * omega ** 2 * (self.nu + 1) + self.E * k ** 2) * self.R ** 2 + self.E * self.n ** 2)
                gamma_v = self.n * ((k ** 2 * (self.nu + 2) * self.E - 2 * self.rho * omega ** 2 * (self.nu + 1)) * self.R ** 2 + self.E * self.n ** 2) * self.E / ((self.nu ** 2 * omega ** 2 * self.rho + self.E * k ** 2 - omega ** 2 * self.rho) * self.R ** 2 + self.E * self.n ** 2) / ((-2 * self.rho * omega ** 2 * (self.nu + 1) + self.E * k ** 2) * self.R ** 2 + self.E * self.n ** 2)
                gamma_w = np.ones_like(gamma_u)
                
        return gamma_u, gamma_v, gamma_w



