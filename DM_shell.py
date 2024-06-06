#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:52:43 2024

@author: lukas
"""

# a function that returns the root based on the input parameters

import numpy as np


def roots(R,nu,h,E,n, rho, omega):
    '''Function that returns the roots/wavenumbers of the system
    Input:
        R = float, radius of the shell midsurface
        nu = float, poissons ratio
        h = float, thickness of the shell
        E = float or complex, Young's modulus of the material
        n = int, modeshape
        rho = float, density of the material
        omega = list, angular frequency
    Output:
        roots = 2D array of size (omega,n_roots)
    '''
    
    roots_list= []
    for ii in range(len(omega)):
        char_poly = [-E ** 3 * R ** 6 * h ** 5 / (1 + nu) / (nu ** 2 - 1) ** 2 / 24
                     ,0
                     ,-(E ** 2 * R ** 8 * h ** 2 * nu ** 2 * omega[ii] ** 2 * rho - 2 * E ** 2 * R ** 8 * h ** 2 * nu * omega[ii] ** 2 * rho - E ** 3 * R ** 7 * h ** 2 * n ** 2 * nu - 3 * E ** 2 * R ** 8 * h ** 2 * omega[ii] ** 2 * rho + E ** 3 * R ** 7 * h ** 2 * n ** 2 + E ** 3 * R ** 6 * h ** 2 * n ** 2 * nu + 3 * E ** 3 * R ** 6 * h ** 2 * n ** 2) * h ** 3 / R ** 2 / (1 + nu) / (nu ** 2 - 1) ** 2 / 24
                     ,0
                     ,-(-2 * E * R ** 8 * h ** 2 * nu ** 3 * omega[ii] ** 4 * rho ** 2 - E ** 2 * R ** 7 * h ** 2 * n ** 2 * nu ** 3 * omega[ii] ** 2 * rho - 2 * E * R ** 8 * h ** 2 * nu ** 2 * omega[ii] ** 4 * rho ** 2 + 3 * E ** 2 * R ** 7 * h ** 2 * n ** 2 * nu ** 2 * omega[ii] ** 2 * rho + E ** 2 * R ** 6 * h ** 2 * n ** 2 * nu ** 3 * omega[ii] ** 2 * rho + 2 * E * R ** 8 * h ** 2 * nu * omega[ii] ** 4 * rho ** 2 + E ** 2 * R ** 7 * h ** 2 * n ** 2 * nu * omega[ii] ** 2 * rho + 2 * E * R ** 8 * h ** 2 * omega[ii] ** 4 * rho ** 2 - 3 * E ** 2 * R ** 7 * h ** 2 * n ** 2 * omega[ii] ** 2 * rho - 7 * E ** 2 * R ** 6 * h ** 2 * n ** 2 * nu * omega[ii] ** 2 * rho - 2 * E ** 3 * R ** 5 * h ** 2 * n ** 4 * nu + 12 * E ** 2 * R ** 8 * nu ** 2 * omega[ii] ** 2 * rho - 6 * E ** 2 * R ** 6 * h ** 2 * n ** 2 * omega[ii] ** 2 * rho + 2 * E ** 3 * R ** 5 * h ** 2 * n ** 4 + 2 * E ** 3 * R ** 4 * h ** 2 * n ** 4 * nu + 4 * E ** 3 * R ** 4 * h ** 2 * n ** 4 - 12 * E ** 2 * R ** 8 * omega[ii] ** 2 * rho - 12 * E ** 3 * R ** 6 * nu ** 2 + 12 * E ** 3 * R ** 6) * h ** 3 / R ** 2 / (1 + nu) / (nu ** 2 - 1) ** 2 / 24
                     ,0
                     ,-(2 * E * R ** 7 * h ** 2 * n ** 2 * nu ** 4 * omega[ii] ** 4 * rho ** 2 - 2 * E * R ** 6 * h ** 2 * n ** 2 * nu ** 4 * omega[ii] ** 4 * rho ** 2 - 4 * E * R ** 7 * h ** 2 * n ** 2 * nu ** 2 * omega[ii] ** 4 * rho ** 2 - 4 * E * R ** 6 * h ** 2 * n ** 2 * nu ** 3 * omega[ii] ** 4 * rho ** 2 - E ** 2 * R ** 5 * h ** 2 * n ** 4 * nu ** 3 * omega[ii] ** 2 * rho + 12 * E * R ** 8 * nu ** 4 * omega[ii] ** 4 * rho ** 2 + 3 * E ** 2 * R ** 5 * h ** 2 * n ** 4 * nu ** 2 * omega[ii] ** 2 * rho + E ** 2 * R ** 4 * h ** 2 * n ** 4 * nu ** 3 * omega[ii] ** 2 * rho - 24 * E * R ** 8 * nu ** 3 * omega[ii] ** 4 * rho ** 2 + 2 * E * R ** 7 * h ** 2 * n ** 2 * omega[ii] ** 4 * rho ** 2 + 4 * E * R ** 6 * h ** 2 * n ** 2 * nu * omega[ii] ** 4 * rho ** 2 + E ** 2 * R ** 5 * h ** 2 * n ** 4 * nu * omega[ii] ** 2 * rho - 48 * E * R ** 8 * nu ** 2 * omega[ii] ** 4 * rho ** 2 + 2 * E * R ** 6 * h ** 2 * n ** 2 * omega[ii] ** 4 * rho ** 2 - 3 * E ** 2 * R ** 5 * h ** 2 * n ** 4 * omega[ii] ** 2 * rho - 7 * E ** 2 * R ** 4 * h ** 2 * n ** 4 * nu * omega[ii] ** 2 * rho + 24 * E * R ** 8 * nu * omega[ii] ** 4 * rho ** 2 - E ** 3 * R ** 3 * h ** 2 * n ** 6 * nu + 24 * E ** 2 * R ** 6 * n ** 2 * nu ** 2 * omega[ii] ** 2 * rho - 6 * E ** 2 * R ** 4 * h ** 2 * n ** 4 * omega[ii] ** 2 * rho + 36 * E * R ** 8 * omega[ii] ** 4 * rho ** 2 + E ** 3 * R ** 3 * h ** 2 * n ** 6 + E ** 3 * R ** 2 * h ** 2 * n ** 6 * nu + 24 * E ** 2 * R ** 6 * nu ** 3 * omega[ii] ** 2 * rho + 3 * E ** 3 * R ** 2 * h ** 2 * n ** 6 - 24 * E ** 2 * R ** 6 * n ** 2 * omega[ii] ** 2 * rho + 36 * E ** 2 * R ** 6 * nu ** 2 * omega[ii] ** 2 * rho - 24 * E ** 2 * R ** 6 * nu * omega[ii] ** 2 * rho - 36 * E ** 2 * R ** 6 * omega[ii] ** 2 * rho) * h ** 3 / R ** 2 / (1 + nu) / (nu ** 2 - 1) ** 2 / 24
                     ,0
                     ,-(-24 * R ** 8 * omega[ii] ** 6 * rho ** 3 * nu ** 5 - 24 * R ** 8 * nu ** 4 * omega[ii] ** 6 * rho ** 3 - 2 * E * R ** 4 * h ** 2 * n ** 4 * nu ** 3 * omega[ii] ** 4 * rho ** 2 + 48 * R ** 8 * nu ** 3 * omega[ii] ** 6 * rho ** 3 + 12 * E * R ** 6 * n ** 2 * nu ** 4 * omega[ii] ** 4 * rho ** 2 - 2 * E * R ** 4 * h ** 2 * n ** 4 * nu ** 2 * omega[ii] ** 4 * rho ** 2 + 48 * R ** 8 * nu ** 2 * omega[ii] ** 6 * rho ** 3 - 24 * E * R ** 6 * n ** 2 * nu ** 3 * omega[ii] ** 4 * rho ** 2 + 2 * E * R ** 4 * h ** 2 * n ** 4 * nu * omega[ii] ** 4 * rho ** 2 - 24 * R ** 8 * nu * omega[ii] ** 6 * rho ** 3 + E ** 2 * R ** 2 * h ** 2 * n ** 6 * nu ** 2 * omega[ii] ** 2 * rho - 48 * E * R ** 6 * n ** 2 * nu ** 2 * omega[ii] ** 4 * rho ** 2 + 2 * E * R ** 4 * h ** 2 * n ** 4 * omega[ii] ** 4 * rho ** 2 - 24 * R ** 8 * omega[ii] ** 6 * rho ** 3 - 2 * E ** 2 * R ** 2 * h ** 2 * n ** 6 * nu * omega[ii] ** 2 * rho + 24 * E * R ** 6 * n ** 2 * nu * omega[ii] ** 4 * rho ** 2 - 24 * E * R ** 6 * nu ** 3 * omega[ii] ** 4 * rho ** 2 + 12 * E ** 2 * R ** 4 * n ** 4 * nu ** 2 * omega[ii] ** 2 * rho - 3 * E ** 2 * R ** 2 * h ** 2 * n ** 6 * omega[ii] ** 2 * rho + 36 * E * R ** 6 * n ** 2 * omega[ii] ** 4 * rho ** 2 - 24 * E * R ** 6 * nu ** 2 * omega[ii] ** 4 * rho ** 2 + 24 * E * R ** 6 * nu * omega[ii] ** 4 * rho ** 2 + E ** 3 * h ** 2 * n ** 8 - 12 * E ** 2 * R ** 4 * n ** 4 * omega[ii] ** 2 * rho + 12 * E ** 2 * R ** 4 * n ** 2 * nu ** 2 * omega[ii] ** 2 * rho + 24 * E * R ** 6 * omega[ii] ** 4 * rho ** 2 - 12 * E ** 2 * R ** 4 * n ** 2 * omega[ii] ** 2 * rho) * h ** 3 / R ** 2 / (1 + nu) / (nu ** 2 - 1) ** 2 / 24]

        roots_temp = np.roots(char_poly)
        
        roots_list.append(roots_temp)
        
    return np.array(roots_list)



def mode_shape(R,nu,h,E,n, rho, omega):
    '''Function that computes the displacement ratios/modeshape 
    Input:
        R = float, radius of the shell midsurface
        nu = float, poissons ratio
        h = float, thickness of the shell
        E = float or complex, Young's modulus of the material
        n = int, modeshape
        rho = float, density of the material
        omega = list, angular frequency
    Output:
        gamma_u = 2D array of the size (omega,n_roots)
        gamma_v = 2D array of the size (omega,n_roots)
        gamma_w = 2D array of the size (omega,n_roots)
    '''
    
    R = R
    nu = nu
    h = h
    E = E
    n = n
    rho = rho 
    omega = omega

        
    gamma_u= np.zeros((len(omega),8),dtype=complex)
    gamma_v= np.zeros((len(omega),8),dtype=complex)
    gamma_w = np.ones((len(omega),8),dtype=complex)
    
    
    wavenumbers = roots(R,nu,h,E,n, rho, omega)
    omega=np.array(omega)
    
    for ii in range(8):
        k = wavenumbers[:,ii]
    
        gamma_u[:,ii] = complex(0, -2) * (1 + nu) * (-2 * R ** 2 * nu ** 2 * omega ** 2 * rho + E * R ** 2 * k ** 2 * nu - E * R ** 2 * k ** 2 + 2 * R ** 2 * omega ** 2 * rho - 2 * E * n ** 2) / (-2 * R ** 4 * omega ** 4 * rho ** 2 * nu ** 3 + E * R ** 4 * k ** 2 * nu ** 2 * omega ** 2 * rho - 2 * R ** 4 * nu ** 2 * omega ** 4 * rho ** 2 - 2 * E * R ** 4 * k ** 2 * nu * omega ** 2 * rho + 2 * R ** 4 * nu * omega ** 4 * rho ** 2 + E ** 2 * R ** 4 * k ** 4 - 3 * E * R ** 4 * k ** 2 * omega ** 2 * rho + E * R ** 2 * n ** 2 * nu ** 2 * omega ** 2 * rho + 2 * R ** 4 * omega ** 4 * rho ** 2 - 2 * E * R ** 2 * n ** 2 * nu * omega ** 2 * rho + 2 * k ** 2 * E ** 2 * n ** 2 * R ** 2 - 3 * E * R ** 2 * n ** 2 * omega ** 2 * rho + E ** 2 * n ** 4) * k * E * nu * R / (2 * nu ** 2 - 2) + complex(0, -2) * (1 + nu) ** 2 * R * n ** 2 * E ** 2 * k / (-2 * R ** 4 * omega ** 4 * rho ** 2 * nu ** 3 + E * R ** 4 * k ** 2 * nu ** 2 * omega ** 2 * rho - 2 * R ** 4 * nu ** 2 * omega ** 4 * rho ** 2 - 2 * E * R ** 4 * k ** 2 * nu * omega ** 2 * rho + 2 * R ** 4 * nu * omega ** 4 * rho ** 2 + E ** 2 * R ** 4 * k ** 4 - 3 * E * R ** 4 * k ** 2 * omega ** 2 * rho + E * R ** 2 * n ** 2 * nu ** 2 * omega ** 2 * rho + 2 * R ** 4 * omega ** 4 * rho ** 2 - 2 * E * R ** 2 * n ** 2 * nu * omega ** 2 * rho + 2 * k ** 2 * E ** 2 * n ** 2 * R ** 2 - 3 * E * R ** 2 * n ** 2 * omega ** 2 * rho + E ** 2 * n ** 4) / (2 * nu ** 2 - 2)
        
        gamma_v[:,ii] =-2 * (1 + nu) ** 2 * R ** 2 * n * E ** 2 * k ** 2 / (-2 * R ** 4 * omega ** 4 * rho ** 2 * nu ** 3 + E * R ** 4 * k ** 2 * nu ** 2 * omega ** 2 * rho - 2 * R ** 4 * nu ** 2 * omega ** 4 * rho ** 2 - 2 * E * R ** 4 * k ** 2 * nu * omega ** 2 * rho + 2 * R ** 4 * nu * omega ** 4 * rho ** 2 + E ** 2 * R ** 4 * k ** 4 - 3 * E * R ** 4 * k ** 2 * omega ** 2 * rho + E * R ** 2 * n ** 2 * nu ** 2 * omega ** 2 * rho + 2 * R ** 4 * omega ** 4 * rho ** 2 - 2 * E * R ** 2 * n ** 2 * nu * omega ** 2 * rho + 2 * k ** 2 * E ** 2 * n ** 2 * R ** 2 - 3 * E * R ** 2 * n ** 2 * omega ** 2 * rho + E ** 2 * n ** 4) * nu / (2 * nu ** 2 - 2) + 2 * (1 + nu) * (2 * R ** 2 * nu ** 2 * omega ** 2 * rho + 2 * E * R ** 2 * k ** 2 - 2 * R ** 2 * omega ** 2 * rho - E * n ** 2 * nu + E * n ** 2) / (-2 * R ** 4 * omega ** 4 * rho ** 2 * nu ** 3 + E * R ** 4 * k ** 2 * nu ** 2 * omega ** 2 * rho - 2 * R ** 4 * nu ** 2 * omega ** 4 * rho ** 2 - 2 * E * R ** 4 * k ** 2 * nu * omega ** 2 * rho + 2 * R ** 4 * nu * omega ** 4 * rho ** 2 + E ** 2 * R ** 4 * k ** 4 - 3 * E * R ** 4 * k ** 2 * omega ** 2 * rho + E * R ** 2 * n ** 2 * nu ** 2 * omega ** 2 * rho + 2 * R ** 4 * omega ** 4 * rho ** 2 - 2 * E * R ** 2 * n ** 2 * nu * omega ** 2 * rho + 2 * k ** 2 * E ** 2 * n ** 2 * R ** 2 - 3 * E * R ** 2 * n ** 2 * omega ** 2 * rho + E ** 2 * n ** 4) * E * n / (2 * nu ** 2 - 2)
        
    
    return gamma_u, gamma_v, gamma_w



mode_shapeu,mode_shapev,mode_shapew = mode_shape(4,0.3,0.02,20e10,1,7580,np.linspace(1,1000,10000))


def DSM_shell(R,nu,h,E,n, rho, omega):
    





