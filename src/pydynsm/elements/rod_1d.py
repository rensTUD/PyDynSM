# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:27:16 2024

@author: rensv
"""

# import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition

@ElementFactory.ElementType('Rod')
class Rod_1D(StructuralElement):

    def __init__(self, rhoA, EA, L, ksi = None):
        
        # define what dofs the rod contributes to and initialise
        dofs = ['x']
        super().__init__(dofs)
        
        # Initialise lcoal rod element with necessary parameters
        self.rhoA = rhoA
        self.EA = EA
        self.L = L
        self.ksi = ksi if ksi is not None else 0.01 # assisgn ksi if given otherwise assign a default value (think of changing this)
    
    def LocalStiffness(self, omega):
        '''
        Determines the stiffness of the rod. As its 1D, stiffness will be 2x2 matrix as:
            [N_left, N_right] = [K_ll, K_lr; K_rl, K_rr].[u_left, u_right]
            
        '''    
        
        # assign local variables for ease of coding
        EA = self.EA
        L = self.L
                
        # determine wavenumber
        beta_r = self.ElementWaveNumbers(omega)
    
        # initialise empty stiffness vector
        K_local = np.empty((2,2),complex)
        
        # assign variabes
        K_local[0,0] = EA * beta_r * np.cos(beta_r * L) / np.sin(beta_r * L)  
        K_local[1,1] = EA * beta_r * np.cos(beta_r * L) / np.sin(beta_r * L)
        K_local[1,0] = - EA * beta_r / np.sin(beta_r * L)
        K_local[0,1] = - EA * beta_r / np.sin(beta_r * L) 
    
        return K_local
    
    def ElementWaveNumbers(self, omega):
        '''
        function that calculates the wavenumbers for frequency omega
        '''
        
        beta_r = omega / (self.EA*(1+2j*self.ksi)/self.rhoA) ** 0.5
        
        return beta_r
    
    def LocalDistributedLoad(self, q, omega):
        '''
        add a distributed load to the local element
        
        '''
        
        # assign local variables for ease of coding
        L = self.L
        

        # determine wavenumber
        beta_r = self.ElementWaveNumbers(omega)
        
        
        el = [(np.cos(beta_r*L) - 1.0)*q/(np.sin(beta_r*L)*beta_r), (np.cos(beta_r*L) - 1.0)*q/(np.sin(beta_r*L)*beta_r)]
        
        return el

    def LocalElementDisplacements(self, u_nodes_global, omega, num_points):
        '''
        function that calculates the local displacement u(s)
        '''
        pass
    
    def Coefficients(self, u_node_local, omega):
        '''
        Calculates the coefficients of the general solution
        '''
        pass
    
    def displacement(self, x, omega, C = None, u_node_local = None):
        '''
        Gets the transverse displacments of the Euler-Bernoulli beam
        
        if C is not given, then calculate it based on u_node_local.
        '''
        pass
    
 