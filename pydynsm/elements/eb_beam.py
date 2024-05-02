# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:27:06 2024

@author: rensv
"""
# import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition

@ElementFactory.ElementType('EB Beam')
class EB_Beam(StructuralElement):
    
    # element_name = 'EB Beam'
    
    def __init__(self, rhoA, EI, L, ksi = None):
        
        # define what dofs the eb beam contributes to and initialise
        dofs = [1,2]
        super().__init__(dofs)
        
        # Initialise lcoal rod element with necessary parameters
        self.ksi = ksi if ksi is not None else 0.01 # assisgn ksi if given otherwise assign a default value (think of changing this)
        self.rhoA = rhoA
        self.EI = EI
        self.L = L
        
        
    def LocalStiffness(self, omega):
        '''
        Determines the stiffness of the rod. As its 1D, stiffness will be 2x2 matrix as:
            [V_left, M_left, V_right, M_right] = K.[w_left, phi_left, w_right, phi_right]
        
        where:
            K = [K_V_ll, K_M_ll, K_V_lr, K_M_lr;
                 K_V_rl, K_M_rl, K_V_rr, K_M_rr]
        '''    
        
        # assign local variables for ease of coding
        EI = self.EI*(1+2j*self.ksi)
        L = self.L
        beta_b = self.ElementWaveNumbers(omega)
        
        
        K_local = np.empty((4,4), complex)
        # we can also create the stiffness matrix like this, probably easier to get direct output from Maple but maybe less clear and less easy to spot mistakes (if present)
        K_local = np.array([ [-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1), EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1), EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1), EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)], 
                           [EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)]])
        
        return K_local        
  
    def ElementWaveNumbers(self,omega):
        '''
        Determines the wavenumbers
        '''
        
        beta_b = (omega**2 * self.rhoA / (self.EI*(1+2j*self.ksi)))
        
        return beta_b
    
    def LocalDistributedLoad(self, q, omega):
        '''
        add a distributed load to the local element
        
        q = [q_z;
             q_phi]
        
        '''
        
        # assign load to itself to keep track
        self.q = q
        
        # assign local variables for ease of coding
        L = self.L

        # determine wavenumber
        beta_b = self.ElementWaveNumbers(omega)
        
        # extract loads
        q_z = q[0]
        q_phi = q[1]
        
        # TODO - check for correctness
        el = [ q_z*((np.cos(beta_b*L) - 1.0)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)),
                        -q_phi*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)),
                        q_z*((np.cos(beta_b*L) - 1)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)), 
                        q_phi*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)) 
                        ]
        
        return el  
    
    def LocalElementDisplacements(self, u_nodes_global, omega, num_points):
        '''
        function that calculates the local displacements w(s) and phi(s)
        '''
        
        # get local axis to evaluate on
        L = self.L
        x = np.linspace ( 0.0, L, num_points )
                
        # determine nodal displacement in local axis
        u_node_local = self.R @ u_nodes_global
        
        # extract only the needed displacements
        u_node_local = u_node_local(self.dofs)
        
                
        # calculate coeficients
        C = self.Coefficients(u_node_local, omega)
        
        # get displacement
        w = self.displacement(x,omega,C)
        
        # get rotations
        phi = self.rotation(x,omega,C)
        
        return [w, phi]
    
    def Coefficients(self, u_node_local, omega):
        '''
        Calculates the coefficients of the general solution
        '''
        # read all the variables
        ksi = self.ksi 
        rhoA = self.rhoA  
        EI = self.EI * (1 + 2j * ksi)
        L = self.L
        beta_b = self.ElementWaveNumbers(omega)
        
        # TODO - should change this based on new derivations..
        
        # get distributed load value
        q = self.q[1]
        # should be like: 
        # q_b, q_m = self.
    
        # calculate the coefficients
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (u_node_local + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        return C
    
    def displacement(self, x, omega, C = None, u_node_local = None):
        '''
        Gets the transverse displacments of the Euler-Bernoulli beam
        
        if C is not given, then calculate it based on u_node_local.
        '''
        
        # TODO - should change this based on new derivations..
        
        ksi = self.ksi
        EI = self.EI*(1 + 2j * ksi)
        q = self.q[1]
        beta_b = self.ElementWaveNumbers(omega)
        
        # check if C is input        
        if C == None:
            C = self.Coefficients(u_node_local, omega)
        
        # displacements
        w = C[3] * np.cos(beta_b * x) + C[2] * np.sin(beta_b * x) + C[0] * np.cosh(beta_b * x) + C[1] * np.sinh(beta_b * x) - q / EI / beta_b**4
        
        return w

    def rotation(self, x, omega, C = None, u_node_local = None):
        '''
        Gets the rotations of the Euler-Bernoulli beam
        '''
        
        # TODO - should change this based on new derivations..
        
        ksi = self.ksi
        EI = self.EI*(1 + 2j * ksi)
        q = self.q[1]
        beta_b = self.ElementWaveNumbers(omega)
        
        # check if C is input        
        if C == None:
            C = self.Coefficients(u_node_local, omega)
        
        # displacements
        phi = -C[3] * beta_b * np.sin(beta_b * x) + C[2] * beta_b * np.cos(beta_b * x) + C[0] * beta_b * np.sinh(beta_b * x) + C[1] * beta_b * np.cosh(beta_b * x)
        
        return -phi   
    

