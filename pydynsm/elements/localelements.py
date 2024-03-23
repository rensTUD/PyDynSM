# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:11:13 2024

@author: rensv
"""

# %% Import dependencies

import numpy as np
# from .baseelement import BaseElement

# %% Class definitions
'''
Considerations:

total stiffness will be 6x6 (i.e. in 2D):
    [N_left, V_left, M_left, N_right, V_right, M_right] = K.[u_left, w_left, phi_left, u_right, w_right, phi_right]

where the contribution of the rod is as:
    K = [K_r_ll, 0, 0, K_r_lr, 0, 0; 
         0     , 0, 0, 0   , 0, 0;
         0     , 0, 0, 0   , 0, 0;
         K_r_rl, 0, 0, K_r_rr, 0, 0;
         0     , 0, 0, 0   , 0, 0;
         0     , 0, 0, 0   , 0, 0]
'''

# %% standard structural element class
class StructuralElement():
    '''
    Parent class of any other element, all functions are inherited by its child classes so add functions if necessary.
    '''
    def __init__(self,dofs):
        
        # leave this here
        self.Ndof = 3
        
        self.dofs = dofs
        self.nodal_dofs = self.GetNodalDofs(self.dofs,self.Ndof)
    
    def GetNodalDofs(self, dofs, Ndof):
        '''
        Translates the dofs of an element to the nodal dofs;
        
        Example:
            Rod element has dof = 0 (axial only)
            This will give: nodal_dofs = [0,3]
        '''
        
        nodal_dofs = []
        # Process each DOF to map it to the current and next node's corresponding DOFs
        for dof in dofs:
            # Add the DOF for the current node
            nodal_dofs.append(dof)
        # After processing each DOF for the current node, add their counterparts in the next node
        for dof in dofs:
            nodal_dofs.append(dof + Ndof)
                
        return nodal_dofs
    
    def FullStiffness(self, omega):
        '''
        Function that assembles the full stiffness matrix based on the local stiffness matrix. 
        
        For example, it will translate the 2x2 K_local of the rod to the full 6x6 matrix which it is in 2D
        
        '''
        
        # initialise 6x6 empty complex matrix
        K_full = np.zeros((6,6), complex)
        
        # calculate local stiffness matrix
        K_local = self.LocalStiffness(omega)
        
        # assign to correct locations in full matrix
        for i, row in enumerate(self.nodal_dofs):
            for j, col in enumerate(self.nodal_dofs):
                K_full[row, col] = K_local[i, j]
                
        return K_full 
    
    def FullDistributedLoad(self, q, omega):
        
        # initialise 6x1 empty complex vector 
        q_full = np.zeros(6,complex)
        
        # calculate local load vector
        q_local = self.LocalDistributedLoad(q, omega)
        
        # assign to correct locations in full vector
        for i, row in enumerate(self.nodal_dofs):
            q_full[row] = q_local[i]
        
        return q_full
    
    def FullDisplacement(self, u_node, omega):
        
        # initialise 6x1 empty complex vector 
        u_full = np.zeros(6,complex)
        
        # calculate local load vector
        u_local = self.LocalDisplacements(u_node, omega)
        
        # assign to correct locations in full vector
        for i, row in enumerate(self.nodal_dofs):
            u_full[row] = u_local[i]
        
        return u_full
    
    def FullElementDisplacements(self, u_nodes_global, omega, num_points):
        
        # intilise full u_elem
        u_elem_full = [None] * self.Ndof
        
        # calculate local u_elem
        u_elem_local = self.LocalElementDisplacements(u_nodes_global, omega, num_points)
        
        # assign
        for i, dof in enumerate(self.dofs):
            u_elem_full[dof] = u_elem_local[i]
        
        return u_elem_full

# %% specific local element classes

class EB_Beam(StructuralElement):
    
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
        beta_b = self.beta_b(omega)
        
        
        K_local = np.empty((4,4), complex)
        # we can also create the stiffness matrix like this, probably easier to get direct output from Maple but maybe less clear and less easy to spot mistakes (if present)
        K_local = np.array([ [-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1), EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1), EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1), EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)], 
                           [EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)]])
        
        return K_local        
  
    def beta_b(self,omega):
        '''
        Determines the wavenumbers
        '''
        
        beta_b = (omega**2 * self.rhoA / (self.EI*(1+2j*self.ksi)))
        
        return beta_b
    
    def LocalDistributedLoad(self, q, omega):
        '''
        add a distributed load to the local element
        
        '''
        
        # assign load to itself to keep track
        self.q = q
        
        # assign local variables for ease of coding
        L = self.L

        # determine wavenumber
        beta_b = self.beta_b(omega)
        
        
        el = np.array([ [q*((np.cos(beta_b*L) - 1.0)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0))],
                        [-q*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0))],
                        [q*((np.cos(beta_b*L) - 1)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0))], 
                        [q*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0))] 
                        ])
        
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
        phi= self.rotation(x,omega,C)
        
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
        beta_b = self.beta_b(omega)
        
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
        beta_b = self.beta_b(omega)
        
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
        beta_b = self.beta_b(omega)
        
        # check if C is input        
        if C == None:
            C = self.Coefficients(u_node_local, omega)
        
        # displacements
        phi = -C[3] * beta_b * np.sin(beta_b * x) + C[2] * beta_b * np.cos(beta_b * x) + C[0] * beta_b * np.sinh(beta_b * x) + C[1] * beta_b * np.cosh(beta_b * x)
        
        return -phi   
    
    
class Rod_1D(StructuralElement):
            
    def __init__(self, rhoA, EA, L, ksi = None):
        
        # define what dofs the rod contributes to and initialise
        dofs = [0]
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
        beta_r = self.beta_r(omega)
    
        # initialise empty stiffness vector
        K_local = np.empty((2,2),complex)
        
        # assign variabes
        K_local[0,0] = EA * beta_r * np.cos(beta_r * L) / np.sin(beta_r * L)  
        K_local[1,1] = EA * beta_r * np.cos(beta_r * L) / np.sin(beta_r * L)
        K_local[1,0] = - EA * beta_r / np.sin(beta_r * L)
        K_local[0,1] = - EA * beta_r / np.sin(beta_r * L) 
    
        return K_local
    
    def beta_r(self, omega):
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
        beta_r = self.beta_r(omega)
        
        
        el = np.array([ [(np.cos(beta_r*L) - 1.0)*q/(np.sin(beta_r*L)*beta_r)], 
                        [(np.cos(beta_r*L) - 1.0)*q/(np.sin(beta_r*L)*beta_r)] 
                        ])
        
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