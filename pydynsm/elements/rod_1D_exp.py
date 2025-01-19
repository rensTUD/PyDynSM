# -*- coding: utf-8 -*-

"""
Created on Thu Apr 18 17:27:16 2024

@author: rensv
"""

# Import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition


@ElementFactory.ElementType('Rod')#, required_parameters = ['rho', 'A', 'E', 'L'])
class Rod1D(StructuralElement):
    """Class for 1D rod element."""

    element_name = 'Rod'

    def __init__(self, rho, A, E, L, ksi=None):
        """
        Initialise a new instance of Rod1D.

        Input:
            rho: value. Density of the element's material [kg/m^3]
            A:   value. Area of the element [m^2]
            E:   value. E-modulus of the element's material [Pa]
            L:   value. Length of element [m]
            ksi: value. Damping of the element's material [-]
        """
        # define what dofs the rod contributes to and initialise
        dofs = ['x']
        super().__init__(dofs)

        # Initialise local rod element with necessary parameters
        self.rho = rho
        self.A = A
<<<<<<< HEAD
        self.ksi = ksi if ksi is not None else 0.01
        self.E = E * (1+2j*self.ksi)
        self.L = L
        # assisgn ksi if given otherwise assign a default value
        
=======
        # self.E = E # please see the line below for complex E value
        self.L = L
        # assisgn ksi if given otherwise assign a default value
        self.ksi = ksi if ksi is not None else 0.01
        self.E = E*(1+2j*self.ksi)
>>>>>>> f7183ef (updated E with a complex valued E for 1D_rod_exp)
        # set q standard to 0
        self.q = np.zeros(len(dofs))

    def LocalStiffness(self, omega):
        """
        Determine the stiffness of the rod.

        As its 1D, stiffness will be 2x2 matrix as:
            [F_left, F_right] = K.[u_left, u_right]
        where:
            K = [K_V_ll, K_V_lr;
                 K_V_rl, K_V_rr]
        Input:
            omega: array. Range of frequencies of analysis
        Output:
            K_local: matrix. Dynamic stiffness matrix (also K_dyn)
        """
        # Assign local variables for ease of coding
        A = self.A
        E = self.E
        L = self.L

        # determine wavenumber
        alpha_1, alpha_2 = self.ElementWaveNumbers(omega)

        # Initialize K matrix
        K_local = np.empty((2, 2), complex)
        # Copy-paste the matrix directly from Maple (also called K_dyn)
        K_local = np.array([[1j*A*E*(np.exp(-1j*alpha_2*L)*alpha_1 - np.exp(-1j*alpha_1*L)*alpha_2)/(np.exp(-1j*alpha_2*L) - np.exp(-1j*alpha_1*L)),
                             1j*A*E*(alpha_1 - alpha_2)/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))],
                            [1j*A*E*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2))/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L)),
                             -1j*A*E*(np.exp(-1j*alpha_2*L)*alpha_2 - np.exp(-1j*alpha_1*L)*alpha_1)/(np.exp(-1j*alpha_2*L) - np.exp(-1j*alpha_1*L))]])

        return K_local

    def ElementWaveNumbers(self, omega):
        """
        Determine the wavenumbers of 1D rod.

        Input:
            omega: array. Range of frequencies of analysis.
        Output:
            alpha_1, alpha_2,: values. Wavenumbers
        """
        rho = self.rho
        E = self.E
        # Copy-paste from Maple document of Wavenumbers
        # 2nd-order derivative in x, so 2 wavenumbers
        alpha_1 = (np.sqrt(rho/E)*omega)
        alpha_2 = (-np.sqrt(rho/E)*omega)

        return alpha_1, alpha_2

    def LocalDistributedLoad(self, q, omega):
        """
        Add a distributed load to the local element.

        Input:
            q:     array. Distributed load. With a shape such as:q = [q_x]
            omega: array. Range of frequencies of analyis
        Output:
            el:    array. Force vector
        """
        # assign load to itself to keep track
        self.q = q

        # assign local variables for ease of coding
        rho = self.rho
        E = self.E
        L = self.L

        # extract loads
        q_x = q[0]

        # determine wavenumber
        alpha_1, alpha_2 = self.ElementWaveNumbers(omega)

        # TODO - check (this is correct from maple, need to be checked later via numerical examples)
        el = np.array([-1j*E*(np.exp(-1j*alpha_2*L)*alpha_1 - np.exp(-1j*alpha_1*L)*alpha_2)/(np.exp(-1j*alpha_2*L) - np.exp(-1j*alpha_1*L))*q_x/rho/omega**2 + -1j*E*(alpha_1 - alpha_2)/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))*q_x/rho/omega**2,
                       -1j*E*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2))/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))*q_x/rho/omega**2 + 1j*E*(np.exp(-1j*alpha_2*L)*alpha_2 - np.exp(-1j*alpha_1*L)*alpha_1)/(np.exp(-1j*alpha_2*L) - np.exp(-1j*alpha_1*L))*q_x/rho/omega**2])

        return el

    def LocalElementDisplacements(self, u_nodes_local, omega, num_points):
        """
        Calcualte the local displacements u(s).

        Input:
            u_nodes_global: array. The nodes in global coordinates
            omega:          array. Range of frequencies of analysis
            num_points:     value. Number of points to divide the element in.
        Output:
            u: array. Amplitude of vertical displacement
        """
        # get local axis to evaluate on
        L = self.L
        x = np.linspace(0.0, L, num_points)

        # calculate coeficients
        C = self.Coefficients(u_nodes_local, omega)

        # get displacement
        u = self.displacement(x, omega, C)

        return [u]
    def LocalElementForces(self, u_nodes_local, omega, num_points):
        """
        Calcualte the local displacements u(s).

        Input:
            u_nodes_global: array. The nodes in global coordinates
            omega:          array. Range of frequencies of analysis
            num_points:     value. Number of points to divide the element in.
        Output:
            u: array. Amplitude of vertical displacement
        """
        # get local axis to evaluate on
        L = self.L
        x = np.linspace(0.0, L, num_points)

        # calculate coeficients
        C = self.Coefficients(u_nodes_local, omega)

        # get displacement
        N = self.axialforce(x, omega, C)

        return [N]
    def LocalElementStresses(self, u_nodes_local, omega, num_points):
        """
        Calcualte the local displacements u(s).

        Input:
            u_nodes_global: array. The nodes in global coordinates
            omega:          array. Range of frequencies of analysis
            num_points:     value. Number of points to divide the element in.
        Output:
            u: array. Amplitude of vertical displacement
        """
        # get local axis to evaluate on
        L = self.L
        x = np.linspace(0.0, L, num_points)

        # calculate coeficients
        C = self.Coefficients(u_nodes_local, omega)

        # get displacement
        sigma = self.axialstress(x, omega, C)

        return [sigma]

    def Coefficients(self, u_nodes_local, omega):
        """
        Calculate the coefficients of the general solution, in this case 2.

        Input:
            u_node_local: local degrees of freedom
            omega:        array. Range of frequencies of analysis
        Output:
            C:            array (2), coefficients of general solution (C1, C2)
        """
        # read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        alpha_1, alpha_2, = self.ElementWaveNumbers(omega)

        # get distributed load value
        q = self.q[0]

        # calculate the coefficients with A_mat copy-pasted from Maple
        A_mat = np.array([[1, 1],
                          [np.exp(-1j*alpha_1*L),
                           np.exp(-1j*alpha_2*L)]])

        # TODO - check
        u_load = np.array([-q/(rho*A*omega**2), -q/(rho*A*omega**2)])
        C = np.linalg.inv(A_mat) @ (u_nodes_local - u_load)

        # return the coefficients, in this case 2
        return C

    def displacement(self, x, omega, C=None, u_node_local=None):
        """
        Get the transverse displacments of the 1D rod.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            u: array. horizontal displacements
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        
        # read all the variables
        rho = self.rho
        A = self.A
        
        # get the wavenumbers
        alpha_1, alpha_2 = self.ElementWaveNumbers(omega)

        # check if C is input
        if C is None:
            C = self.Coefficients(u_node_local, omega)
            
        # get distributed load value
        q = self.q[0]
            
        # displacements
        u_load = np.array([-q/(rho*A*omega**2)])
        u = C[0]*np.exp(-1j*alpha_1*x) + C[1]*np.exp(-1j*alpha_2*x) + u_load

        return u
    
    def axialforce(self, x, omega, C=None, u_node_local=None):
        """
        Get the transverse displacments of the 1D rod.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            u: array. horizontal displacements
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        
        # read all the variables
        # rho = self.rho
        E = self.E
        A = self.A
        
        # get the wavenumbers
        alpha_1, alpha_2 = self.ElementWaveNumbers(omega)

        # check if C is input
        if C is None:
            C = self.Coefficients(u_node_local, omega)
            
        # get distributed load value
        # q = self.q[0]
            
        # displacements
        # u_load = np.array([-q/(rho*A*omega**2)])
        # u = C[0]*np.exp(-1j*alpha_1*x) + C[1]*np.exp(-1j*alpha_2*x) + u_load
        dudx = C[0]*np.exp(-1j*alpha_1*x)*(-1j*alpha_1)+C[1]*np.exp(-1j*alpha_2*x)*(-1j*alpha_2)
        axialforce = E*A*dudx

        return axialforce
    
    def axialstress(self, x, omega, C=None, u_node_local=None):
        """
        Get the transverse displacments of the 1D rod.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            u: array. horizontal displacements
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        
        # read all the variables
        # rho = self.rho
        E = self.E
        
        # get the wavenumbers
        alpha_1, alpha_2 = self.ElementWaveNumbers(omega)

        # check if C is input
        if C is None:
            C = self.Coefficients(u_node_local, omega)
            
        # get distributed load value
        # q = self.q[0]
            
        # displacements
        # u_load = np.array([-q/(rho*A*omega**2)])
        # u = C[0]*np.exp(-1j*alpha_1*x) + C[1]*np.exp(-1j*alpha_2*x) + u_load
        dudx = C[0]*np.exp(-1j*alpha_1*x)*(-1j*alpha_1)+C[1]*np.exp(-1j*alpha_2*x)*(-1j*alpha_2)
        axialstress = E*dudx

        return axialstress
