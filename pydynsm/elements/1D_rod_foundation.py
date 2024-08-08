# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:27:16 2024

@author: rensv
"""

# import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition


@ElementFactory.ElementType('Rod over viscoelastic foundation')
class Rod1D_foundation(StructuralElement):
    """Class for 1D rod element over viscoelastic foundation."""

    element_name = 'Rod over viscoelastic foundation'

    def __init__(self, rho, A, E, L, k, ksi=None):
        """
        Initialise the Rod1D_foundation class.

        Input:
            rho: value. Density of the element's material [kg/m^3]
            A:   value. Area of the element [m^2]
            E:   value. E-modulus of the element's material [Pa]
            L:   value. Length of element [m]
            k:   value. Soil stiffness [N/m^2]
            ksi: value. Damping of the element's material [-]
        """
        # define what dofs the rod contributes to and initialise
        dofs = [0]
        super().__init__(dofs)

        # Initialise local rod element with necessary parameters
        self.rho = rho
        self.A = A
        self.E = E
        self.L = L
        self.k = k
        # assisgn ksi if given otherwise assign a default value
        self.ksi = ksi if ksi is not None else 0.01

    def local_stiffness(self, omega):
        """
        Determine the stiffness of the rod.

        As its 1D, stiffness will be 2x2 matrix as:
            [F_left, F_right] = K.[u_left, u_right]
        where: ASK THIIIIIS
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

        # Initialize K matrix, 2x2 matrix
        K_local = np.empty((2, 2), complex)
        # Copy-paste the matrix directly from Maple (also called K_dyn)
        K_local = np.array([[1j*A*E*(np.exp(-1j*alpha_2*L)*alpha_1 - np.exp(-1j*alpha_1*L)*alpha_2)/(np.exp(-1j*alpha_2*L) - np.exp(-1j*alpha_1*L)),
                             1j*A*E*(alpha_1 - alpha_2)/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))]
                            [1j*A*E*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2))/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))
                             -1j*A*E*(np.exp(-1j*alpha_2*L)*alpha_2 - np.exp(-1j*alpha_1*L)*alpha_1)/(np.exp(-1j*alpha_2*L) - np.exp(-1j*alpha_1*L))]])

        return K_local

    def element_wavenumbers(self, omega):
        """
        Determine the wavenumbers of 1D rod.

        Input:
            omega: array. Range of frequencies of analysis.
        Output:
            alpha_1, alpha_2,: values. Wavenumbers
        """
        rho = self.rho
        E = self.E
        A = self.A
        k = self.k

        # Copy-paste from Maple document of Wavenumbers
        # 2nd-order derivative in x, so 2 wavenumbers
        alpha_1 = np.sqrt((rho*A*omega**2 + k)/(A*E))
        alpha_2 = - np.sqrt((rho*A*omega**2 + k)/(A*E))

        return alpha_1, alpha_2

    def local_distributed_load(self, q, omega):
        """
        Add a distributed load to the local element.

        Input:
            q: array. Distributed load. With a shape such as:
            q = [q_z, q_phi]
            omega: array. Range of frequencies of analyis
        Output:
            el:    array. Force vector
        """
        # assign load to itself to keep track
        self.q = q

        # assign local variables for ease of coding
        rho = self.rho
        E = self.E
        A = self.A
        L = self.L
        k = self.k

        # extract loads
        q_x = q[0]

        # determine wavenumber
        alpha_1, alpha_2 = self.ElementWaveNumbers(omega)

        # Force matrix, from Maple (called F_q)
        el = np.array([1j*A*E*(np.exp(-1j*alpha_1*L)*alpha_2 - np.exp(-1j*alpha_2*L)*alpha_1)/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))*q_x/(rho*A*omega**2 + k)
                       + 1j*A*E*(alpha_1 - alpha_2)/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))*q_x/(rho*A*omega**2 + k),
                       1j*A*E*(alpha_1 - alpha_2)/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))*np.exp(-1j*L*(alpha_1 + alpha_2))*q_x/(rho*A*omega**2 + k)
                       + -1j*A*E*(np.exp(-1j*alpha_1*L)*alpha_1 - np.exp(-1j*alpha_2*L)*alpha_2)/(-np.exp(-1j*alpha_2*L) + np.exp(-1j*alpha_1*L))*q_x/(rho*A*omega**2 + k)])

        return el

    def local_element_displacements(self, u_nodes_global, omega, num_points):
        """
        Calculate local displacements u(s).

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

        # determine nodal displacement in local axis
        u_node_local = self.R @ u_nodes_global

        # extract only the needed displacements
        u_node_local = u_node_local(self.dofs)

        # calculate coeficients
        C = self.Coefficients(u_node_local, omega)

        # get displacement
        u = self.displacement(x, omega, C)

        return [u]

    def coefficients(self, u_node_local, omega):
        """
        Calculate the coefficients of the general solution, in this case 2.

        Input:
            u_node_local: local degrees of freedom
            omega: array. Range of frequencies of analysis
        Output:
            C: array (2), coefficients of general solution (C1, C2)
        """
        # read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        k = self.k
        alpha_1, alpha_2, = self.ElementWaveNumbers(omega)

        # get distributed load value
        q = self.q[0]

        # calculate the coefficients with A_mat copy-pasted from Maple
        A_mat = np.array([[1, 1],
                          [np.exp(-1j*alpha_1*L),
                           np.exp(-1j*alpha_2*L)]])

        # from Maple u_part.
        u_load = np.array([q/(rho*A*omega**2 + k), q/(rho*A*omega**2 + k)])
        # Calculate C using local nodes and A matrix
        C = np.linalg.inv(A_mat) @ (u_node_local + u_load)

        return C

    def displacement(self, x, omega, C=None, u_node_local=None):
        """
        Get the transverse displacments of the 1D rod with viscoel foundation.

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
        alpha_1, alpha_2 = self.ElementWaveNumbers(omega)

        # check if C is input
        if C is None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        u = C[0]*np.exp(-1j*alpha_1*x) + C[1]*np.exp(-1j*alpha_2*x)

        return u
