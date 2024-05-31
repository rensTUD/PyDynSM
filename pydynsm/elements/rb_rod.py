# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:27:06 2024

@author: rensv
"""
# import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition


@ElementFactory.ElementType('Rayleigh-Bishop Rod')
class RayleighBishopRod(StructuralElement):

    element_name = 'Rayleigh-Bishop Rod'

    def __init__(self, rho, A, E, Ir, L, nu, G, ksi=None):
        """
        Input:
            rho: value. Density of the element's material [kg/m^3]
            A:   value. Area of the element [m^2]
            E:   value. E-modulus of the element's material [Pa]
            Ir:  value. Cross-sectional moment of inertia of the element [m^4]
            L:   value. Length of element [m]
            nu:  value. Poisson's ratio [-]
            ksi: value. Damping of the element's material [-]
        """
        # define what dofs the RB rod contributes to and initialise
        # TODO - WHAT ARE THE DEGREES OF FREEDOM, maybe degree of freedom 4?
        dofs = [0]
        super().__init__(dofs)

        # Initialise local rod element with necessary parameters
        self.rho = rho
        self.A = A
        self.E = E
        self.Ir = Ir
        self.L = L
        self.nu = nu
        self.G = G
        # assisgn ksi if given otherwise assign a default value
        self.ksi = ksi if ksi is not None else 0.01

    def local_stiffness(self, omega):
        '''
        Determines the stiffness of the RB-rod.
        As its 2D, stiffness will be 4x4 matrix as:
            [V_left, M_left, V_right, M_right] =
            K.[w_left, phi_left, w_right, phi_right]
        where:
            K = [K_V_ll, K_M_ll, K_V_lr, K_M_lr;
                 K_V_rl, K_M_rl, K_V_rr, K_M_rr]
        Input:
            omega: array. Range of frequencies of analysis
        Output:
            K_local: matrix. Dynamic stiffness matrix (also K_dyn)
        '''

        # Assign local variables for ease of coding
        rho = self.rho
        A = self.A
        E = self.E
        Ir = self.Ir
        L = self.L
        nu = self.nu
        G = self.G
        # Obtain wavenumbers
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # Initialise K_local
        K_local = np.empty((4, 4), complex)
        # Copy-paste the matrix directly from Maple (also called K_dyn)
        K_local = np.array([[-1j*(-alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - alpha_1*(alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)))*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             ((alpha_1 - alpha_2)*(alpha_3 - alpha_4)*(nu**2*G*Ir*alpha_3**2 + G*Ir*alpha_3*alpha_4*nu**2 + G*Ir*alpha_4**2*nu**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*(nu**2*G*Ir*alpha_2**2 + G*Ir*alpha_2*alpha_4*nu**2 + G*Ir*alpha_4**2*nu**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*(nu**2*G*Ir*alpha_2**2 + G*Ir*alpha_2*alpha_3*nu**2 + nu**2*G*Ir*alpha_3**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*(nu**2*G*Ir*alpha_1**2 + G*Ir*alpha_1*alpha_4*nu**2 + G*Ir*alpha_4**2*nu**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*(nu**2*G*Ir*alpha_1**2 + G*Ir*alpha_1*alpha_3*nu**2 + nu**2*G*Ir*alpha_3**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_1 - alpha_2)*(alpha_3 - alpha_4)*(nu**2*G*Ir*alpha_1**2 + G*Ir*alpha_1*alpha_2*nu**2 + nu**2*G*Ir*alpha_2**2 - Ir*nu**2*omega**2*rho + A*E))/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -1j*(-alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) + alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) - ((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*alpha_3*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3))*(alpha_1 - alpha_2))*nu**2*Ir*G/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -nu**2*Ir*(-(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) + (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) - ((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3))*(alpha_1 - alpha_2))*G/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [-(-alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - alpha_1*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*(-(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4))*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -(-alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) + alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) - ((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*alpha_3*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3))*(alpha_1 - alpha_2))*nu**2*Ir*G/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*(-(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) + (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) - ((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3))*(alpha_1 - alpha_2))*nu**2*Ir*G/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [-1j*nu**2*Ir*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)))*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)))*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*(-alpha_1*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_1*alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_1*alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + ((alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*alpha_2*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*alpha_4)*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             (-(alpha_1 - alpha_2)*(alpha_3 - alpha_4)*(nu**2*G*Ir*alpha_1**2 + G*Ir*alpha_1*alpha_2*nu**2 + nu**2*G*Ir*alpha_2**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*(nu**2*G*Ir*alpha_1**2 + G*Ir*alpha_1*alpha_3*nu**2 + nu**2*G*Ir*alpha_3**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*(nu**2*G*Ir*alpha_1**2 + G*Ir*alpha_1*alpha_4*nu**2 + G*Ir*alpha_4**2*nu**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*(nu**2*G*Ir*alpha_2**2 + G*Ir*alpha_2*alpha_3*nu**2 + nu**2*G*Ir*alpha_3**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*(nu**2*G*Ir*alpha_2**2 + G*Ir*alpha_2*alpha_4*nu**2 + G*Ir*alpha_4**2*nu**2 - Ir*nu**2*omega**2*rho + A*E)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_1 - alpha_2)*(alpha_3 - alpha_4)*(nu**2*G*Ir*alpha_3**2 + G*Ir*alpha_3*alpha_4*nu**2 + G*Ir*alpha_4**2*nu**2 - Ir*nu**2*omega**2*rho + A*E))/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (-alpha_1 + alpha_4)*(-alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (-alpha_2 + alpha_4)*(-alpha_1 + alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [-nu**2*Ir*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)))*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)))*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             (-alpha_1*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_1*alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_1*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + alpha_4*(alpha_2*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -1j*(-(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*nu**2*Ir*G/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))]])

        return K_local

    def element_wavenumbers(self, omega):
        '''
        Determines the wavenumbers of Rayleigh-Bishop rod
        Input:
            omega: array. Range of frequencies of analysis.
        Output:
            alpha_1, alpha_2, alpha_3, alpha_4: values. Wavenumbers
        '''
        rho = self.rho
        A = self.A
        E = self.E
        Ir = self.Ir
        nu = self.nu
        G = self.G

        alpha_1 = 1/G/Ir/nu*np.sqrt(-2*G*Ir*(-Ir*nu**2*omega**2*rho + A*E - np.sqrt(Ir**2*nu**4*omega**4*rho**2 - 2*A*E*Ir*nu**2*omega**2*rho + 4*A*G*Ir*nu**2*omega**2*rho + A**2*E**2)))/2
        alpha_2 = -1/G/Ir/nu*np.sqrt(-2*G*Ir*(-Ir*nu**2*omega**2*rho + A*E - np.sqrt(Ir**2*nu**4*omega**4*rho**2 - 2*A*E*Ir*nu**2*omega**2*rho + 4*A*G*Ir*nu**2*omega**2*rho + A**2*E**2)))/2
        alpha_3 = 1/G/Ir/nu*np.sqrt(-2*G*Ir*(-Ir*nu**2*omega**2*rho + A*E + np.sqrt(Ir**2*nu**4*omega**4*rho**2 - 2*A*E*Ir*nu**2*omega**2*rho + 4*A*G*Ir*nu**2*omega**2*rho + A**2*E**2)))/2
        alpha_4 = -1/G/Ir/nu*np.sqrt(-2*G*Ir*(-Ir*nu**2*omega**2*rho + A*E + np.sqrt(Ir**2*nu**4*omega**4*rho**2 - 2*A*E*Ir*nu**2*omega**2*rho + 4*A*G*Ir*nu**2*omega**2*rho + A**2*E**2)))/2

        return alpha_1, alpha_2, alpha_3, alpha_4

    def local_distributed_load(self, q, omega):
        '''
        Add a distributed load to the local element
        Input:
            q: array. Distributed load. With a shape such as:
            q = [q_z, q_theta] (kind of torsion)
            omega: array. Range of frequencies of analyis
        Output:
        '''
        # assign load to itself to keep track
        self.q = q

        # assign local variables for ease of coding
        rho = self.rho
        A = self.A
        Ir = self.Ir
        L = self.L
        nu = self.nu
        G = self.G

        # determine wavenumber
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # extract loads
        q_x = q[0]
        q_theta = q[1]

        # TODO - check for correctness
        el = np.array([-1j*G*(-alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4))*alpha_1)*nu**2*Ir/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_x/rho/A/omega**2 + -1j*G*(-alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) + alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) - (alpha_1 - alpha_2)*((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*alpha_3*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))*nu**2*Ir/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_x/rho/A/omega**2,
                       -G*(-alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*alpha_1)*nu**2*Ir/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_theta/rho/A/omega**2 - G*nu**2*Ir*(-alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) + alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) - (alpha_1 - alpha_2)*((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*alpha_3*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_theta/rho/A/omega**2,
                       -1j*G*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4))*(alpha_3 - alpha_4))*nu**2*Ir/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_x/rho/A/omega**2 + 1j*G*nu**2*(-alpha_1*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_1*alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_1*alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + alpha_4*((alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*alpha_2*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)))*Ir/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_x/rho/A/omega**2,
                       -G*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3))*(alpha_3 - alpha_4))*nu**2*Ir/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_theta/rho/A/omega**2 + G*nu**2*Ir*(-alpha_1*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + alpha_1*alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - alpha_1*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + alpha_4*(alpha_2*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))/(-(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) - (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) + (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_theta/rho/A/omega**2])

        return el

    def local_element_displacements(self, u_nodes_global, omega, num_points):
        '''
        This function calculates the coefficients C, local displacements w(s)
        and rotational displacement phi(s).
        Input:
            u_nodes_global: array. The nodes in global coordinates
            omega: array. Range of frequencies of analysis
            num_points: value. Number of points to divide the element in.
        Output:
            w: array. Amplitude of vertical displacement
            phi: array. Amplitude of rotational displacement
        '''

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

        # get rotations
        phi = self.rotation(x, omega, C)

        return [u, phi]

    def coefficients(self, u_node_local, omega):
        '''
        Calculates the coefficients of the general solution, in this case 4
        Input:
            u_node_local: local degrees of freedom
            omega: array. Range of frequencies of analysis
        Output:
            C: array (4), coefficients of general solution (C1, C2, C3, C4)
        '''
        # read all the variables
        rho = self.rho
        A = self.A
        E = self.E
        Ir = self.Ir
        L = self.L
        nu = self.nu
        G = self.G
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # get distributed load value
        q = self.q[1]
        # should be like:
        # q_b, q_m = self.

        # calculate the coefficients with A_mat copy-pasted from Maple
        A_mat = np.array([[1, 1, 1, 1],
                          [-1j*alpha_1, -1j*alpha_2, -1j*alpha_3, -1j*alpha_4],
                          [np.exp(-1j*alpha_1*L), np.exp(-1j*alpha_2*L),
                           np.exp(-1j*alpha_3*L), np.exp(-1j*alpha_4*L)],
                          [-1j*np.exp(-1j*alpha_1*L)*alpha_1,
                           -1j*np.exp(-1j*alpha_2*L)*alpha_2,
                           -1j*np.exp(-1j*alpha_3*L)*alpha_3,
                           -1j*np.exp(-1j*alpha_4*L)*alpha_4]])

        # TODO - check if this is correct
        u_part = np.array([q/rho/A/omega**2, 0, q/rho/A/omega**2, 0])
        C = np.linalg.inv(A_mat) @ (u_node_local + u_part)

        return C

    def displacement(self, x, omega, C=None, u_node_local=None):
        """
        Gets the transverse displacments of the Rayleigh-Bishop rod
        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            w: array. Transverse displacements
        Note:
            if C is not given, then calculate it based on u_node_local.
        """

        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # check if C is input
        if C == None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        u = (C[0]*np.exp(-1j*alpha_1*x) +
             C[1]*np.exp(-1j*alpha_2*x) +
             C[2]*np.exp(-1j*alpha_3*x) +
             C[3]*np.exp(-1j*alpha_4*x))

        return u

    def rotation(self, x, omega, C=None, u_node_local=None):
        '''
        Gets the rotations of the Rayleigh-Bishop rod
        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            w: array. Transverse displacements
        Note:
            if C is not given, then calculate it based on u_node_local.
        '''

        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # check if C is input
        if C == None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        phi = (-1j*C[0]*alpha_1*np.exp(-1j*alpha_1*x) -
               1j*C[1]*alpha_2*np.exp(-1j*alpha_2*x) -
               1j*C[2]*alpha_3*np.exp(-1j*alpha_3*x) -
               1j*C[3]*alpha_4*np.exp(-1j*alpha_4*x))

        return -phi
