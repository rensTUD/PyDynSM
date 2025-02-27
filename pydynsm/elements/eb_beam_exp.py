# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:27:06 2024

@author: rensv
"""
# import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition


@ElementFactory.ElementType('EulerBernoulli Beam')#, required_parameters = ['rho', 'A', 'E', 'Ib', 'L'])
class EulerBernoulliBeam(StructuralElement):
    """Class for Euler-Bernoulli beam element."""

    def __init__(self, rho, A, E, Ib, Wb, L, ksi=None):
        """
        Initialise EulerBernoulliBeam class.

        Input:
            rho : value. Density of element [kg/m^3]
            A   : value. Area of cross section [m^2]
            E   : value. Young's modulus [Pa]
            Ib  : value. Cross-sectional inertia [m^4]
            L   : value. Length of the element [m]
            ksi : value. Material damping [-], default: none
        """
        # define the of dofs the EulerBernoulli beam and initialise
        dofs = ['z','phi_y']
        super().__init__(dofs)

        # Initialise local beam element with necessary parameters
        self.rho = rho
        self.A = A
        # self.E = E # please see the line below for complex E value
        self.Ib = Ib
        self.L = L
        # Assign ksi if given,  otherwise assign a default value
        self.ksi = ksi if ksi is not None else 0.01
        self.E = E*(1+2j*self.ksi)
        self.Wb = Wb 
        # set q standard to 0
        self.q = np.zeros(len(dofs))

    def LocalStiffness(self, omega):
        """
        Determine the stiffness of the EB-beam.

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
        """
        # Assign local variables for ease of coding
        E = self.E
        Ib = self.Ib
        L = self.L
        # Obtain the wavenumbers
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # Initialise K_local
        K_local = np.empty((4, 4), complex)
        # Copy-paste the matrix directly from Maple (also called K_dyn)
        K_local = np.array([[-1j*Ib*E*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4))*alpha_1)/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             Ib*(((alpha_3**3 - alpha_4**3)*alpha_1 + (-alpha_3**3 + alpha_4**3)*alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) + ((-alpha_2**3 + alpha_4**3)*alpha_1 + alpha_3*alpha_2**3 - alpha_4**3*alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + ((alpha_2**3 - alpha_3**3)*alpha_1 - alpha_4*alpha_2**3 + alpha_4*alpha_3**3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + ((alpha_2 - alpha_3)*alpha_1**3 - alpha_4**3*(alpha_2 - alpha_3))*np.exp(-1j*L*(alpha_2 + alpha_3)) + ((-alpha_2 + alpha_4)*alpha_1**3 + alpha_3**3*(alpha_2 - alpha_4))*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_1 - alpha_2)*(alpha_1**2 + alpha_2*alpha_1 + alpha_2**2)*(alpha_3 - alpha_4))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*Ib*E*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -Ib*((alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [Ib*E*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*Ib*((alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -Ib*E*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -1j*Ib*((alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + ((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3))*(alpha_1 - alpha_2))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [1j*Ib*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4))*(alpha_3 - alpha_4))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + ((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4))*(alpha_3 - alpha_4))*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*Ib*E*(alpha_1*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*alpha_4)/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -Ib*E*(((alpha_3 - alpha_4)*alpha_1**3 - alpha_2**3*(alpha_3 - alpha_4))*np.exp(-1j*L*(alpha_1 + alpha_2)) + ((-alpha_2 + alpha_4)*alpha_1**3 + alpha_3**3*(alpha_2 - alpha_4))*np.exp(-1j*L*(alpha_1 + alpha_3)) + ((alpha_2 - alpha_3)*alpha_1**3 - alpha_4**3*(alpha_2 - alpha_3))*np.exp(-1j*L*(alpha_1 + alpha_4)) + ((alpha_2**3 - alpha_3**3)*alpha_1 - alpha_4*alpha_2**3 + alpha_4*alpha_3**3)*np.exp(-1j*L*(alpha_2 + alpha_3)) + ((-alpha_2**3 + alpha_4**3)*alpha_1 + alpha_3*alpha_2**3 - alpha_4**3*alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_3**2 + alpha_3*alpha_4 + alpha_4**2)*(alpha_1 - alpha_2))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [-Ib*E*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3))*(alpha_3 - alpha_4))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -1j*Ib*((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + ((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3))*(alpha_3 - alpha_4))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -Ib*E*(alpha_1*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_4*(alpha_2*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             -1j*Ib*((alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))]])

        return K_local

    def ElementWaveNumbers(self, omega):
        """
        Determine the wavenumbers of Euler-Bernoulli beam.

        Input:
            omega: array. Range of frequencies of analysis.
        Output:
            alpha_1, alpha_2, alpha_3, alpha_4: values. Wavenumbers
        """
        rho = self.rho
        A = self.A
        E = self.E
        Ib = self.Ib
        # Copy-paste from Maple document of Wavenumbers
        # 4th-order derivative in x, so 4 wavenumbers
        alpha_1 = 1/E/Ib*(A*omega**2*rho*E**3*Ib**3)**(1/4)
        alpha_2 = 1j/E/Ib*(A*omega**2*rho*E**3*Ib**3)**(1/4)
        alpha_3 = -1/E/Ib*(A*omega**2*rho*E**3*Ib**3)**(1/4)
        alpha_4 = -1j/E/Ib*(A*omega**2*rho*E**3*Ib**3)**(1/4)

        return alpha_1, alpha_2, alpha_3, alpha_4

    def LocalDistributedLoad(self, q, omega):
        """
        Add a distributed load to the local element.

        Input:
            q: array. Distributed load. With a shape such as:
            q = [q_z, q_phi]
            omega: array. Range of frequencies of analyis
        Output:
            el: array. Force vector
        """
        # assign load to itself to keep track
        self.q = q

        # assign local variables for ease of coding
        rho = self.rho
        A = self.A
        E = self.E
        Ib = self.Ib
        L = self.L

        # determine wavenumber
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # extract loads
        q_z = q[0]
        q_phi = q[1]

        # TODO - check for correctness
        el = np.array([-1j*E*Ib*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A + 1j*E*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))*Ib/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A,
                       E*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))*Ib/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A - E*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)))*Ib/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A,
                       1j*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4))*(alpha_3 - alpha_4))*E*Ib/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A + 1j*E*(alpha_1*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_4*(alpha_2*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)))*Ib/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A,
                       -E*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3))*(alpha_3 - alpha_4))*Ib/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A - E*(alpha_1*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*alpha_4)*Ib/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A])

        return el

    def LocalElementDisplacements(self, u_nodes_local, omega, num_points):
        """
        Calculate local displacements w(s) and rotational displacement phi(s).

        Input:
            u_nodes_global: array. The nodes in global coordinates
            omega: array. Range of frequencies of analysis
            num_points: value. Number of points to divide the element in.
        Output:
            w: array. Amplitude of vertical displacement
            phi: array. Amplitude of rotational displacement
        """
        # get local axis to evaluate on
        L = self.L
        x = np.linspace(0.0, L, num_points)

        # calculate coeficients
        C = self.Coefficients(u_nodes_local, omega)

        # get displacement
        w = self.displacement(x, omega, C)

        # get rotations
        phi = self.rotation(x, omega, C)

        return [w, phi]
    
    def LocalElementForces(self,u_nodes_local,omega,num_points):
        """
        Calculate local Element shear force V(s) and moment M(s).

        Input:
            u_nodes_global: array. The nodes in global coordinates
            omega: array. Range of frequencies of analysis
            num_points: value. Number of points to divide the element in.
        Output:
            V: array. Amplitude of elemental shear force
            M: array. Amplitude of elemental bending moment
        """
        # get local axis to evaluate on
        L = self.L
        x = np.linspace(0.0, L, num_points)
        
        # calculate coeficients
        C = self.Coefficients(u_nodes_local, omega)
        
        # get shear force
        V = self.shearforce(x, omega, C)
        # get bending moment
        M = self.moment(x, omega, C)
        
        return [V, M]
    
    def LocalElementStresses(self,u_nodes_local,omega,num_points):
        """
        Calculate local Element shear force V(s) and moment M(s).

        Input:
            u_nodes_global: array. The nodes in global coordinates
            omega: array. Range of frequencies of analysis
            num_points: value. Number of points to divide the element in.
        Output:
            V: array. Amplitude of elemental shear force
            M: array. Amplitude of elemental bending moment
        """
        # get local axis to evaluate on
        L = self.L
        x = np.linspace(0.0, L, num_points)
        
        # calculate coeficients
        C = self.Coefficients(u_nodes_local, omega)
        
        # get shear force
        tau = self.shearstress(x, omega, C)
        # get bending moment
        sigma_yy = self.bendingstress(x, omega, C)
        
        return [tau, sigma_yy]

    def Coefficients(self, u_nodes_local, omega):
        """
        Calculate the coefficients of the general solution, in this case 4.

        Input:
            u_node_local: local degrees of freedom
            omega: array. Range of frequencies of analysis
        Output:
            C: array (4), coefficients of general solution (C1, C2, C3, C4)
        """
        # Read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # get distributed load value
        # extract loads
        q_z = self.q[0]
        q_phi = self.q[1] 
        # should be like:
        # q_b, q_m = self.

        # calculate the coefficients with A_mat copy-pasted from Maple
        A_mat = np.array([[1, 1, 1, 1],
                          [1j*alpha_1, 1j*alpha_2, 1j*alpha_3, 1j*alpha_4],
                          [np.exp(-1j*alpha_1*L), np.exp(-1j*alpha_2*L),
                           np.exp(-1j*alpha_3*L), np.exp(-1j*alpha_4*L)],
                          [1j*np.exp(-1j*alpha_1*L)*alpha_1,
                           1j*np.exp(-1j*alpha_2*L)*alpha_2,
                           1j*np.exp(-1j*alpha_3*L)*alpha_3,
                           1j*np.exp(-1j*alpha_4*L)*alpha_4]])

        # TODO - check for correctness
        u_load = np.array([q_z/omega**2/rho/A, 0, q_z/omega**2/rho/A, 0])
        C = np.linalg.inv(A_mat) @ (u_nodes_local - u_load)
        # + np.array([1/(E*Ib*beta_b**4),0,1/(E*Ib*beta_b**4),0]) * q)

        return C

    def displacement(self, x, omega, C=None, u_node_local=None):
        """
        Get the transverse displacments of the EB beam.

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
        # Read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)
        
        # get distributed load value
        # extract loads
        q_z = self.q[0]
        q_phi = self.q[1] 
        # should be like:
        # q_b, q_m = self.
        
        # Check if C is input, otherwise, calculate it
        if C is None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        u_load = np.array([q_z/omega**2/rho/A])
        w = (C[0]*np.exp(-1j*alpha_1*x) +
             C[1]*np.exp(-1j*alpha_2*x) +
             C[2]*np.exp(-1j*alpha_3*x) +
             C[3]*np.exp(-1j*alpha_4*x)) + u_load

        return w

    def rotation(self, x, omega, C=None, u_node_local=None):
        """
        Get the rotations of the EB beam.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            phi: array. Rotational displacements
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        # Read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)
        
        # get distributed load value
        # extract loads
        q_z = self.q[0]
        q_phi = self.q[1] 
        # should be like:
        # q_b, q_m = self.

        # check if C is input
        if C is None:
            C = self.Coefficients(u_node_local, omega)

        # Rotations
        u_load = 0
        phi = (-1j*C[0]*alpha_1*np.exp(-1j*alpha_1*x) -
               1j*C[1]*alpha_2*np.exp(-1j*alpha_2*x) -
               1j*C[2]*alpha_3*np.exp(-1j*alpha_3*x) -
               1j*C[3]*alpha_4*np.exp(-1j*alpha_4*x)) + u_load

        return -phi
    
    def moment(self, x, omega, C=None, u_node_local=None):
        """
        Get the moment field of the EB beam.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            mmt: array. moment field of the given element
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        # Read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        E = self.E
        I = self.Ib
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)
        
        if C is None:
            C = self.Coefficients(u_node_local, omega)
        
        kappa = ((-1j*alpha_1)**2*C[0]*np.exp(-1j*alpha_1*x) +
                (-1j*alpha_2)**2*C[1]*np.exp(-1j*alpha_2*x) +
                (-1j*alpha_3)**2*C[2]*np.exp(-1j*alpha_3*x) +
                (-1j*alpha_4)**2*C[3]*np.exp(-1j*alpha_4*x))
        
        mmt = -E*I*kappa
        return mmt
    
    def bendingstress(self, x, omega, C=None, u_node_local=None):
        """
        Get the moment field of the EB beam.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            mmt: array. moment field of the given element
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        # Read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        E = self.E
        I = self.Ib
        W = self.Wb
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)
        
        if C is None:
            C = self.Coefficients(u_node_local, omega)
        
        kappa = ((-1j*alpha_1)**2*C[0]*np.exp(-1j*alpha_1*x) +
                (-1j*alpha_2)**2*C[1]*np.exp(-1j*alpha_2*x) +
                (-1j*alpha_3)**2*C[2]*np.exp(-1j*alpha_3*x) +
                (-1j*alpha_4)**2*C[3]*np.exp(-1j*alpha_4*x))
        
        bs = -E*I*kappa/W
        return bs
    
    def shearforce(self, x, omega, C=None, u_node_local=None):
        """
        Get the shear force field of the EB beam.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            sf: array. moment field of the given element
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        # Read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        E = self.E
        I = self.Ib
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)
        
        if C is None:
            C = self.Coefficients(u_node_local, omega)
            
        dudx3 = ((-1j*alpha_1)**3*C[0]*np.exp(-1j*alpha_1*x) +
                (-1j*alpha_2)**3*C[1]*np.exp(-1j*alpha_2*x) +
                (-1j*alpha_3)**3*C[2]*np.exp(-1j*alpha_3*x) +
                (-1j*alpha_4)**3*C[3]*np.exp(-1j*alpha_4*x))
        
        sf = -E*I*dudx3
        return sf
    
    def shearstress(self, x, omega, C=None, u_node_local=None):
        """
        Get the shear force field of the EB beam.

        Input:
            x: array. Points along element
            omega: array. Range of frequencies of analysis
            C: array. Values of coefficients of general solution
            u_node_local: local nodes
        Ouput:
            sf: array. moment field of the given element
        Note:
            if C is not given, then calculate it based on u_node_local.
        """
        # Read all the variables
        rho = self.rho
        A = self.A
        L = self.L
        E = self.E
        I = self.Ib
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)
        
        if C is None:
            C = self.Coefficients(u_node_local, omega)
            
        dudx3 = ((-1j*alpha_1)**3*C[0]*np.exp(-1j*alpha_1*x) +
                (-1j*alpha_2)**3*C[1]*np.exp(-1j*alpha_2*x) +
                (-1j*alpha_3)**3*C[2]*np.exp(-1j*alpha_3*x) +
                (-1j*alpha_4)**3*C[3]*np.exp(-1j*alpha_4*x))
        
        ss = -E*I*dudx3/self.A
        return ss
    
