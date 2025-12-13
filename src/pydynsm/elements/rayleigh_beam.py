# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:27:06 2024

@author: rensv
"""
# import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition


@ElementFactory.ElementType('Rayleigh Beam')
class RayleighBeam(StructuralElement):
    """Class for Rayleigh beam."""

    element_name = 'Rayleigh Beam'

    def __init__(self, section, rho, E, L, ksi=None):
        """
        Initialise RayleighBeam class.

        Parameters
        ----------
        section : Section
            Section object providing geometric properties (A, I_y)
        rho : float
            Density of the element's material [kg/m^3]
        E : float
            E-modulus of the element's material [Pa]
        L : float
            Length of element [m]
        ksi : float, optional
            Damping of the element's material [-], default: 0.01
        """
        # define what dofs the rayleigh beam contributes to and initialise
        # the dofs are w and phi, so 1 and 2
        dofs = ['z','phi_y']
        super().__init__(dofs)

        # Extract geometric properties from section
        self.A = section.A
        self.Ib = section.I_y  # For 2D x-z plane bending

        # Initialise local beam element with necessary parameters
        self.rho = rho
        self.E = E
        self.L = L
        # assisgn ksi if given otherwise assign a default value
        self.ksi = ksi if ksi is not None else 0.01

    def LocalStiffness(self, omega):
        """
        Determine the stiffness of the beam.

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
        rho = self.rho
        E = self.E
        Ib = self.Ib
        L = self.L
        # Calculate wavenumber with corresponding function
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # Initialize K_local
        K_local = np.empty((4, 4), complex)
        # Copy-paste the matrix directly from Maple (also called K_dyn)
        K_local = np.array([[-1j*Ib*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             Ib*((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_3 ** 2 + E*alpha_3*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_2 ** 2 + E*alpha_2*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_2 ** 2 + E*alpha_2*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_1 ** 2 + E*alpha_1*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_1 ** 2 + E*alpha_1*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_1 ** 2 + E*alpha_1*alpha_2 + E*alpha_2 ** 2 + omega ** 2*rho))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*Ib*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -Ib*((alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [Ib*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*((alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4))*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -Ib*E*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -1j*Ib*((alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + ((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3))*(alpha_1 - alpha_2))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [1j*Ib*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)))*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                             1j*(alpha_1*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*alpha_4)*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -Ib*((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_1 ** 2 + E*alpha_1*alpha_2 + E*alpha_2 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_1 ** 2 + E*alpha_1*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_1 ** 2 + E*alpha_1*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_2 ** 2 + E*alpha_2*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_2 ** 2 + E*alpha_2*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_3 ** 2 + E*alpha_3*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))],
                            [-Ib*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -1j*Ib*((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -Ib*E*(alpha_1*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*alpha_4)/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)),
                            -1j*((alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))]])

        return K_local

    def ElementWaveNumbers(self, omega):
        """
        Determine the wavenumbers of Rayleigh beam.

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
        alpha_1 = 1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega + np.sqrt(Ib**2*rho**2*omega**2 + 4*E*Ib*A*rho))*omega)/2
        alpha_2 = -1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega + np.sqrt(Ib**2*omega**2*rho**2 + 4*A*E*Ib*rho))*omega)/2
        alpha_3 = 1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega - np.sqrt(Ib**2*omega**2*rho**2 + 4*A*E*Ib*rho))*omega)/2
        alpha_4 = -1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega - np.sqrt(Ib**2*omega**2*rho**2 + 4*A*E*Ib*rho))*omega)/2

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
        L = self.L
        rho = self.rho
        A = self.A
        E = self.E
        Ib = self.Ib

        # determine wavenumber
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # extract loads
        q_z = q[0]
        q_phi = q[1]

        # TODO - check if it is correct
        el = np.array([-1j*Ib*E*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A + 1j*Ib*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3))*(alpha_1 - alpha_2))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A,
                       Ib*E*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*alpha_1)/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_phi/omega**2/rho/A - Ib*E*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + (alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3))*(alpha_1 - alpha_2))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_phi/omega**2/rho/A,
                       1j*Ib*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A + 1j*Ib*(alpha_1*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*alpha_4)*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_z/omega**2/rho/A,
                       -Ib*E*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_phi/omega**2/rho/A - Ib*(alpha_1*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*alpha_4)*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*q_phi/omega**2/rho/A])
        return el  

    def LocalElementDisplacements(self, u_nodes_global, omega, num_points):
        """
        Calculate local tranlsational w(s) and rotational displacement phi(s).

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

        # determine nodal displacement in local axis
        u_node_local = self.R @ u_nodes_global

        # extract only the needed displacements
        u_node_local = u_node_local(self.dofs)

        # calculate coeficients
        C = self.Coefficients(u_node_local, omega)

        # get displacement
        w = self.displacement(x, omega, C)

        # get rotations
        phi = self.rotation(x, omega, C)

        return [w, phi]

    def Coefficients(self, u_node_local, omega):
        """
        Calculate the coefficients of the general solution, in this case 4.

        Input:
            u_node_local: local degrees of freedom
            omega: array. Range of frequencies of analysis
        Output:
            C: array (4), coefficients of general solution (C1, C2, C3, C4)
        """
        # Read all the variables
        L = self.L
        rho = self.who
        A = self.A
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # get distributed load value
        q = self.q[1]
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

        # TODO - check correctness
        u_load = np.array([q/(omega**2*rho*A), 0, q/(omega**2*rho*A), 0])
        C = np.linalg.inv(A_mat) @ (u_node_local + u_load)
        return C

    def displacement(self, x, omega, C=None, u_node_local=None):
        """
        Get the transverse displacments of the Rayleigh beam.

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
        if C is None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        w = (C[0]*np.exp(-1j*alpha_1*x) +
             C[1]*np.exp(-1j*alpha_2*x) +
             C[2]*np.exp(-1j*alpha_3*x) +
             C[3]*np.exp(-1j*alpha_4*x))

        return w

    def rotation(self, x, omega, C=None, u_node_local=None):
        """
        Get the rotations of the Rayleigh beam.

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
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # check if C is input, otherwise calculate
        if C is None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        phi = (-1j*C[0]*alpha_1*np.exp(-1j*alpha_1*x) -
               1j*C[1]*alpha_2*np.exp(-1j*alpha_2*x) -
               1j*C[2]*alpha_3*np.exp(-1j*alpha_3*x) -
               1j*C[3]*alpha_4*np.exp(-1j*alpha_4*x))

        return -phi
