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
class Rayleigh_Beam(StructuralElement):

    element_name = 'Rayleigh Beam'

    def __init__(self, rho, A, E, Ib, L, ksi=None):

        # define what dofs the rayleigh beam contributes to and initialise
        # the dofs are w and phi, so 1 and 2
        dofs = [1, 2]
        super().__init__(dofs)

        # Initialise local beam element with necessary parameters
        # assisgn ksi if given otherwise assign a default value (think of changing this)
        self.ksi = ksi if ksi is not None else 0.01
        self.rho = rho
        self.A = A
        self.E = E
        self.Ib = Ib
        self.L = L

    def LocalStiffness(self, omega):
        '''
        Determines the stiffness of the beam.
        As its 2D, stiffness will be 4x4 matrix as:
        [V_left, M_left, V_right, M_right] = K.[w_left, phi_left, w_right, phi_right]

        where: ASK THIS
            K = [K_V_ll, K_M_ll, K_V_lr, K_M_lr;
                 K_V_rl, K_M_rl, K_V_rr, K_M_rr]
        '''

        # assign local variables for ease of coding
        E = self.E
        Ib = self.Ib
        rho = self.rho
        L = self.L
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # Initialize K_local
        K_local = np.empty((4, 4), complex)
        # we can also create the stiffness matrix like this,
        # probably easier to get direct output from Maple but maybe less clear
        # and less easy to spot mistakes (if present)
        K_local = np.array([[-1j*Ib*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), Ib*((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_3 ** 2 + E*alpha_3*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_2 ** 2 + E*alpha_2*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_2 ** 2 + E*alpha_2*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_1 ** 2 + E*alpha_1*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_1 ** 2 + E*alpha_1*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_1 ** 2 + E*alpha_1*alpha_2 + E*alpha_2 ** 2 + omega ** 2*rho))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), 1j*Ib*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -Ib*((alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)*np.exp(-1j*alpha_1*L) - (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 + alpha_2 + alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))], [Ib*(alpha_3*alpha_4*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_2*alpha_4*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_1*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), 1j*((alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4))*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -Ib*E*(alpha_1*(alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + (alpha_1 - alpha_2)*(alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -1j*Ib*((alpha_3 - alpha_4)*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*alpha_1*L) - (alpha_3 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*alpha_2*L) + ((alpha_2 - alpha_4)*(alpha_1 - alpha_4)*np.exp(-1j*alpha_3*L) - np.exp(-1j*alpha_4*L)*(alpha_2 - alpha_3)*(alpha_1 - alpha_3))*(alpha_1 - alpha_2))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))], [1j*Ib*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*(alpha_1 + alpha_2 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*(alpha_1 + alpha_3 + alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)*(alpha_2 + alpha_3 + alpha_4)))*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), 1j*(alpha_1*alpha_2*(alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*alpha_4)*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -Ib*((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_1 ** 2 + E*alpha_1*alpha_2 + E*alpha_2 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_1 ** 2 + E*alpha_1*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_1 ** 2 + E*alpha_1*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*(E*alpha_2 ** 2 + E*alpha_2*alpha_3 + E*alpha_3 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*(E*alpha_2 ** 2 + E*alpha_2*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*(E*alpha_3 ** 2 + E*alpha_3*alpha_4 + E*alpha_4 ** 2 + omega ** 2*rho))/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))], [-Ib*(alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*(alpha_2*(alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*alpha_1*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -1j*Ib*((alpha_2 - alpha_3)*(alpha_1 - alpha_3)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2 + alpha_4)) + (alpha_3 - alpha_4)*((alpha_1 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3 + alpha_4)) - np.exp(-1j*L*(alpha_2 + alpha_3 + alpha_4))*(alpha_2 - alpha_4)*(alpha_2 - alpha_3)))*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -Ib*E*(alpha_1*alpha_2*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - alpha_1*alpha_3*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + alpha_1*alpha_4*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + alpha_2*alpha_3*(alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2*(alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) - np.exp(-1j*L*(alpha_3 + alpha_4))*alpha_3*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))*alpha_4)/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2)), -1j*((alpha_1 - alpha_2)*(alpha_1 + alpha_2)*(alpha_3 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_1 - alpha_3)*(alpha_1 + alpha_3)*(alpha_2 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_1 - alpha_4)*(alpha_1 + alpha_4)*(alpha_2 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_2 + alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_2 + alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_3 + alpha_4)*(alpha_1 - alpha_2))*Ib*E/((alpha_3 - alpha_4)*(alpha_1 - alpha_2)*np.exp(-1j*L*(alpha_1 + alpha_2)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_1 + alpha_3)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_1 + alpha_4)) + (alpha_2 - alpha_3)*(alpha_1 - alpha_4)*np.exp(-1j*L*(alpha_2 + alpha_3)) - (alpha_2 - alpha_4)*(alpha_1 - alpha_3)*np.exp(-1j*L*(alpha_2 + alpha_4)) + np.exp(-1j*L*(alpha_3 + alpha_4))*(alpha_3 - alpha_4)*(alpha_1 - alpha_2))]])

        return K_local

    def ElementWaveNumbers(self, omega):
        '''
        Determines the wavenumbers of Rayleigh beam
        '''
        E = self.E
        Ib = self.Ib
        rho = self.rho
        A = self.A

        alpha_1 = 1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega + np.sqrt(Ib**2*rho**2*omega**2 + 4*E*Ib*A*rho))*omega)/2
        alpha_2 = -1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega + np.sqrt(Ib**2*omega**2*rho**2 + 4*A*E*Ib*rho))*omega)/2
        alpha_3 = 1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega - np.sqrt(Ib**2*omega**2*rho**2 + 4*A*E*Ib*rho))*omega)/2
        alpha_4 = -1/E/Ib*np.sqrt(2)*np.sqrt(E*Ib*(Ib*rho*omega - np.sqrt(Ib**2*omega**2*rho**2 + 4*A*E*Ib*rho))*omega)/2

        return alpha_1, alpha_2, alpha_3, alpha_4

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
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # extract loads
        q_z = q[0]
        q_phi = q[1]

        # TODO - change for model, ask how to do it
        # el = [ q_z*((np.cos(beta_b*L) - 1.0)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)),
        #                 -q_phi*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)),
        #                 q_z*((np.cos(beta_b*L) - 1)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)), 
        #                 q_phi*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)) 
        #                 ]

        return el  

    def LocalElementDisplacements(self, u_nodes_global, omega, num_points):
        '''
        function that calculates:
            coefficients C
            local displacements w(s)
            phi(s)
        '''

        # get local axis to evaluate on
        L = self.L
        x = np.linspace(0.0, L, num_points)

        # TODO - what is self.R
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
        '''
        Calculates the coefficients of the general solution
        '''
        # read all the variables
        L = self.L
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # TODO - should change this based on new derivations..(?)

        # get distributed load value
        q = self.q[1]
        # should be like:
        # q_b, q_m = self.

        # calculate the coefficients
        A_mat = np.array([[1, 1, 1, 1],
                          [1j*alpha_1, 1j*alpha_2, 1j*alpha_3, 1j*alpha_4],
                          [np.exp(-1j*alpha_1*L), np.exp(-1j*alpha_2*L),
                           np.exp(-1j*alpha_3*L), np.exp(-1j*alpha_4*L)],
                          [1j*np.exp(-1j*alpha_1*L)*alpha_1,
                           1j*np.exp(-1j*alpha_2*L)*alpha_2,
                           1j*np.exp(-1j*alpha_3*L)*alpha_3,
                           1j*np.exp(-1j*alpha_4*L)*alpha_4]])

        # TODO - include external forcing?
        C = A_mat @ (u_node_local)  # + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)

        return C

    def displacement(self, x, omega, C=None, u_node_local=None):
        '''
        Gets the transverse displacments of the Rayleigh beam

        if C is not given, then calculate it based on u_node_local.
        '''

        ksi = self.ksi
        EI = self.E * self.I*(1 + 2j * ksi)
        q = self.q[0]  # TODO not used? remove?
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # check if C is input        
        if C == None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        w = (C[0]*np.exp(-1j*alpha_1*x) +
             C[1]*np.exp(-1j*alpha_2*x) +
             C[2]*np.exp(-1j*alpha_3*x) +
             C[3]*np.exp(-1j*alpha_4*x))

        return w

    def rotation(self, x, omega, C=None, u_node_local=None):
        '''
        Gets the rotations of the Rayleigh beam
        '''


        ksi = self.ksi
        EI = self.E*self.I*(1 + 2j * ksi)
        q = self.q[1]  # TODO - not used - remove?
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

# %%% Register class with the ElementFactory

ElementFactory.RegisterElement(Rayleigh_Beam)
