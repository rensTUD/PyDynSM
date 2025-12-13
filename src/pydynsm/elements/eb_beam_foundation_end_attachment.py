# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:27:06 2024

@author: rensv
"""
# import dependencies
import numpy as np
from .structuralelement import StructuralElement, ElementFactory

# %% class definition


@ElementFactory.ElementType('EulerBernoulli Beam foundation attachments')
class EulerBernoulliBeamFoundationEndAttachment(StructuralElement):
    """Class for Euler-Bernoulli beam element over viscoelastic foundation with end attachments."""


    def __init__(self, section, rho, E, L, kd=0, cd=0, Pm1=0, Pm2=0, J1=0, J2=0, K1=0, K2=0, Cd1=0, Cd2=0, Kr1=0, Kr2=0, Cr1=0, Cr2=0, ksi=None ):
        """
        Initialise EB beam over viscoelastic foundation.

        Parameters
        ----------
        section : Section
            Section object providing geometric properties (A, I_y, W_y)
        rho : float
            Density of element [kg/m^3]
        E : float
            Young's modulus [Pa]
        L : float
            Length of the element [m]
        kd : float, optional
            Stiffness of viscoelastic foundation [N/m^2], default: 0
        cd : float, optional
            Damping of viscoelastic foundation [N/s/m^2], default: 0
        ksi : float, optional
            Material damping [-], default: 0
        Pm1, Pm2 : float, optional
            Point Mass at the left (1) and right (2) end of the beam [kg], default: 0
        J1, J2 : float, optional
            Point Mass moment of inertia at the left (1) and right (2) end of the beam [kg*m^2], default: 0
        K1, K2 : float, optional
            Point spring stiffness at the left (1) and right (2) end of the beam [N/m], default: 0
        Cd1, Cd2 : float, optional
            Point dashpot damping at the left (1) and right (2) end of the beam [N*s/m], default: 0
        Kr1, Kr2 : float, optional
            Point rotational spring stiffness at the left (1) and right (2) end of the beam [N*m/rad], default: 0
        Cr1, Cr2 : float, optional
            Point rotational dashpot damping at the left (1) and right (2) end of the beam [N*m*s/rad], default: 0
        """
        # define the of dofs the EulerBernoulli beam and initialise
        dofs = ['z','phi_y']
        super().__init__(dofs)

        # Extract geometric properties from section
        # For 2D x-z plane bending: Ib maps to I_y, Wb maps to W_y
        self.A = section.A
        self.Ib = section.I_y
        self.Wb = section.W_y

        # Initialise local beam element with necessary parameters
        self.rho = rho
        self.L = L
        self.kd = kd
        self.cd = cd
        # Assign ksi if given,  otherwise assign a default value
        self.ksi = ksi if ksi is not None else 0
        self.E = E*(1+2j*self.ksi)
        self.q = np.zeros(len(dofs))
        self.Pm1 = Pm1
        self.Pm2 = Pm2
        self.J1 = J1
        self.J2 = J2
        self.K1 = K1
        self.K2 = K2
        self.Cd1 = Cd1
        self.Cd2 = Cd2
        self.Kr1 = Kr1
        self.Kr2 = Kr2
        self.Cr1 = Cr1
        self.Cr2 = Cr2 
        

    def LocalStiffness(self, omega):
        """
        Determine the stiffness of the EB-beam over viscoelastic foundation.

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
        I_b = self.Ib
        L = self.L
        Pm1 = self.Pm1
        Pm2 = self.Pm2
        J1 = self.J1
        J2 = self.J2
        K1 = self.K1
        K2 = self.K2
        Cd1 = self.Cd1
        Cd2 = self.Cd2
        Kr1 = self.Kr1
        Kr2 = self.Kr2
        Cr1 = self.Cr1
        Cr2 = self.Cr2


        # Obtain the wavenumbers
        alpha1, alpha2, alpha3, alpha4 = self.ElementWaveNumbers(omega)

        # Initialise K_local
        K_local = np.empty((4, 4), complex)
        # Copy-paste the matrix directly from Maple (also called K_dyn)
        
        K_local = np.array([[((-alpha3 + alpha4) * (alpha1 - alpha2) * (1j * alpha4 ** 2 * alpha3 * I_b * E + 1j * alpha3 ** 2 * alpha4 * I_b * E + -1j * Cd1 * omega + Pm1 * omega ** 2 - K1) * np.exp(-1j * L * (alpha1 + alpha2)) - (1j * alpha4 ** 2 * alpha2 * I_b * E + 1j * alpha2 ** 2 * alpha4 * I_b * E + -1j * Cd1 * omega + Pm1 * omega ** 2 - K1) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (1j * alpha2 ** 2 * alpha3 * I_b * E + 1j * alpha3 ** 2 * alpha2 * I_b * E + -1j * Cd1 * omega + Pm1 * omega ** 2 - K1) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * (1j * alpha4 ** 2 * alpha1 * I_b * E + 1j * alpha1 ** 2 * alpha4 * I_b * E + -1j * Cd1 * omega + Pm1 * omega ** 2 - K1) * np.exp(-1j * L * (alpha2 + alpha3)) - (1j * alpha1 ** 2 * alpha3 * I_b * E + 1j * alpha3 ** 2 * alpha1 * I_b * E + -1j * Cd1 * omega + Pm1 * omega ** 2 - K1) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + (-alpha3 + alpha4) * (1j * alpha1 ** 2 * alpha2 * I_b * E + 1j * alpha2 ** 2 * alpha1 * I_b * E + -1j * Cd1 * omega + Pm1 * omega ** 2 - K1) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha3 + alpha4))) / ((alpha1 - alpha2) * (-alpha3 + alpha4) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha1 - alpha2) * (-alpha3 + alpha4)),-E * I_b * (((alpha3 ** 3 - alpha4 ** 3) * alpha1 + (-alpha3 ** 3 + alpha4 ** 3) * alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) + ((-alpha2 ** 3 + alpha4 ** 3) * alpha1 + alpha2 ** 3 * alpha3 - alpha4 ** 3 * alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + ((alpha2 ** 3 - alpha3 ** 3) * alpha1 - alpha2 ** 3 * alpha4 + alpha3 ** 3 * alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 ** 3 - alpha4 ** 3) * np.exp(-1j * L * (alpha2 + alpha3)) + ((-alpha2 + alpha4) * alpha1 ** 3 + alpha3 ** 3 * (alpha2 - alpha4)) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha1 - alpha2) * (alpha1 ** 2 + alpha1 * alpha2 + alpha2 ** 2) * (alpha3 - alpha4)) / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),-1j * (alpha1 * (alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * (alpha2 + alpha3 + alpha4) * np.exp(-1j * alpha1 * L) - alpha2 * (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * (alpha1 + alpha3 + alpha4) * np.exp(-1j * alpha2 * L) + (alpha1 - alpha2) * (alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 + alpha2 + alpha4) * np.exp(-1j * alpha3 * L) - np.exp(-1j * alpha4 * L) * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 + alpha2 + alpha3))) * E * I_b / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),E * I_b * ((alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * (alpha2 + alpha3 + alpha4) * np.exp(-1j * alpha1 * L) - (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * (alpha1 + alpha3 + alpha4) * np.exp(-1j * alpha2 * L) + ((alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 + alpha2 + alpha4) * np.exp(-1j * alpha3 * L) - np.exp(-1j * alpha4 * L) * (alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 + alpha2 + alpha3)) * (alpha1 - alpha2)) / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2))],[-E * I_b * (alpha3 * alpha4 * (alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - alpha2 * alpha4 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + alpha2 * alpha3 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + alpha1 * (alpha4 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - alpha3 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * alpha2 * (alpha3 - alpha4) * (alpha1 - alpha2))) / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),(-(1j * alpha4 * I_b * E + 1j * alpha3 * I_b * E + 1j * omega * Cr1 - omega ** 2 * J1 + Kr1) * (-alpha3 + alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) + (1j * alpha4 * I_b * E + 1j * alpha2 * I_b * E + 1j * omega * Cr1 - omega ** 2 * J1 + Kr1) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) - (1j * alpha2 * I_b * E + 1j * alpha3 * I_b * E + 1j * omega * Cr1 - omega ** 2 * J1 + Kr1) * (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) - (1j * alpha4 * I_b * E + 1j * alpha1 * I_b * E + 1j * omega * Cr1 - omega ** 2 * J1 + Kr1) * (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) + (1j * alpha1 * I_b * E + 1j * alpha3 * I_b * E + 1j * omega * Cr1 - omega ** 2 * J1 + Kr1) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) - (1j * alpha1 * I_b * E + 1j * alpha2 * I_b * E + 1j * omega * Cr1 - omega ** 2 * J1 + Kr1) * (-alpha3 + alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha3 + alpha4))) / ((alpha1 - alpha2) * (-alpha3 + alpha4) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha1 - alpha2) * (-alpha3 + alpha4)),(alpha1 * (alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * np.exp(-1j * alpha1 * L) - alpha2 * (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * alpha2 * L) + (alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * np.exp(-1j * alpha3 * L) - np.exp(-1j * alpha4 * L) * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3)) * (alpha1 - alpha2)) * E * I_b / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),1j * E * I_b * ((alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * np.exp(-1j * alpha1 * L) - (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * alpha2 * L) + ((alpha2 - alpha4) * (alpha1 - alpha4) * np.exp(-1j * alpha3 * L) - np.exp(-1j * alpha4 * L) * (alpha2 - alpha3) * (alpha1 - alpha3)) * (alpha1 - alpha2)) / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2))],[-1j * E * (alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 - alpha2) * (alpha1 + alpha2 + alpha3) * np.exp(-1j * L * (alpha1 + alpha2 + alpha3)) - alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha2) * (alpha1 + alpha2 + alpha4) * np.exp(-1j * L * (alpha1 + alpha2 + alpha4)) + (alpha3 - alpha4) * (alpha2 * (alpha1 - alpha4) * (alpha1 - alpha3) * (alpha1 + alpha3 + alpha4) * np.exp(-1j * L * (alpha1 + alpha3 + alpha4)) - np.exp(-1j * L * (alpha2 + alpha3 + alpha4)) * alpha1 * (alpha2 - alpha4) * (alpha2 - alpha3) * (alpha2 + alpha3 + alpha4))) * I_b / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),E * ((alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 - alpha2) * (alpha1 + alpha2 + alpha3) * np.exp(-1j * L * (alpha1 + alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha2) * (alpha1 + alpha2 + alpha4) * np.exp(-1j * L * (alpha1 + alpha2 + alpha4)) + ((alpha1 - alpha4) * (alpha1 - alpha3) * (alpha1 + alpha3 + alpha4) * np.exp(-1j * L * (alpha1 + alpha3 + alpha4)) - np.exp(-1j * L * (alpha2 + alpha3 + alpha4)) * (alpha2 - alpha4) * (alpha2 - alpha3) * (alpha2 + alpha3 + alpha4)) * (alpha3 - alpha4)) * I_b / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),(-(-alpha3 + alpha4) * (1j * alpha1 ** 2 * alpha2 * I_b * E + 1j * alpha2 ** 2 * alpha1 * I_b * E + 1j * omega * Cd2 - omega ** 2 * Pm2 + K2) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) + (1j * alpha1 ** 2 * alpha3 * I_b * E + 1j * alpha3 ** 2 * alpha1 * I_b * E + 1j * omega * Cd2 - omega ** 2 * Pm2 + K2) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) - (alpha2 - alpha3) * (-alpha1 + alpha4) * (1j * alpha4 ** 2 * alpha1 * I_b * E + 1j * alpha1 ** 2 * alpha4 * I_b * E + 1j * omega * Cd2 - omega ** 2 * Pm2 + K2) * np.exp(-1j * L * (alpha1 + alpha4)) - (1j * alpha2 ** 2 * alpha3 * I_b * E + 1j * alpha3 ** 2 * alpha2 * I_b * E + 1j * omega * Cd2 - omega ** 2 * Pm2 + K2) * (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) + (1j * alpha4 ** 2 * alpha2 * I_b * E + 1j * alpha2 ** 2 * alpha4 * I_b * E + 1j * omega * Cd2 - omega ** 2 * Pm2 + K2) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) - (-alpha3 + alpha4) * (1j * alpha4 ** 2 * alpha3 * I_b * E + 1j * alpha3 ** 2 * alpha4 * I_b * E + 1j * omega * Cd2 - omega ** 2 * Pm2 + K2) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha3 + alpha4))) / ((alpha1 - alpha2) * (-alpha3 + alpha4) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha1 - alpha2) * (-alpha3 + alpha4)),E * I_b * ((alpha3 - alpha4) * (alpha1 ** 3 - alpha2 ** 3) * np.exp(-1j * L * (alpha1 + alpha2)) + ((-alpha2 + alpha4) * alpha1 ** 3 + alpha3 ** 3 * (alpha2 - alpha4)) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 ** 3 - alpha4 ** 3) * np.exp(-1j * L * (alpha1 + alpha4)) + ((alpha2 ** 3 - alpha3 ** 3) * alpha1 - alpha2 ** 3 * alpha4 + alpha3 ** 3 * alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) + ((-alpha2 ** 3 + alpha4 ** 3) * alpha1 + alpha2 ** 3 * alpha3 - alpha4 ** 3 * alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha3 ** 2 + alpha3 * alpha4 + alpha4 ** 2) * (alpha1 - alpha2)) / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2))],[E * I_b * (alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2 + alpha3)) - alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2 + alpha4)) + (alpha2 * (alpha1 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3 + alpha4)) - np.exp(-1j * L * (alpha2 + alpha3 + alpha4)) * alpha1 * (alpha2 - alpha4) * (alpha2 - alpha3)) * (alpha3 - alpha4)) / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),1j * E * ((alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2 + alpha4)) + ((alpha1 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3 + alpha4)) - np.exp(-1j * L * (alpha2 + alpha3 + alpha4)) * (alpha2 - alpha4) * (alpha2 - alpha3)) * (alpha3 - alpha4)) * I_b / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),(alpha1 * alpha2 * (alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - alpha1 * alpha3 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + alpha1 * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + alpha2 * alpha3 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) - np.exp(-1j * L * (alpha3 + alpha4)) * alpha3 * (alpha3 - alpha4) * (alpha1 - alpha2)) * alpha4) * E * I_b / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)),((-alpha3 + alpha4) * (1j * alpha1 * I_b * E + 1j * alpha2 * I_b * E + -1j * omega * Cr2 + omega ** 2 * J2 - Kr2) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (1j * alpha1 * I_b * E + 1j * alpha3 * I_b * E + -1j * omega * Cr2 + omega ** 2 * J2 - Kr2) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (1j * alpha4 * I_b * E + 1j * alpha1 * I_b * E + -1j * omega * Cr2 + omega ** 2 * J2 - Kr2) * (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (1j * alpha2 * I_b * E + 1j * alpha3 * I_b * E + -1j * omega * Cr2 + omega ** 2 * J2 - Kr2) * (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (1j * alpha4 * I_b * E + 1j * alpha2 * I_b * E + -1j * omega * Cr2 + omega ** 2 * J2 - Kr2) * (-alpha2 + alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + (-alpha3 + alpha4) * (1j * alpha4 * I_b * E + 1j * alpha3 * I_b * E + -1j * omega * Cr2 + omega ** 2 * J2 - Kr2) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha3 + alpha4))) / ((alpha1 - alpha2) * (-alpha3 + alpha4) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (-alpha1 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha1 - alpha3) * (-alpha2 + alpha4) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha1 - alpha2) * (-alpha3 + alpha4))]])

        return K_local

    def ElementWaveNumbers(self, omega):
        """
        Determine the wavenumbers of EB beam over viscoelastic foundation.

        Input:
            omega: array. Range of frequencies of analysis.
        Output:
            alpha_1, alpha_2, alpha_3, alpha_4: values. Wavenumbers
        """
        rho = self.rho
        A = self.A
        E = self.E
        Ib = self.Ib
        kd = self.kd
        cd = self.cd
        # Copy-paste from Maple document of Wavenumbers
        # 4th-order derivative in x, so 4 wavenumbers
        alpha_1 = 1/E/Ib*((A*omega**2*rho + -1j * cd * omega - kd)*E**3*Ib**3)**(1/4)
        alpha_2 = 1j/E/Ib*((A*omega**2*rho + -1j * cd * omega - kd)*E**3*Ib**3)**(1/4)
        alpha_3 = -1/E/Ib*((A*omega**2*rho + -1j * cd * omega - kd)*E**3*Ib**3)**(1/4)
        alpha_4 = -1j/E/Ib*((A*omega**2*rho + -1j * cd * omega - kd)*E**3*Ib**3)**(1/4)

        return alpha_1, alpha_2, alpha_3, alpha_4

    def LocalDistributedLoad(self, q, omega):
        """
        Add a distributed load to the local element.

        Input:
            q: array. Distributed load. With a shape such as:
            q = [q_z, q_phi]
            omega: array. Range of frequencies of analyis
        Output:
            el: array. Force vector.
        """
        # assign load to itself to keep track
        self.q = q

        # assign local variables for ease of coding
        rho = self.rho
        A = self.A
        E = self.E
        I_b = self.Ib
        L = self.L
        kd = self.kd
        cd = self.cd
        Pm1 = self.Pm1
        Pm2 = self.Pm2
        J1 = self.J1
        J2 = self.J2
        K1 = self.K1
        K2 = self.K2
        Cd1 = self.Cd1
        Cd2 = self.Cd2
        Kr1 = self.Kr1
        Kr2 = self.Kr2
        Cr1 = self.Cr1
        Cr2 = self.Cr2

        # determine wavenumber
        alpha1, alpha2, alpha3, alpha4 = self.ElementWaveNumbers(omega)

        # extract loads
        q_z = q[0]
        q_phi = q[1]

        # TODO - check for correctness
        el = np.array([0.1e1 / (-1j * omega ** 2 * rho * A + 1j * kd - cd * omega) * q_z / (-(alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) + (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) - (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) - (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) + (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) - np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)) * ((alpha1 - alpha2) * (alpha3 - alpha4) * np.exp(-1j * L * (alpha1 + alpha2)) * (alpha3 ** 2 * alpha4 * I_b * E + alpha4 ** 2 * alpha3 * I_b * E + -1j * Pm1 * omega ** 2 + 1j * K1 - Cd1 * omega) + (-alpha1 + alpha3) * (alpha2 - alpha4) * np.exp(-1j * L * (alpha1 + alpha3)) * (alpha2 ** 2 * alpha4 * I_b * E + alpha4 ** 2 * alpha2 * I_b * E + -1j * Pm1 * omega ** 2 + 1j * K1 - Cd1 * omega) + (alpha1 - alpha4) * (alpha2 - alpha3) * np.exp(-1j * L * (alpha1 + alpha4)) * (alpha2 ** 2 * alpha3 * I_b * E + alpha3 ** 2 * alpha2 * I_b * E + -1j * Pm1 * omega ** 2 + 1j * K1 - Cd1 * omega) + (alpha1 - alpha4) * (alpha2 - alpha3) * np.exp(-1j * L * (alpha2 + alpha3)) * (alpha1 ** 2 * alpha4 * I_b * E + alpha4 ** 2 * alpha1 * I_b * E + -1j * Pm1 * omega ** 2 + 1j * K1 - Cd1 * omega) + (-alpha1 + alpha3) * (alpha2 - alpha4) * np.exp(-1j * L * (alpha2 + alpha4)) * (alpha1 ** 2 * alpha3 * I_b * E + alpha3 ** 2 * alpha1 * I_b * E + -1j * Pm1 * omega ** 2 + 1j * K1 - Cd1 * omega) + (alpha1 - alpha2) * (alpha3 - alpha4) * np.exp(-1j * L * (alpha3 + alpha4)) * (alpha1 ** 2 * alpha2 * I_b * E + alpha2 ** 2 * alpha1 * I_b * E + -1j * Pm1 * omega ** 2 + 1j * K1 - Cd1 * omega) + (-alpha1 * (alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * (alpha2 + alpha3 + alpha4) * np.exp(-1j * alpha1 * L) + alpha2 * (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * (alpha1 + alpha3 + alpha4) * np.exp(-1j * alpha2 * L) - (alpha1 - alpha2) * (alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 + alpha2 + alpha4) * np.exp(-1j * alpha3 * L) - np.exp(-1j * alpha4 * L) * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 + alpha2 + alpha3))) * I_b * E),I_b * E * (alpha3 * alpha4 * (alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - alpha2 * alpha4 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + alpha2 * alpha3 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + alpha1 * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - alpha1 * alpha3 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + alpha1 * alpha2 * (alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha3 + alpha4)) - alpha1 * (alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * np.exp(-1j * alpha1 * L) + alpha2 * (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * alpha2 * L) - (alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * np.exp(-1j * alpha3 * L) - np.exp(-1j * alpha4 * L) * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3)) * (alpha1 - alpha2)) * q_z / (-(alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) + (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) - (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) - (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) + (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) - np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)) / (-omega ** 2 * rho * A + kd + 1j * cd * omega),0.1e1 / (1j * omega ** 2 * rho * A + -1j * kd + cd * omega) * q_z / (-(alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) + (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) - (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) - (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) + (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) - np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)) * (E * I_b * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 - alpha2) * (alpha1 + alpha2 + alpha3) * np.exp(-1j * L * (alpha1 + alpha2 + alpha3)) - E * I_b * alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha2) * (alpha1 + alpha2 + alpha4) * np.exp(-1j * L * (alpha1 + alpha2 + alpha4)) + alpha2 * E * I_b * (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * (alpha1 + alpha3 + alpha4) * np.exp(-1j * L * (alpha1 + alpha3 + alpha4)) - alpha1 * E * I_b * (alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * (alpha2 + alpha3 + alpha4) * np.exp(-1j * L * (alpha2 + alpha3 + alpha4)) + (alpha1 ** 2 * alpha2 * I_b * E + alpha2 ** 2 * alpha1 * I_b * E + 1j * Pm2 * omega ** 2 + -1j * K2 + omega * Cd2) * (alpha1 - alpha2) * (alpha3 - alpha4) * np.exp(-1j * L * (alpha1 + alpha2)) + (-alpha1 + alpha3) * (alpha2 - alpha4) * (alpha1 ** 2 * alpha3 * I_b * E + alpha3 ** 2 * alpha1 * I_b * E + 1j * Pm2 * omega ** 2 + -1j * K2 + omega * Cd2) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha1 - alpha4) * (alpha1 ** 2 * alpha4 * I_b * E + alpha4 ** 2 * alpha1 * I_b * E + 1j * Pm2 * omega ** 2 + -1j * K2 + omega * Cd2) * (alpha2 - alpha3) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha1 - alpha4) * (alpha2 - alpha3) * (alpha2 ** 2 * alpha3 * I_b * E + alpha3 ** 2 * alpha2 * I_b * E + 1j * Pm2 * omega ** 2 + -1j * K2 + omega * Cd2) * np.exp(-1j * L * (alpha2 + alpha3)) + (-alpha1 + alpha3) * (alpha2 ** 2 * alpha4 * I_b * E + alpha4 ** 2 * alpha2 * I_b * E + 1j * Pm2 * omega ** 2 + -1j * K2 + omega * Cd2) * (alpha2 - alpha4) * np.exp(-1j * L * (alpha2 + alpha4)) + (alpha3 ** 2 * alpha4 * I_b * E + alpha4 ** 2 * alpha3 * I_b * E + 1j * Pm2 * omega ** 2 + -1j * K2 + omega * Cd2) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4)),0.1e1 / ((alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) + np.exp(-1j * L * (alpha3 + alpha4)) * (alpha3 - alpha4) * (alpha1 - alpha2)) * I_b * (alpha4 * (alpha2 - alpha3) * (alpha1 - alpha3) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2 + alpha3)) - alpha3 * (alpha2 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2 + alpha4)) + alpha2 * (alpha3 - alpha4) * (alpha1 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3 + alpha4)) - alpha1 * (alpha3 - alpha4) * (alpha2 - alpha4) * (alpha2 - alpha3) * np.exp(-1j * L * (alpha2 + alpha3 + alpha4)) + alpha1 * alpha2 * (alpha3 - alpha4) * (alpha1 - alpha2) * np.exp(-1j * L * (alpha1 + alpha2)) - alpha1 * alpha3 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha1 + alpha3)) + alpha1 * alpha4 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha1 + alpha4)) + alpha2 * alpha3 * (alpha2 - alpha3) * (alpha1 - alpha4) * np.exp(-1j * L * (alpha2 + alpha3)) - (alpha2 * (alpha2 - alpha4) * (alpha1 - alpha3) * np.exp(-1j * L * (alpha2 + alpha4)) - np.exp(-1j * L * (alpha3 + alpha4)) * alpha3 * (alpha3 - alpha4) * (alpha1 - alpha2)) * alpha4) * E * q_z / (-omega ** 2 * rho * A + kd + 1j * cd * omega)])



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
        rho = self.rho
        A = self.A
        L = self.L
        kd = self.kd
        cd = self.cd
        alpha_1, alpha_2, alpha_3, alpha_4 = self.ElementWaveNumbers(omega)

        # get distributed load value
        q_z = self.q[0]
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
        u_load = np.array([-q_z/(A*omega**2*rho - kd - 1j*cd*omega),
                           0,
                           -q_z/(A*omega**2*rho - kd - 1j*cd*omega),
                           0])
        C = np.linalg.inv(A_mat) @ (u_node_local - u_load)
        # + np.array([1/(E*Ib*beta_b**4),0,1/(E*Ib*beta_b**4),0]) * q)

        return C

    def displacement(self, x, omega, C=None, u_node_local=None):
        """
        Get the transverse displacments of the EB beam over viscoel foundation.

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

        # Check if C is input, otherwise, calculate it
        if C is None:
            C = self.Coefficients(u_node_local, omega)

        # displacements
        w = (C[0]*np.exp(-1j*alpha_1*x) +
             C[1]*np.exp(-1j*alpha_2*x) +
             C[2]*np.exp(-1j*alpha_3*x) +
             C[3]*np.exp(-1j*alpha_4*x)) - self.q[0]/(self.rho*self.A*omega**2 - self.kd - 1j*self.cd*omega)

        return w

    def rotation(self, x, omega, C=None, u_node_local=None):
        """
        Get the rotations of the EB beam over viscoelastic foundation.

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

        # check if C is input
        if C is None:
            C = self.Coefficients(u_node_local, omega)

        # Rotations
        phi = (-1j*C[0]*alpha_1*np.exp(-1j*alpha_1*x) -
               1j*C[1]*alpha_2*np.exp(-1j*alpha_2*x) -
               1j*C[2]*alpha_3*np.exp(-1j*alpha_3*x) -
               1j*C[3]*alpha_4*np.exp(-1j*alpha_4*x))

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
        
        dudx2 = ((-1j*alpha_1)**2*C[0]*np.exp(-1j*alpha_1*x) +
                (-1j*alpha_2)**2*C[1]*np.exp(-1j*alpha_2*x) +
                (-1j*alpha_3)**2*C[2]*np.exp(-1j*alpha_3*x) +
                (-1j*alpha_4)**2*C[3]*np.exp(-1j*alpha_4*x))
        
        mmt = -E*I*dudx2
        return mmt    

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