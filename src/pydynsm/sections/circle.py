# -*- coding: utf-8 -*-
"""
Created on [Date]

@author: rensv
"""

import numpy as np
from .section import Section, SectionFactory


@SectionFactory.SectionType('circle')
class Circle(Section):
    """
    Circular solid cross-section.
    
    Parameters
    ----------
    diameter : float
        Diameter of the circle [m]
    
    Coordinate system:
        - x: horizontal/axial direction
        - z: vertical direction
        - y: in-plane horizontal direction
    """
    
    section_name = 'circle'
    
    def __init__(self, diameter):
        """
        Initialize circular section.
        
        Parameters
        ----------
        diameter : float
            Diameter of the circle [m]
        """
        if diameter <= 0:
            raise ValueError("Diameter must be positive.")
        
        self.diameter = diameter
        
        super().__init__()
    
    def compute_properties(self):
        """
        Compute geometric properties for circular section.
        
        Formulas:
            A = π * (diameter/2)^2
            I_y = I_z = π * diameter^4 / 64  (same for any axis due to symmetry)
            W_y = W_z = π * diameter^3 / 32  (section modulus)
        """
        d = self.diameter
        r = d / 2  # radius
        
        # Cross-sectional area
        self._A = np.pi * r**2
        
        # Second moments of area (same for any axis due to circular symmetry)
        self._I_y = np.pi * d**4 / 64
        self._I_z = self._I_y  # Same due to symmetry
        
        # Section moduli
        # W = I / c, where c = r = d/2 (distance to extreme fiber)
        self._W_y = self._I_y / r
        self._W_z = self._W_y  # Same due to symmetry
