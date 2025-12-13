# -*- coding: utf-8 -*-
"""
Created on [Date]

@author: rensv
"""

import numpy as np
from .section import Section, SectionFactory


@SectionFactory.SectionType('hollow_circle')
class HollowCircle(Section):
    """
    Circular hollow cross-section (tube).
    
    Parameters
    ----------
    outer_diameter : float
        Outer diameter of the circle [m]
    inner_diameter : float
        Inner diameter of the circle [m]
    
    Coordinate system:
        - x: horizontal/axial direction
        - z: vertical direction
        - y: in-plane horizontal direction
    """
    
    section_name = 'hollow_circle'
    
    def __init__(self, outer_diameter, inner_diameter):
        """
        Initialize hollow circular section.
        
        Parameters
        ----------
        outer_diameter : float
            Outer diameter [m]
        inner_diameter : float
            Inner diameter [m]
        """
        if outer_diameter <= 0 or inner_diameter <= 0:
            raise ValueError("Outer and inner diameters must be positive.")
        if inner_diameter >= outer_diameter:
            raise ValueError("Inner diameter must be less than outer diameter.")
        
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter
        
        super().__init__()
    
    def compute_properties(self):
        """
        Compute geometric properties for hollow circular section.
        
        Formulas:
            A = π * (outer_diameter^2 - inner_diameter^2) / 4
            I_y = I_z = π * (outer_diameter^4 - inner_diameter^4) / 64
            W_y = W_z = I / (outer_diameter/2) = π * (outer_diameter^4 - inner_diameter^4) / (32 * outer_diameter)
        """
        d_o = self.outer_diameter
        d_i = self.inner_diameter
        r_o = d_o / 2  # outer radius
        r_i = d_i / 2  # inner radius
        
        # Cross-sectional area
        self._A = np.pi * (r_o**2 - r_i**2)
        
        # Second moments of area (same for any axis due to circular symmetry)
        self._I_y = np.pi * (d_o**4 - d_i**4) / 64
        self._I_z = self._I_y  # Same due to symmetry
        
        # Section moduli
        # W = I / c, where c = r_o = d_o/2 (distance to extreme fiber)
        self._W_y = self._I_y / r_o
        self._W_z = self._W_y  # Same due to symmetry
