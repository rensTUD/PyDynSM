# -*- coding: utf-8 -*-
"""
Created on [Date]

@author: rensv
"""

import numpy as np
from .section import Section, SectionFactory


@SectionFactory.SectionType('rectangle')
class Rectangle(Section):
    """
    Rectangular solid cross-section.
    
    Parameters
    ----------
    width : float
        Width of the rectangle in y-direction [m]
    height : float
        Height of the rectangle in z-direction [m]
    
    Coordinate system:
        - x: horizontal/axial direction
        - z: vertical direction (height)
        - y: in-plane horizontal direction (width)
    """
    
    section_name = 'rectangle'
    
    def __init__(self, width, height):
        """
        Initialize rectangular section.
        
        Parameters
        ----------
        width : float
            Width in y-direction [m]
        height : float
            Height in z-direction [m]
        """
        if width <= 0 or height <= 0:
            raise ValueError("Width and height must be positive.")
        
        self.width = width
        self.height = height
        
        super().__init__()
    
    def compute_properties(self):
        """
        Compute geometric properties for rectangular section.
        
        Formulas:
            A = width * height
            I_y = width * height^3 / 12  (bending in x-z plane)
            I_z = height * width^3 / 12  (bending in x-y plane)
            W_y = width * height^2 / 6  (section modulus about y-axis)
            W_z = height * width^2 / 6  (section modulus about z-axis)
        """
        width = self.width
        height = self.height
        
        # Cross-sectional area
        self._A = width * height
        
        # Second moments of area
        # I_y: bending in x-z plane (about y-axis)
        self._I_y = width * height**3 / 12
        
        # I_z: bending in x-y plane (about z-axis)
        self._I_z = height * width**3 / 12
        
        # Section moduli
        # W_y = I_y / c_y, where c_y = height/2 (distance to extreme fiber)
        self._W_y = self._I_y / (height / 2)
        
        # W_z = I_z / c_z, where c_z = width/2 (distance to extreme fiber)
        self._W_z = self._I_z / (width / 2)
