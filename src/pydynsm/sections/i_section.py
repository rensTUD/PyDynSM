# -*- coding: utf-8 -*-
"""
Created on [Date]

@author: rensv
"""

import numpy as np
from .section import Section, SectionFactory


@SectionFactory.SectionType('i_section')
class ISection(Section):
    """
    Symmetric I or T cross-section.
    
    For symmetric I-section: top and bottom flanges have same width and thickness.
    For T-section: only one flange (treated as symmetric with zero bottom flange).
    
    Parameters
    ----------
    flange_width : float
        Width of each flange in y-direction [m]
    flange_thickness : float
        Thickness of each flange in z-direction [m]
    web_height : float
        Height of the web (between flanges) in z-direction [m]
    web_thickness : float
        Thickness of the web in y-direction [m]
    
    Coordinate system:
        - x: horizontal/axial direction
        - z: vertical direction (height)
        - y: in-plane horizontal direction (width)
    """
    
    section_name = 'i_section'
    
    def __init__(self, flange_width, flange_thickness, web_height, web_thickness):
        """
        Initialize I-section.
        
        Parameters
        ----------
        flange_width : float
            Width of each flange [m]
        flange_thickness : float
            Thickness of each flange [m]
        web_height : float
            Height of the web [m]
        web_thickness : float
            Thickness of the web [m]
        """
        if flange_width <= 0 or flange_thickness <= 0 or web_height <= 0 or web_thickness <= 0:
            raise ValueError("All dimensions must be positive.")
        
        self.flange_width = flange_width
        self.flange_thickness = flange_thickness
        self.web_height = web_height
        self.web_thickness = web_thickness
        
        super().__init__()
    
    def compute_properties(self):
        """
        Compute geometric properties for I-section using parallel axis theorem.
        
        The section consists of:
            - Top flange: width × flange_thickness
            - Web: web_thickness × web_height
            - Bottom flange: width × flange_thickness (same as top)
        
        Due to symmetry, the centroid is at mid-height.
        """
        b_f = self.flange_width
        t_f = self.flange_thickness
        h_w = self.web_height
        t_w = self.web_thickness
        
        # Total height
        h_total = 2 * t_f + h_w
        
        # Centroid location (from bottom, at mid-height due to symmetry)
        z_c = h_total / 2
        
        # Cross-sectional area
        A_top_flange = b_f * t_f
        A_web = t_w * h_w
        A_bottom_flange = b_f * t_f
        self._A = A_top_flange + A_web + A_bottom_flange
        
        # Second moment of area about y-axis (bending in x-z plane)
        # Using parallel axis theorem: I = I_centroid + A * d^2
        
        # Top flange
        I_top_flange_own = b_f * t_f**3 / 12  # About its own centroid
        z_top_centroid = h_total - t_f / 2  # Centroid of top flange from bottom
        d_top = z_top_centroid - z_c  # Distance from section centroid
        I_top_flange = I_top_flange_own + A_top_flange * d_top**2
        
        # Web
        I_web_own = t_w * h_w**3 / 12  # About its own centroid
        z_web_centroid = t_f + h_w / 2  # Centroid of web from bottom
        d_web = z_web_centroid - z_c  # Distance from section centroid
        I_web = I_web_own + A_web * d_web**2
        
        # Bottom flange
        I_bottom_flange_own = b_f * t_f**3 / 12  # About its own centroid
        z_bottom_centroid = t_f / 2  # Centroid of bottom flange from bottom
        d_bottom = z_bottom_centroid - z_c  # Distance from section centroid
        I_bottom_flange = I_bottom_flange_own + A_bottom_flange * d_bottom**2
        
        # Total I_y
        self._I_y = I_top_flange + I_web + I_bottom_flange
        
        # Second moment of area about z-axis (bending in x-y plane)
        # For I-section, this is simpler (no parallel axis needed if symmetric about z)
        I_z_top_flange = t_f * b_f**3 / 12
        I_z_web = h_w * t_w**3 / 12
        I_z_bottom_flange = t_f * b_f**3 / 12
        self._I_z = I_z_top_flange + I_z_web + I_z_bottom_flange
        
        # Section moduli
        # W_y = I_y / c_y, where c_y is distance to extreme fiber from centroid
        c_y_max = z_c  # Maximum distance to extreme fiber
        self._W_y = self._I_y / c_y_max
        
        # W_z = I_z / c_z, where c_z is distance to extreme fiber from centroid
        c_z_max = b_f / 2  # Maximum distance to extreme fiber (half of flange width)
        self._W_z = self._I_z / c_z_max
