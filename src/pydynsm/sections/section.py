# -*- coding: utf-8 -*-
"""
Created on [Date]

@author: rensv
"""

import numpy as np
from abc import ABC, abstractmethod
import inspect


class Section(ABC):
    """
    Abstract base class for cross-sectional geometry.
    
    Sections store only geometric properties computed from dimensions.
    Material properties (E, rho, ksi, etc.) are provided separately when setting element type.
    
    Coordinate system:
        - x: horizontal/axial direction
        - z: vertical direction  
        - y: in-plane horizontal direction
    
    Geometric properties:
        - A: cross-sectional area [m^2]
        - I_y: second moment of area about y-axis (bending in x-z plane) [m^4]
        - I_z: second moment of area about z-axis (bending in x-y plane) [m^4]
        - W_y: section modulus about y-axis (W_y = I_y / c_y) [m^3]
        - W_z: section modulus about z-axis (W_z = I_z / c_z) [m^3]
    """
    
    section_name = None  # Should be overridden in child classes
    
    def __init__(self):
        """Initialize the section and compute geometric properties."""
        if self.section_name is None:
            raise ValueError(f"Class {self.__class__.__name__} must define a 'section_name' class attribute.")
        
        # Compute geometric properties from dimensions
        self.compute_properties()
    
    @abstractmethod
    def compute_properties(self):
        """
        Compute geometric properties (A, I_y, I_z, W_y, W_z) from section dimensions.
        
        This method should set:
            self.A: cross-sectional area
            self.I_y: second moment of area about y-axis
            self.I_z: second moment of area about z-axis
            self.W_y: section modulus about y-axis
            self.W_z: section modulus about z-axis
        """
        pass
    
    @property
    def A(self):
        """Cross-sectional area [m^2]"""
        return self._A
    
    @property
    def I_y(self):
        """Second moment of area about y-axis (bending in x-z plane) [m^4]"""
        return self._I_y
    
    @property
    def I_z(self):
        """Second moment of area about z-axis (bending in x-y plane) [m^4]"""
        return self._I_z
    
    @property
    def W_y(self):
        """Section modulus about y-axis [m^3]"""
        return self._W_y
    
    @property
    def W_z(self):
        """Section modulus about z-axis [m^3]"""
        return self._W_z


class SectionFactory:
    """
    Factory class for creating section instances by type name.
    Similar to ElementFactory pattern.
    """
    
    # Initialize empty dictionary that contains the section classes
    sections = {}
    
    @classmethod
    def RegisterSection(cls, section_class):
        """
        Registers a new section class in the factory under its 'section_name' attribute.
        Also checks the args and kwargs of the section type, and sets the args as required parameters.
        """
        if not hasattr(section_class, 'section_name') or not isinstance(section_class.section_name, str):
            raise ValueError("Section classes must have a 'section_name' attribute of type str before registration.")
        
        # Inspect the __init__ method and extract parameter information
        init_signature = inspect.signature(section_class.__init__)
        
        # Get required parameters (those without default values)
        required_params = [
            param.name for param in init_signature.parameters.values()
            if param.default == param.empty and param.name != 'self'
        ]
        
        # Get all parameters (both required and optional)
        all_params = [
            param.name for param in init_signature.parameters.values() if param.name != 'self'
        ]
        
        # Store the class along with its required and all parameters in the factory
        cls.sections[section_class.section_name] = {
            'class': section_class,
            'required_params': required_params,
            'all_params': all_params
        }
    
    @classmethod
    def CreateSection(cls, section_type, **props):
        """
        Create an instance of a section by type name. Validates that all required parameters are provided.
        
        Parameters
        ----------
        section_type : str
            Name of the section type (e.g., 'rectangle', 'circle')
        **props : dict
            Section dimensions (e.g., width, height, diameter, etc.)
        
        Returns
        -------
        Section
            Instance of the requested section type
        """
        if section_type not in cls.sections:
            available_types = ", ".join(cls.sections.keys())
            raise ValueError(f"No section registered with name: {section_type}. Available types are: {available_types}")
        
        # Retrieve the stored class and parameter info
        section_info = cls.sections[section_type]
        section_class = section_info['class']
        required_params = section_info['required_params']
        all_params = section_info['all_params']
        
        # Validate the provided props against required and all params
        missing_params = [param for param in required_params if param not in props]
        extra_params = [param for param in props if param not in all_params]
        
        if missing_params:
            raise ValueError(f"Missing required parameters for section '{section_type}': {', '.join(missing_params)}")
        if extra_params:
            raise ValueError(f"Unexpected parameters provided for section '{section_type}': {', '.join(extra_params)}")
        
        return section_class(**props)
    
    @classmethod
    def ListSectionTypes(cls):
        """
        Return a list of currently registered section types,
        printed neatly in the terminal.
        """
        print("Available Section Types:")
        for name in sorted(cls.sections.keys()):
            print(f"  - {name}")
    
    @classmethod
    def SectionType(cls, name: str):
        """
        Decorator to register new section types by setting the 'section_name' attribute
        on the section class and then registering it.
        
        Parameters
        ----------
        name : str
            Name of the section type (e.g., 'rectangle', 'circle')
        
        Returns
        -------
        decorator : function
            Decorator function that registers the section class
        """
        def decorator(section_class):
            # Set section class name
            setattr(section_class, 'section_name', name)
            
            # Register the section class
            cls.RegisterSection(section_class)
            return section_class
        return decorator
