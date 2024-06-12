# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:26:31 2024

@author: rensv
"""

# %% Import dependencies

import numpy as np
from abc import ABC, abstractmethod

# %% Class definitions
'''
Considerations:

total stiffness will be 6x6 (i.e. in 2D):
    [N_left, V_left, M_left, N_right, V_right, M_right] = K.[u_left, w_left, phi_left, u_right, w_right, phi_right]

where the contribution of the rod is as:
    K = [K_r_ll, 0, 0, K_r_lr, 0, 0; 
         0     , 0, 0, 0   , 0, 0;
         0     , 0, 0, 0   , 0, 0;
         K_r_rl, 0, 0, K_r_rr, 0, 0;
         0     , 0, 0, 0   , 0, 0;
         0     , 0, 0, 0   , 0, 0]
'''

# %% standard structural element class

class StructuralElement(ABC):
    '''
    Parent class of any other element, all functions are inherited by its child classes so add functions if necessary.
    '''
   
    element_name = None # initialise element name which should be overrided within the child class
    
    def __init__(self,dofs):
        
        
        self.dofs = dofs
        
        # check if an element name is assigned
        if self.element_name is None:
            raise ValueError(f"Class {self.__class__.__name__} must define a 'name' class attribute.")
            
    @abstractmethod
    def ElementWaveNumbers(self, omega):
        '''
        Must return the wavenumbers of the element
        '''
        pass
    
    @abstractmethod
    def LocalStiffness(self, omega):
        '''
        Must return the local stiffness matrix for the element.
        '''
        pass
    
    @abstractmethod
    def LocalElementDisplacements(self, u_nodes_global, omega, num_points):
        '''
        Must return the local element displacements as a function of global displacements, frequency, and number of evaluation points.
        
        Return should be:
            [displacement_dof1 
             displacement_dof2
             ...]
        '''
        pass
    
    @abstractmethod
    def Coefficients(self, u_node_local, omega):
        '''
        Must return the solution coefficients of the general solution of the element
        '''
        pass
    
    @abstractmethod
    def LocalDistributedLoad(self, q, omega):
        '''
        Must return the local distributed load vector.
        '''
        pass

# %% Element factory class
    
class ElementFactory:
    '''
    class that keeps track of all added elements such that they can easily be added to the Element class
    '''
    
    # initialise empty dictionary that contains the classes
    elements = {}

    @classmethod
    def RegisterElement(cls, element_class):
        """
        Registers a new element class in the factory under its 'name' attribute.
        """
        if not hasattr(element_class, 'element_name') or not isinstance(element_class.element_name, str):
            raise ValueError("Element classes must have a 'name' attribute of type str before registration.")
            
        cls.elements[element_class.element_name] = element_class
    
    @classmethod
    def CreateElement(cls, element_name, **kwargs):
        """
        Create an instance of an element by name. If the element type is not found,
        raise a ValueError with a message listing available element types.
        """
        if element_name not in cls.elements:
            available_types = ", ".join(cls.elements.keys())
            raise ValueError(f"No element registered with name: {element_name}. Available types are: {available_types}")
        element_class = cls.elements[element_name]
        return element_class(**kwargs)
    
    @classmethod
    def ListElementTypes(cls):
        """
        Return a list of currently registered element types.
        """
        return list(cls.elements.keys())
    
    @classmethod
    def ElementType(cls, name):
        """
        Decorator to register new element types by setting the 'name' attribute
        on the element class and then registering it.
        """
        def decorator(element_class):
            setattr(element_class, 'element_name', name)  # Explicitly set the name attribute
            cls.RegisterElement(element_class)   # Register the class using the newly set name
            return element_class
        return decorator
    
    