# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:26:31 2024

@author: rensv
"""

# %% Import dependencies

import numpy as np
from abc import ABC, abstractmethod
import inspect
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
        
        # assign dofs in text
        self.dofs = dofs
        
        # check if an element name is assigned
        if self.element_name is None:
            raise ValueError(f"Class {self.__class__.__name__} must define a 'name' class attribute.")
        
        
    
    def GetNodalDofs(self, dofs, Ndof):
        '''
        Translates the dofs of an element to the nodal dofs;
        
        Example:
            Rod element has dof = 0 (axial only)
            This will give: nodal_dofs = [0,3]
        '''
        
        nodal_dofs = []
        # Process each DOF to map it to the current and next node's corresponding DOFs
        for dof in dofs:
            # Add the DOF for the current node
            nodal_dofs.append(dof)
        # After processing each DOF for the current node, add their counterparts in the next node
        for dof in dofs:
            nodal_dofs.append(dof + Ndof)
                
        return nodal_dofs
    
    def FullStiffness(self, omega):
        '''
        Function that assembles the full stiffness matrix based on the local stiffness matrix. 
        
        For example, it will translate the 2x2 K_local of the rod to the full 6x6 matrix which it is in 2D
        
        '''
        
        # initialise 6x6 empty complex matrix
        K_full = np.zeros((6,6), complex)
        
        # calculate local stiffness matrix
        K_local = self.LocalStiffness(omega)
        
        # assign to correct locations in full matrix
        for i, row in enumerate(self.nodal_dofs):
            for j, col in enumerate(self.nodal_dofs):
                K_full[row, col] = K_local[i, j]
                
        return K_full 
    
    def FullDistributedLoad(self, q, omega):
        
        # initialise 6x1 empty complex vector 
        q_full = np.zeros(6,complex)
        
        # calculate local load vector
        q_local = self.LocalDistributedLoad(q, omega)
        
        # assign to correct locations in full vector
        for i, row in enumerate(self.nodal_dofs):
            q_full[row] = q_local[i]
        
        return q_full
    
    def FullDisplacement(self, u_node, omega):
        
        # initialise 6x1 empty complex vector 
        u_full = np.zeros(6,complex)
        
        # calculate local load vector
        u_local = self.LocalDisplacements(u_node, omega)
        
        # assign to correct locations in full vector
        for i, row in enumerate(self.nodal_dofs):
            u_full[row] = u_local[i]
        
        return u_full
    
    def FullElementDisplacements(self, u_nodes_global, omega, num_points):
        
        # intilise full u_elem
        u_elem_full = [None] * self.Ndof
        
        # calculate local u_elem
        u_elem_local = self.LocalElementDisplacements(u_nodes_global, omega, num_points)
        
        # assign
        for i, dof in enumerate(self.dofs):
            u_elem_full[dof] = u_elem_local[i]
        
        return u_elem_full
    
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
        Also checks the args and kwargs of the element type, and sets the args as required parameters.
        """
        if not hasattr(element_class, 'element_name') or not isinstance(element_class.element_name, str):
            raise ValueError("Element classes must have a 'name' attribute of type str before registration.")
        
        # Inspect the __init__ method and extract parameter information
        init_signature = inspect.signature(element_class.__init__)
        
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
        cls.elements[element_class.element_name] = {
            'class': element_class,
            'required_params': required_params,
            'all_params': all_params
        }
        
    
    @classmethod
    def CreateElement(cls, element_name, **kwargs):
        """
        Create an instance of an element by name. Validates that all required parameters are provided.
        """
        if element_name not in cls.elements:
            available_types = ", ".join(cls.elements.keys())
            raise ValueError(f"No element registered with name: {element_name}. Available types are: {available_types}")
               
        # Retrieve the stored class and parameter info
        element_info = cls.elements[element_name]
        element_class = element_info['class']
        required_params = element_info['required_params']
        all_params = element_info['all_params']
        
        # Validate the provided kwargs against required and all params
        missing_params = [param for param in required_params if param not in kwargs]
        extra_params = [param for param in kwargs if param not in all_params]

        if missing_params:
            raise ValueError(f"Missing required parameters for element '{element_name}': {', '.join(missing_params)}")
        if extra_params:
            raise ValueError(f"Unexpected parameters provided for element '{element_name}': {', '.join(extra_params)}")
                
        return element_class(**kwargs)
        
    
    @classmethod
    def ListElementTypes(cls):
        """
        Return a list of currently registered element types.
        """
        return list(cls.elements.keys())
    
    @classmethod
    def ElementType(cls, name: str):#, required_parameters: list | None = None):
        """
        Decorator to register new element types by setting the 'name' attribute
        on the element class and then registering it.
        """
        def decorator(element_class):
            # set element class name
            setattr(element_class, 'element_name', name)  # Explicitly set the name attribute
            # # set required parameters if present
            # if required_parameters is not None:
            #     setattr(element_class, 'required_parameters', required_parameters)
            # register the element class
            cls.RegisterElement(element_class)   # Register the class using the newly set name
            return element_class
        return decorator
    
    