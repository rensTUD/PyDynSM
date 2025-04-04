# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:34:56 2024

@author: rensv
"""

# %% Import dependencies
import numpy as np
from collections import defaultdict
from typing import Optional, Dict
from .dofs import DOF
from .dofs import DOFContainer


# %% class definition
class Node:
    """
    A class to represent a node in a structural system.

    Attributes
    ----------
    nn : int
        Class variable to count the number of nodes.
    dof_configurations : dict
        Dictionary containing DOF configurations for different node types.
    x : float
        x-coordinate of the node.
    z : float
        z-coordinate of the node.
    y : float
        y-coordinate of the node, applicable for 3D configurations.
    config : str
        Configuration type of the node (e.g., '2D', '3D', '2D_torsion').
    dofs : dict
        Dictionary of degrees of freedom for the node with DOF names as keys and their statuses as values.
    nodal_forces : list
        List to store forces applied to the node.
    id : int
        Unique identifier for the node.
    connected_elements : list
        List of elements connected to this node.
    """
    
    # initialise number of nodes (nn) as class-variable, to be updated by every new class object
    nn   = 0
    
    # Node configurations are ordered as: 'Name': [['linear dofs'],['rotational dofs']]. Linear dofs also define the amount of coordinates we can assign (thus 2D or 3D basically)
    dof_configurations = {
        '2D': ['x', 'z','phi_y'],
        '3D': ['x', 'z', 'y', 'phi_x', 'phi_z', 'phi_y'],
        '2D_torsion': ['x', 'z', 'phi_x', 'phi_y'] 
        }
    
    def __init__(self,x,z,y=0,config='2D', dof_config=None):
        """
        Initializes the Node with its coordinates and DOF configuration.

        Parameters
        ----------
        x : float
            x-coordinate of the node.
        z : float
            z-coordinate of the node.
        y : float, optional
            y-coordinate of the node, default is 0 (for 2D configurations).
        config : str, optional
            Configuration type of the node, default is '2D'.
        """
        
        # dofs of the node: {dof: value} - None = Free, 0 = Fixed, value = presribed
        self.config = config
        
        # Load any default dof_config or apply custom one SHOULD THINK WHETHER THIS IS USEFULL OR JUST LEAVE IT FULLY AUTOMATIC BASED ON APPLIED SECTIONS
        self.dof_config = dof_config if dof_config else self.dof_configurations[config]
        
        # initialise the DOFcontainer
        self.dof_container = DOFContainer()
        
        # Set up DOFs based on configuration
        for dof_name in self.dof_config:
            self.dof_container.set_dof(dof_name)
                
        # set the location in space
        self.x     = x
        self.z     = z
        self.y     = y
        
        # node name
        self.id = Node.nn
        Node.nn   += 1 # increment the class variables
        
        # Initialise empty lists needed for later
        self.nodal_loads = defaultdict(dict)
        self.connected_elements = []
        
    def get_dof_indices(self):
        """Returns a list of global DOF indices for the node."""
        return [dof.index for dof in self.dof_container.dofs.values() if dof.index is not None]

    def update_element_dof(self,changes=None):
        """
        Updates the DOFs of connected elements based on changes in the node's DOFs.

        Parameters
        ----------
        changes : dict, optional
            Dictionary of changes in DOFs, default is None.
        """
        if changes:
            for element in self.connected_elements:
                element.update_node_dofs(self,changes)
                
    def apply_dof_change_to_elements(self, dof_name, value):
        """
        Applies the change in a specific DOF value to the connected elements,
        mapping the global DOF to the corresponding local and global DOFs in each element.
    
        Parameters
        ----------
        dof_name : str
            The name of the global DOF that has changed (e.g., 'x', 'phi_y').
        value : Any
            The value to set for the specified DOF.
        """
        # Get the DOF from the container
        dof = self.dof_container.get_dof(dof_name)
        
        # Ensure the DOF exists and is part of the node's current configuration
        if dof is None:
            raise ValueError(f"DOF '{dof_name}' is not part of the node's current configuration.")
        
        # Propagate the change to all connected elements if present
        if not self.connected_elements:
            return
        else:
            for element in self.connected_elements:
                element.apply_global_dof_change(self, dof_name, value)
    
        print(f"Global DOF '{dof_name}' with value '{value}' applied to connected elements.")               

    def connect_element(self, element):
        """
        Connects an element to this node.

        Parameters
        ----------
        element : Element
            The element to be connected to this node.
        """
        if element not in self.connected_elements:
            self.connected_elements.append(element)
            
    def fix_node(self, *dofs):
        """
        Fixes specified DOFs of the node to zero.

        Parameters
        ----------
        *dofs : str
            DOFs to be fixed.
        """
    
        self.prescribe_node(**{dof: 0 for dof in dofs})
        
    def free_node(self, *dofs):
        """
        Frees specified DOFs of the node (sets them to None).

        Parameters
        ----------
        *dofs : str
            DOFs to be freed.
        """
        
        self.prescribe_node(**{dof: None for dof in dofs})
       
    def prescribe_node(self, **dofs):
        """
        Sets specified DOFs of the node to given values.

        Parameters
        ----------
        **dofs : dict
            DOFs to be prescribed with their values.
        
        Raises
        ------
        ValueError
            If any of the specified DOFs are not in the current configuration.
        """

        try:
            for dof_name, value in dofs.items():
                if not self.dof_container.has_dof(dof_name):
                    raise ValueError(f"DOF '{dof_name}' is not available in the current configuration. Available DOFs: {self.dof_config}")
                
                dof = self.dof_container.get_dof(dof_name)
                if dof.value != value:
                    # set value if is not the same
                    dof.value = value
                    # also apply the change to the connected elements
                    self.apply_dof_change_to_elements(dof_name,value)
        except ValueError as e:
            print(e)
            return
                
    def add_load(self, **loads):
        '''
        Adds distributed loads to the Node.
    
        Parameters
        ----------
        **loads : dict
            Keyword arguments where the key is the DOF (degree of freedom) and the value is the load magnitude. 
            e.g. node.add_load(x=load_x, y=load_y)
        '''
        for dof, load in loads.items():
            if self.dof_container.has_dof(dof):
                if dof in self.nodal_loads:
                    print(f'Load for DOF {dof} is already present and will be overwritten. Old value: {self.nodal_loads[dof]}, New value: {load}')
                self.nodal_loads[dof] = load
            else:
                print(f'Node {self.node_id} does not contain DOF {dof}. Load has not been added.')
                
    def remove_load(self, dof):
        '''
        Removes the applied load for a specific degree of freedom (DOF) from the node.
    
        Parameters
        ----------
        dof : str
            The degree of freedom (e.g., 'x', 'y', 'z', 'phi_x', etc.) for which the load should be removed.
        '''
        if dof in self.nodal_loads:
            removed_load = self.nodal_loads.pop(dof)
            print(f'Load for DOF {dof} (value: {removed_load}) has been removed.')
        else:
            print(f'No load found for DOF {dof}. Nothing was removed.')
                                   
    def get_coords(self):
        """
        Returns the coordinates of the node based on its DOF configuration.

        Returns
        -------
        tuple
            Coordinates of the node (x, z, y) if applicable.
        """
        coords = [self.x, self.z, self.y]

        return tuple(coords)
