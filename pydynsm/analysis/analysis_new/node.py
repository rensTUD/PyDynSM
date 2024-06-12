# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:34:56 2024

@author: rensv
"""

# %% Import dependencies
import numpy as np

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
        '2D': [['x', 'z'],['phi_y']],
        '3D': [['x', 'z', 'y'], ['phi_x', 'phi_z', 'phi_y']],
        '2D_torsion': [['x', 'z'], ['phi_x', 'phi_y']] 
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
        # self.dofs = {dof: None for dof in self.dof_configurations[config]}  # Initialize all DOFs to None
        self.dofs = {dof: None for dof_list in self.dof_configurations[config] for dof in dof_list}

        # location in space
        self.x     = x
        self.z     = z
        self.y     = y
        
        # initialise empty forces lists
        self.nodal_forces = []
        
        # node name
        self.id = Node.nn
        
        # empty connected_elements list
        self.connected_elements = []
        
        # increment the class variables
        Node.nn   += 1

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
        changes = {}
        try:
            for dof, value in dofs.items():
                if dof not in self.dofs:
                    raise ValueError(f"DOF '{dof}' is not available in the current configuration. Available DOFs: {self.dof_configurations[self.config]}")
                if self.dofs[dof] != value:
                    self.dofs[dof] = value
                    changes[dof] = value
        except ValueError as e:
            print(e)
            return
        
        if changes:
            for element in self.connected_elements:
                element.update_node_dofs(self,changes)
    
    def get_coords(self):
        """
        Returns the coordinates of the node based on its DOF configuration.

        Returns
        -------
        tuple
            Coordinates of the node (x, z, y) if applicable.
        """
        coords = [self.x, self.z]
        if 'y' in self.dof_configurations[self.config][0]:
            coords.append(self.y)
        return tuple(coords)
