# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:34:56 2024

@author: rensv
"""

# %% Import dependencies
import numpy as np
from collections import defaultdict

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
        self.dofs = {dof: None for dof_list in self.dof_config for dof in dof_list}
        
        self.dof_indices = {}
        
        # location in space
        self.x     = x
        self.z     = z
        self.y     = y
        
        # initialise empty forces dictionary
        self.element_loads = defaultdict(dict)
        
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
                
    def apply_dof_change_to_elements(self, changed_dof):
            """
            Applies the change in a specific global DOF to the connected elements,
            mapping the global DOF to the corresponding local and global DOFs in each element.
    
            Parameters
            ----------
            changed_dof : str
                The name of the global DOF that has changed (e.g., 'x', 'phi_y').
            """
            # Check if the changed_dof is part of the node's current DOFs
            if changed_dof not in self.dofs:
                raise ValueError(f"DOF '{changed_dof}' is not part of the node's current configuration.")
            
            # Get the changed DOF value
            changed_value = self.dofs[changed_dof]
    
            # Propagate the change to all connected elements
            for element in self.connected_elements:
                element.apply_global_dof_change(self, changed_dof, changed_value)
    
            print(f"Global DOF '{changed_dof}' change applied to connected elements.")                

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
                    raise ValueError(f"DOF '{dof}' is not available in the current configuration. Available DOFs: {self.dof_config}")
                if self.dofs[dof] != value:
                    self.dofs[dof] = value
                    changes[dof] = value
        except ValueError as e:
            print(e)
            return
        
        if changes:
            for element in self.connected_elements:
                element.update_node_dofs(self,changes)

    # def AddLoad(self, p):
    #     """
    #     Adds the given loads to the node.

    #     The load is a vector p, which includes the load in the x and y direction and a moment. 
    #     These loads are added to the existing loads of the node.

    #     Parameters:
    #         p (numpy.array): A vector containing the load in the x direction, the load in the y direction, 
    #                          and the moment. 
    #     """
    #     self.nodal_forces.append(p)
        
    def AddLoad(self,**loads):
        '''
        Adds distributed loads to the element.
    
        Parameters
        ----------
        **loads : dict
            Keyword arguments where the key is the DOF (degree of freedom) and the value is the load magnitude.
        '''
        
        # Assign unpacked loads to the element's list of loads
        for dof, load in loads.items():
            if dof not in self.element_loads.keys():
                self.nodal_loads[dof] = load
            else:
                print(f'Load already added for DOF {dof}')        

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
