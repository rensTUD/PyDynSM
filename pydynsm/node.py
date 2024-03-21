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
    The Node class is used to store node information and keep track of the total number of 
    Degrees of Freedom (DOFs) of the problem. It introduces automatic bookkeeping in its 
    initialization, which efficiently keeps track of which DOFs belong to which nodes. This 
    makes it easier to assemble matrices from multiple elements.

    Attributes:
        x (float): The x-coordinate of the node.
        z (float): The z-coordinate of the node.
        p (numpy.array):  The load vector of the node.
        dofs (list): The Degrees of Freedom associated with the node.

    Methods:
        clear(): Clears the counting of degrees of freedom and number of nodes.
        __init__(x, z): The constructor for Node class.
        add_load(p): Adds the given loads to the node.
        get_coords(): Returns the coordinates of the node.
        __str__(): Returns a string representation of the node.
    """
    
    # determine class variables shared through all instances of the Node class
    ndof = 0
    nn   = 0
    constrained_dofs_global = []
    
    # mapping of dofs
    dof_map = {'x': 0, 'z': 1, 'phi': 2}
    
    def clear():
        """
        Clears the counting of degrees of freedom and number of nodes.

        This method resets the class-level counters for degrees of freedom and number of nodes. 
        It should be used when you want to start a new problem from scratch.
        """
        Node.ndof = 0
        Node.nn = 0
        Node.constrained_dofs_global.clear()
        
    def __init__(self, x, z, assembler=None, x_fixed=False, z_fixed=False, phi_fixed=False): 
        """
        The constructor for Node class.

        Parameters:
            x (float):        The x-coordinate of the node.
            z (float):        The z-coordinate of the node.
            p (numpy.array):  The load vector of the node.
            dofs (list):      The Degrees of Freedom (u (in direction of x), w (in direction of z), phi (from z to x)) associated with the node.
        """
        
        # location in space
        self.x     = x
        self.z     = z
        
        # initialise empty forces lists
        self.p     = np.zeros(3,complex)
        
        # numbers of dofs [x,z,phi]
        self.dofs  = [Node.ndof, Node.ndof+1, Node.ndof+2]
        
        # empty list with the fixed dofs and their values: [[dof1, value], [dof2, value]]
        self.constrained_dofs = []        
        
        
        # check whether node is being initialised with some fixation
        if x_fixed:
            self.fix_dof('x')
        if z_fixed:
            self.fix_dof('z')
        if phi_fixed:
            self.fix_dof('phi')
        
        # increment the class variables
        Node.ndof += 3
        Node.nn   += 1
        
        # node name
        self.name = f"Node {Node.nn}"
        
        # add the node to the assembler if given
        if assembler is not None:
            assembler.RegisterNode(self)

    def add_load(self, p):
        """
        Adds the given loads to the node.

        The load is a vector p, which includes the load in the x and y direction and a moment. 
        These loads are added to the existing loads of the node.

        Parameters:
            p (numpy.array): A vector containing the load in the x direction, the load in the y direction, 
                             and the moment. 
        """
        self.p += p

    def get_coords(self):
        """
        Returns the coordinates of the node.

        Returns:
           numpy.array: An array containing the x and z coordinates of the node.
        """
        return np.array([self.x, self.z])
    
    def fix_dofs(self, *args):
        """
        Fixes specific degrees of freedom for this node, with given values or default.
    
        Parameters:
        *args: A mix of strings and tuples, where strings represent the DOF to be fixed with 
               a default value, and tuples represent the DOF and its specific value to be fixed.
        """
        for arg in args:
            if isinstance(arg, str):  # Single DOF with default value
                dof = arg
                value = 0  # Default value
                self.fix_dof(dof, value)
            elif isinstance(arg, list) and len(arg) == 2:
                dof, value = arg
                self.fix_dof(dof, value)
            else:
                raise ValueError("Invalid argument format. Use 'dof' or ('dof', value).")
    
    def fix_dof(self, dof, value):
        """
        Helper method to fix a single degree of freedom with a given value.
    
        Parameters:
        dof (str): The identifier of the degree of freedom to be fixed ('x', 'z', or 'phi').
        value (float): The value to fix the degree of freedom to.
        """
        if dof in Node.dof_map:
            dof_index = self.dofs[Node.dof_map[dof]]
            if (dof_index, value) not in Node.constrained_dofs_global:
                Node.constrained_dofs_global.append([dof_index, value])
                self.constrained_dofs.append([dof_index, value])
                print(f"DOF {dof} for node at ({self.x}, {self.z}) is now fixed with value {value}.")
            else:
                print(f"DOF {dof} for node at ({self.x}, {self.z}) is already fixed.")
        else:
            raise ValueError(f"Invalid DOF '{dof}'. Please use 'x', 'z', or 'phi'.")

    def fix_node(self):
        """
        Constrain all degrees of freedom for this node.
        """
        self.fix_dofs('x', 'z', 'phi')

    def __str__(self):
        """
        Returns a string representation of the node.

        Returns:
            str: A string representation of the node.
        """
        return f"This node has:\n - x coordinate={self.x},\n - z coordinate={self.z},\n - degrees of freedom={self.dofs},\n - load vector={self.p})"
