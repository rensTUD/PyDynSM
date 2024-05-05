# -*- coding: utf-8 -*-
"""
Created on Sun May  5 18:21:53 2024

@author: rensv
"""



import numpy as np
from scipy.linalg import null_space
from scipy.linalg import inv

# %% universal decorator:
def decorate_after(function):
    """Decorator factory that takes a post-method function to call after the decorated method."""
    def decorator(method):
        def wrapper(self, *args, **kwargs):
            result = method(self, *args, **kwargs)
            if result is not None:
                function(self, result)
            else:
                function(self)
            return result
        return wrapper
    return decorator

# %% node class
class Node:
    
    nn   = 0
    
    def __init__(self,x,z):
        
        # location in space
        self.x     = x
        self.z     = z
        
        # initialise empty forces lists
        self.nodal_forces     = []
        
        # node name
        self.id = Node.nn
        
        # empty connected_elements list
        self.connected_elements = []
        
        # increment the class variables
        Node.nn   += 1
        
        # dofs of the node: {dof: value} - None = Free, 0 = Fixed, value = presribed
        self.dofs = {'x': None, 'y': None, 'z': None} # initialise as free nodes

    def update_element_dof(self,changes=None):
        if changes:
            for element in self.connected_elements:
                element.update_node_dofs(self,changes)
        else:
            print('No changes to update')                

    def connect_element(self, element):
        if element not in self.connected_elements:
            self.connected_elements.append(element)
            
    @decorate_after(update_element_dof)
    def fix_node(self, *dofs):
        """
        Fix specified DOFs to zero and return changes.
        """
        changes = {}
        for dof in dofs:
            if self.dofs[dof] != 0:
                self.dofs[dof] = 0
                changes[dof] = 0
        return changes if changes else None  # Return None if no changes
        
    @decorate_after(update_element_dof)
    def free_node(self, *dofs):
        """
        Set specified DOFs to free (None) and return changes.
        """
        changes = {}
        for dof in dofs:
            if self.dofs[dof] != None:
                self.dofs[dof] = None
                changes[dof] = None
        return changes if changes else None
    
    @decorate_after(update_element_dof)     
    def prescribe_node(self, **dofs):
        """
        Set specified DOFs to given values and return changes.
        """
        changes = {}
        for dof, value in dofs.items():
            if self.dofs[dof] != value:
                self.dofs[dof] = value
                changes[dof] = value
        return changes if changes else None
    

    
    
nodes = [Node(i,0) for i in range(3)]      
nodes.append(Node(1,1))
node1 = nodes[0]
# # text fix 'x' and 'y' of node 1
# node1.fix_node('x','y')
# print(node1.dofs)
# # text free 'x'
# node1.free_node('x')
# print(node1.dofs)
# # test set value
# node1.prescribe_node(x=1,y=2)
# print(node1.dofs)

#%% element class
class Element:
    
    ne = 0
    
    def __init__(self, nodes):
        
        self.id = Element.ne
        
        self.nodes = nodes
        
        for node in nodes:
            node.connect_element(self)  # Connect this element to its nodes
        
        # assign same dofs of node on initialisation
        self.dofs = {node.id: node.dofs for node in nodes}
        
        # increment element number
        Element.ne += 1

    def fix_dof(self, node, *dofs):
        '''
        Fixed dof at a certain node on element level
        '''
        if node in self.nodes:
            for dof in dofs:
                self.dofs[node.id][dof] = 0
                
    def free_dof(self, node, *dofs):
        '''
        Fixed dof at a certain node on element level
        '''
        if node in self.nodes:
            for dof in dofs:
                self.dofs[node.id][dof] = None 
                
    def prescribe_dof(self, node, **dofs):
        
        if node in self.nodes:
            for dof, value in dofs.items():
                self.dofs[node.id][dof] = value                

    def update_node_dofs(self, node, changes):
        if node in self.nodes:    
            for dof, value in changes.items():           
                self.dofs[node.id][dof] = value

# test create elements        
elements = [Element([nodes[0], nodes[1]]),
            Element([nodes[1], nodes[2]]),
            Element([nodes[3],nodes[1]])]
element1 = elements[0]

# # test fix local 'x' and 'y' of the element
# element1.fix_dof(node1,'x','y')
# print(element1.dofs)
# # test free local 'x' 
# element1.free_dof(node1,'x')
# print(element1.dofs)
# # test set value
# element1.prescribe_dof(node1, x=1,y=2)
# print(element1.dofs)

# set x,z of first and z of last node fixed
node1.fix_node('x','z')
print(node1.dofs)

node3 = nodes[2]
node3.fix_node('z')
print(node3.dofs)

node4 = nodes[3]
node4.fix_node('z')
print(node4.dofs)

# %%

def find_unique_redundant_dofs(B):
    num_dofs = B.shape[1]
    redundant_dofs = set()
    unique_dofs = set()

    for row in B:
        positives = np.where(row > 0)[0]
        negatives = np.where(row < 0)[0]

        for pos in positives:
            unique_dofs.add(pos)
        for neg in negatives:
            redundant_dofs.add(neg)

    # Only columns that are completely zero and not identified in negatives are truly unique
    zero_dofs = {i for i in range(num_dofs) if np.all(B[:, i] == 0)}
    unique_dofs = (unique_dofs | zero_dofs) - redundant_dofs

    return sorted(unique_dofs), sorted(redundant_dofs)

def calculate_L(B, unique_dofs, redundant_dofs):
    B_r = B[:, redundant_dofs]
    B_u = B[:, unique_dofs]

    try:
        B_r_inv = inv(B_r)
        L_lower = -B_r_inv @ B_u
        L = np.vstack((np.eye(len(unique_dofs)), L_lower))
        return L
    except np.linalg.LinAlgError:
        raise ValueError("Matrix B_r is singular and cannot be inverted.")

def assign_dof_indices(elements):
    dof_indices = {}
            
    for element in elements:
        base_index = 6 * element.id
        for idx, node in enumerate(element.nodes):
            # Assign indices for each DOF of this node in the current element
            offset = base_index + 3 * idx # Each element contributes 6 DOFs, 3 per node
            dof_indices[(node.id, element.id)] = {
                'x': offset,
                'y': offset + 1,
                'z': offset + 2
            }
            
    return dof_indices

def build_matrix_B(nodes, elements, dof_indices):
    num_dof = 6 * len(elements)  # Assuming 2 nodes per element, each with 3 DOFs
    B = np.zeros((0, num_dof))

    # Enforce compatibility conditions at each node
    for node in nodes:
        if len(node.connected_elements) > 1:  # Only consider nodes connected to more than one element
            for dof in ['x', 'y', 'z']:
                all_dofs = []
                # Collect all DOF indices and values for this node from each connected element
                for element in node.connected_elements:
                    if (node.id, element.id) in dof_indices:
                        index = dof_indices[(node.id, element.id)][dof]
                        dof_value = element.dofs[node.id][dof]  # Retrieve the DOF value from the element's settings
                        all_dofs.append((index, dof_value))

                # Use the first DOF as the reference for compatibility
                if all_dofs:
                    primary_index, primary_value = all_dofs[0]
                    for index, value in all_dofs[1:]:
                        if value == primary_value:  # Compare values, enforce compatibility if they match
                            row = np.zeros(num_dof)
                            row[primary_index] = 1  # Positive sign for the primary dof
                            row[index] = -1  # Negative sign for the redundant dof
                            B = np.vstack([B, row])  # Add this new row to the matrix B

    return B

def connectivity(nodes, elements):
    dof_indices = assign_dof_indices(elements)
    B = build_matrix_B(nodes, elements, dof_indices)
    unique_dofs, redundant_dofs = find_unique_redundant_dofs(B)
    L = calculate_L(B, unique_dofs, redundant_dofs)

    return B, L, dof_indices, unique_dofs, redundant_dofs



# Usage
B, L, dof_indices, unique_dofs, redundant_dofs = connectivity(nodes, elements)
print("B:", B)
print("L:", L)



