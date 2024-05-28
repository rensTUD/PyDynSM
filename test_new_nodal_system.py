# -*- coding: utf-8 -*-
"""
Created on Sun May  5 18:21:53 2024

@author: rensv
"""


import numpy as np
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
    
    def __init__(self,x,z,y=0,config='2D'):
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
        # else:
        #     print('No changes to update')                

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
            
    # @decorate_after(update_element_dof)
    def fix_node(self, *dofs):
        """
        Fixes specified DOFs of the node to zero.

        Parameters
        ----------
        *dofs : str
            DOFs to be fixed.
        """
    
        self.prescribe_node(**{dof: 0 for dof in dofs})
        
    # @decorate_after(update_element_dof)
    def free_node(self, *dofs):
        """
        Frees specified DOFs of the node (sets them to None).

        Parameters
        ----------
        *dofs : str
            DOFs to be freed.
        """
        
        self.prescribe_node(**{dof: None for dof in dofs})
    
    # @decorate_after(update_element_dof)     
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
    

    
    
nodes = [Node(i,0) for i in range(3)]      
nodes.append(Node(1,1))
node1 = nodes[0]
node2 = nodes[1]
node3 = nodes[2]
node4 = nodes[3]
#%%

#%%% test fix 'x' and 'y' of node 1, should throw error
node1.fix_node('x','y')
print(node1.dofs)

# now just add z
node1.fix_node('z')
print(node1.dofs)

#%%% test free 'x'
node1.free_node('x')
print(node1.dofs)

#%%% test set value
node1.prescribe_node(x=1,y=2)
print(node1.dofs)

#%% element class
class Element:
    """
    A class to represent an element in a structural system.

    Attributes
    ----------
    ne : int
        Class variable to count the number of elements.
    id : int
        Unique identifier for the element.
    nodes : list of Node
        List of nodes connected to this element.
    dofs : dict
        Dictionary of DOFs for each node in the element.
    """
    
    # initialise number of elements (ne) as class-variable, to be updated by every new class object
    ne = 0
    
    def __init__(self, nodes):
        """
        Initializes the Element with its connected nodes.

        Parameters
        ----------
        nodes : list of Node
            List of nodes connected to this element.
        """
        
        # Element number
        self.id = Element.ne
        Element.ne += 1 # increment element number
        
        # reference the nodes
        self.nodes = nodes
        
        # connect the element to the node
        for node in nodes:
            node.connect_element(self)  # Connect this element to its nodes
        
        # assign same dofs of node on initialisation
        self.dofs = {node.id: node.dofs.copy() for node in nodes}
        
    def fix_dof(self, node, *dofs):
        """
        Fixes specified DOFs at a given node in the element.

        Parameters
        ----------
        node : Node
            The node where DOFs are to be fixed.
        *dofs : str
            DOFs to be fixed at the node.

        Raises
        ------
        ValueError
            If any of the specified DOFs are not in the node's current configuration.
        """
        
        self.prescribe_dof(node, **{dof: 0 for dof in dofs})
                
    def free_dof(self, node, *dofs):
        """
        Frees specified DOFs at a given node in the element.

        Parameters
        ----------
        node : Node
            The node where DOFs are to be freed.
        *dofs : str
            DOFs to be freed at the node.

        Raises
        ------
        ValueError
            If any of the specified DOFs are not in the node's current configuration.
        """
        
        self.prescribe_dof(node, **{dof: None for dof in dofs})
                
    def prescribe_dof(self, node, **dofs):
        """
        Prescribes specified DOFs at a given node in the element to given values.

        Parameters
        ----------
        node : Node
            The node where DOFs are to be prescribed.
        **dofs : dict
            DOFs and their values to be prescribed at the node.

        Raises
        ------
        ValueError
            If any of the specified DOFs are not in the node's current configuration.
        """
        try:
            for dof, value in dofs.items():
                if dof not in node.dofs:
                    raise ValueError(f"DOF '{dof}' is not available in the node's current configuration. Available DOFs: {node.dof_configurations[node.config]}")
                self.dofs[node.id][dof] = value
        except ValueError as e:
            print(e)              

    def update_node_dofs(self, node, changes):
        """
        Updates the DOFs of a given node in the element based on changes.

        Parameters
        ----------
        node : Node
            The node whose DOFs are to be updated.
        changes : dict
            Dictionary of DOFs and their new values.
        """
        if node in self.nodes:    
            for dof, value in changes.items():           
                self.dofs[node.id][dof] = value

# test create elements        
elements = [Element([nodes[0], nodes[1]]),
            Element([nodes[1], nodes[2]]),
            Element([nodes[3],nodes[1]])]
element1 = elements[0]
element2 = elements[1]
element3 = elements[2]


#%%
#%%% test fix local 'x' and 'y' of the element
element1.fix_dof(node1,'x','y')
print(element1.dofs)
#%%% test free local 'x' 
element1.free_dof(node1,'x')
print(element1.dofs)
#%%% test set value
element1.prescribe_dof(node1, x=1,y=2)
print(element1.dofs)
#%%
# set x,z of first and z of last node fixed
node1.fix_node('x','z')
print(node1.dofs)


node3.fix_node('z')
print(node3.dofs)


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
    global_index = 0  # Start a global index counter for all DOFs in the system

    for node in nodes:
        for element in node.connected_elements:
            # Initialize an entry for this node within this element
            dof_indices[(node.id, element.id)] = {}
            # Use the DOF configuration for this particular node (based on the node's own configuration)
            for dof in node.dofs.keys():
                dof_indices[(node.id, element.id)][dof] = global_index
                global_index += 1  # Increment global index for each DOF

    return dof_indices

def build_matrix_B(nodes, elements, dof_indices):
    # Calculate the total number of DOFs
    num_dof = sum(len(node.dofs) for element in elements for node in element.nodes )  # Dynamic count of all DOFs
    B = np.zeros((0, num_dof))

    # Enforce compatibility conditions at each node
    for node in nodes:
        if len(node.connected_elements) > 1:  # Only consider nodes connected to more than one element
            # We need to handle this correctly using (node.id, element.id)
            connected_indices = []
            for element in node.connected_elements:
                if (node.id, element.id) in dof_indices:
                    connected_indices.append((element, dof_indices[(node.id, element.id)]))

            # Proceed only if there are connected elements with valid DOF indices
            if connected_indices:
                # Assuming all elements have the same DOFs for nodes
                for dof in node.dofs.keys():
                    all_dofs = []
                    for element, indices in connected_indices:
                        if dof in indices:
                            index = indices[dof]
                            dof_value = element.dofs[node.id][dof] # Node's DOF value from the element itself
                            all_dofs.append((index, dof_value))
                    
                    # Use the first DOF as the reference for compatibility
                    if all_dofs:
                        primary_index, primary_value = all_dofs[0]
                        for index, value in all_dofs[1:]:
                            if value == primary_value:  # Enforce compatibility if values match
                                row = np.zeros(num_dof)
                                row[primary_index] = 1  # Positive sign for the primary DOF
                                row[index] = -1  # Negative sign for the compared DOF
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


#%% classify u_unique and q nodes in terms of free, fixed, and prescribed

def classify_u_q(nodes, unique_dofs, dof_indices):
    """
    Classifies unique DOFs as free, fixed, or prescribed and extracts their indices.
    
    Args:
    nodes (list of Node): List of all nodes.
    unique_dof_indices (list): List of indices considered unique.
    dof_indices (dict): Maps (node_id, element_id) to a dict of DOFs and their indices.
    
    Returns:
    
    q_classified: dict: A dictionary with keys 'free', 'fixed', 'prescribed', and 'values',
          containing the indices in `unique_dof_indices` and values for prescribed DOFs in the original u space
    
    U_classified: dict: A dictionary with keys 'free', 'fixed', 'prescribed', and 'values',
      containing the indices in `unique_dof_indices` and values for prescribed DOFs in q_space          
    """
    reverse_dof_lookup = {index: (node_id, dof_name) for (node_id, element_id), dofs in dof_indices.items() for dof_name, index in dofs.items()}
    
    q_classified = {
        'free': [],
        'fixed': [],
        'prescribed': [],
        'values': []
    }
    
    U_classified = {
    'free': [],
    'fixed': [],
    'prescribed': [],
    'values': []
    }

    for index, dof_index in enumerate(unique_dofs):
        if dof_index not in reverse_dof_lookup:
            continue  # Skip if no matching node-dof pair is found
        node_id, dof_name = reverse_dof_lookup[dof_index]
        node = next((n for n in nodes if n.id == node_id), None)
        
        if not node:
            continue  # Skip if no node is found
        
        dof_value = node.dofs[dof_name]
        if dof_value is None:
            q_classified['free'].append(index)
            U_classified['free'].append(dof_index)
        elif dof_value == 0:
            q_classified['fixed'].append(index)
            U_classified['fixed'].append(dof_index)
        else:
            q_classified['prescribed'].append(index)
            q_classified['values'].append((index, dof_value))
            U_classified['prescribed'].append(dof_index)
            U_classified['values'].append(dof_index)

    return q_classified, U_classified

q_classified, U_classified = classify_u_q(nodes, unique_dofs, dof_indices)

print("U Free DOF Indices:", U_classified['free'])
print("U Fixed DOF Indices:", U_classified['fixed'])
print("U Prescribed DOF Indices:", U_classified['prescribed'])
print("U Values of Prescribed DOFs:", U_classified['values'])

print("q Free DOF Indices:", q_classified['free'])
print("q Fixed DOF Indices:", q_classified['fixed'])
print("q Prescribed DOF Indices:", q_classified['prescribed'])
print("q Values of Prescribed DOFs:", q_classified['values'])