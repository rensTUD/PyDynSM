# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:52:01 2024

@author: rensv
"""

# Import dependencies
import numpy as np
from scipy.linalg import inv
from ...elements import ElementFactory
from collections import defaultdict

# Class definition
class Analysis:
    """
    A class to handle structural analysis computations such as stiffness matrices, force vectors, and displacement calculations.
    """

    def GlobalStiffness(self, nodes, elements, omega):
        """
        Assembles the global unconstrained stiffness matrix.

        Parameters
        ----------
        nodes : list of Node
            List of nodes in the structural system.
        elements : list of Element
            List of elements in the structural system.
        omega : float
            Frequency parameter.

        Returns
        -------
        k_global : numpy.ndarray
            Global stiffness matrix.
        """
        # load k_global
        k_global = np.zeros((self.num_dofs, self.num_dofs), complex)
        
        # set up block matrix per element
        for e in elements:
            elmat = e.Stiffness(omega)
            idofs = e.get_element_node_dof_indices_global()
            k_global[np.ix_(idofs, idofs)] = elmat
            
        # sort by unique / redundant dofs to comply with the structure of L
        k_uu = k_global[np.ix_(self.unique_dofs,self.unique_dofs)]
        k_ur = k_global[np.ix_(self.unique_dofs,self.redundant_dofs)]
        k_ru = k_global[np.ix_(self.redundant_dofs,self.unique_dofs)]
        k_rr = k_global[np.ix_(self.redundant_dofs,self.redundant_dofs)]
        
        k_global = np.block([[k_uu, k_ur],
                                    [k_ru, k_rr]])
        
        # constrain by using L
        K_global = self.L.T @ k_global @ self.L
    
        return K_global

    def GlobalForce(self, nodes, elements, omega):
        """
        Assembles the global unconstrained force vector.

        Parameters
        ----------
        nodes : list of Node
            List of nodes in the structural system.
        elements : list of Element
            List of elements in the structural system.
        omega : float
            Frequency parameter.

        Returns
        -------
        f_global : numpy.ndarray
            Global force vector.
        """
        
        # set the global force vector
        f_global = np.zeros(self.num_dofs, complex)
        
        # get all element loads
        for element in elements:
            # if no loads present continue loop
            if not element.element_loads:
                continue
            
            # get dofs of the element
            dofs = element.get_element_node_dof_indices_global()
            
            # for element_load in element.element_loads:
            f_global[np.ix_(dofs)] += element.EvaluateDistributedLoad(element.element_loads, omega)
        
        # get all nodal loads - ASSUMING EVERYTHING IS ALREADY WITHIN GLOBAL DIRECTIONS 
        for node in nodes:
            if not node.nodal_loads:
                continue
            
            # load the dof container
            node_dof_container = node.dof_container
            
            # loop over nodal loads
            for dof_name, load in node.nodal_loads.items():
                # add the load to the dof if dof is present in the current dof config of that node (probably should add this check in the nodal class..)
                if node_dof_container.has_dof(dof_name):
                    dof = node_dof_container.dofs[dof_name]
                    index = dof.index
                    if index is not None:
                        value = load(omega) if callable(load) else load
                        f_global[index] += value
        
        # sort by unique / redundant dofs to comply with the structure of L
        f_global_u = f_global[np.ix_(self.unique_dofs)]
        f_global_r = f_global[np.ix_(self.redundant_dofs)]
        
        f_global = np.hstack((f_global_u,f_global_r))
        
        # constrain by using L
        F_global = self.L.T @ f_global 
        
        return F_global

    def GlobalConstrainedStiffness(self, nodes, elements, omega):
        """
        Constrains the global stiffness matrix.

        Parameters
        ----------
        nodes : list of Node
            List of nodes in the structural system.
        elements : list of Element
            List of elements in the structural system.
        omega : float
            Frequency parameter.

        Returns
        -------
        k_constrained : numpy.ndarray
            Constrained global stiffness matrix.
        """
        k_global = self.GlobalStiffness(nodes, elements, omega)
        return k_global[np.ix_(self.free_dofs, self.free_dofs)]

    def GlobalConstrainedForce(self, nodes, elements, omega):
        """
        Constrains the global force vector.

        Parameters
        ----------
        nodes : list of Node
            List of nodes in the structural system.
        elements : list of Element
            List of elements in the structural system.
        omega : float
            Frequency parameter.

        Returns
        -------
        f_constrained : numpy.ndarray
            Constrained global force vector.
        """
        k_global = self.GlobalStiffness(nodes, elements, omega)
        f_global = self.GlobalForce(nodes, elements, omega)
        
        constrained_indices = list(self.constrained_dofs.keys())
        constrained_values = list(self.constrained_dofs.values())
        
        K_free_constrained = k_global[np.ix_(self.free_dofs, constrained_indices)]
        F_free = f_global[self.free_dofs]
        
        return F_free - K_free_constrained @ constrained_values

    def SolveUfree(self, Kc_global, fc_global):
        """
        Solves the free displacements.

        Parameters
        ----------
        Kc_global : numpy.ndarray
            Constrained global stiffness matrix.
        fc_global : numpy.ndarray
            Constrained global force vector.

        Returns
        -------
        u_free : numpy.ndarray
            Free displacements.
        """
        return np.linalg.inv(Kc_global) @ fc_global

    def SupportReactions(self, k_global, u_free, f_global):
        """
        Calculates the support reactions.

        Parameters
        ----------
        k_global : numpy.ndarray
            Global stiffness matrix.
        u_free : numpy.ndarray
            Free displacements.
        f_global : numpy.ndarray
            Global force vector.

        Returns
        -------
        reactions : numpy.ndarray
            Support reactions.
        """
        constrained_indices = list(self.constrained_dofs.keys())
        constrained_values = list(self.constrained_dofs.values())
        
        Kcf = k_global[np.ix_(constrained_indices, self.free_dofs)]
        Kcc = k_global[np.ix_(constrained_indices, constrained_indices)]
        
        return (Kcf @ u_free) + (Kcc @ constrained_values) - f_global[np.ix_(constrained_indices)]

    def FullDisplacement(self, u_free):
        """
        Returns the full displacement vector for the entire structure.

        Parameters
        ----------
        u_free : numpy.ndarray
            Free displacements.

        Returns
        -------
        u_full : numpy.ndarray
            Full displacement vector.
        """
        constrained_indices = list(self.constrained_dofs.keys())
        constrained_values = list(self.constrained_dofs.values())
        
        u_full = np.zeros(len(self.free_dofs) + len(constrained_indices), dtype=complex)
        
        u_full[np.ix_(self.free_dofs)] = u_free
        u_full[np.ix_(constrained_indices)] = constrained_values
        
        return u_full
    
    def ElementDisplacements(self, elements, u_nodes_global, omega, num_points=20):
        """
        Gathers the displacements for each element from the element-level calculations.
    
        Parameters
        ----------
        elements : list of Element
            List of elements in the structural system.
        u_nodes_global : numpy.ndarray
            Global displacements of nodes.
        omega : float
            Frequency parameter.
        num_points : int, optional
            Number of points along the element to calculate displacements, default is 20.
    
        Returns
        -------
        element_displacements : dict
            Dictionary with element IDs as keys and their global displacements (as numpy arrays) as values.
        """
        
        # Apply L matrix to get displacements of redundant nodes as well
        u_nodes_global_all = self.L @ u_nodes_global
        # Get the current order
        current_order = self.unique_dofs + self.redundant_dofs
        
        # Create an index map to reorder the elements to match sorted_order
        sorted_order = sorted(current_order)
        index_map = [current_order.index(idx) for idx in sorted_order]
        
        # Reorder u_nodes_global_all using the index map
        u_nodes_global_all_sorted = np.array(u_nodes_global_all)[index_map]
        
        element_displacements = {}
    
        for element in elements:
            # Retrieve global nodal DOF indices for the current element
            global_dof_indices = element.get_node_dof_indices_global()
    
            # Extract the relevant global displacements for this element
            u_element_global = u_nodes_global_all_sorted[global_dof_indices]
    
            # Calculate the displacements for the element using its Displacements method
            u_elem = element.Displacements(u_element_global, omega, num_points)
    
            # Store the displacements in the dictionary using the element's ID
            element_displacements[element.id] = u_elem
    
        return element_displacements
    
    def SolveEigenvector(self, nodes, elements, omega, fixed_index=0):
            """
            Computes an eigenvector for the system using the constrained global stiffness matrix.
    
            Parameters
            ----------
            nodes : list of Node
                List of nodes in the structural system.
            elements : list of Element
                List of elements in the structural system.
            omega : float
                Frequency parameter.
            fixed_index : int, optional
                Index of the eigenvector component to be set to 1 (default is 0).
    
            Returns
            -------
            U : numpy.ndarray
                Computed eigenvector normalized by fixing one component to 1.
            """
            # Retrieve the constrained global stiffness matrix
            K_constrained = self.GlobalConstrainedStiffness(nodes, elements, omega)
    
            # Ensure the matrix is not singular
            if np.linalg.matrix_rank(K_constrained) < K_constrained.shape[0]:
                raise ValueError("Singular stiffness matrix, the system may be improperly constrained.")
    
            # Define the number of DOFs
            num_dofs = K_constrained.shape[0]
    
            # Partition the matrix explicitly
            free_dofs = [i for i in range(num_dofs) if i != fixed_index]
    
            # Extract submatrices
            K_rr = K_constrained[np.ix_(free_dofs, free_dofs)]  # Reduced stiffness matrix
            K_rp = K_constrained[np.ix_(free_dofs, [fixed_index])]  # Column corresponding to fixed index
    
            # Solve for the unknown components
            try:
                U_r = -np.linalg.solve(K_rr, K_rp).flatten()  # Solve for U_r
            except np.linalg.LinAlgError:
                raise ValueError("Numerical issue: unable to solve for eigenvector.")
                
            # Reconstruct the full eigenvector
            U = np.zeros(num_dofs, dtype=complex)
            U[fixed_index] = 1  # Set fixed value
            U[free_dofs] = U_r  # Assign solved components
    
            return U                
    
    def ElementForces(self, elements, u_nodes_global, omega, num_points=20):
        """
        Gathers the displacements for each element from the element-level calculations.
    
        Parameters
        ----------
        elements : list of Element
            List of elements in the structural system.
        u_nodes_global : numpy.ndarray
            Global displacements of nodes.
        omega : float
            Frequency parameter.
        num_points : int, optional
            Number of points along the element to calculate displacements, default is 20.
    
        Returns
        -------
        element_forces : dict
            Dictionary with element IDs as keys and their global displacements (as numpy arrays) as values.
        """
        
        # Apply L matrix to get displacements of redundant nodes as well
        u_nodes_global_all = self.L @ u_nodes_global
        # Get the current order
        current_order = self.unique_dofs + self.redundant_dofs
        
        # Create an index map to reorder the elements to match sorted_order
        sorted_order = sorted(current_order)
        index_map = [current_order.index(idx) for idx in sorted_order]
        
        # Reorder u_nodes_global_all using the index map
        u_nodes_global_all_sorted = np.array(u_nodes_global_all)[index_map]
        
        element_forces = {}
    
        for element in elements:
            # Retrieve global nodal DOF indices for the current element
            global_dof_indices = element.get_node_dof_indices_global()
    
            # Extract the relevant global displacements for this element
            u_element_global = u_nodes_global_all_sorted[global_dof_indices]
    
            # Calculate the displacements for the element using its Displacements method
            f_elem = element.Forces(u_element_global, omega, num_points)
    
            # Store the displacements in the dictionary using the element's ID
            element_forces[element.id] = f_elem
    
        return element_forces
    

    

    def ElementStresses(self, elements, u_nodes_global, omega, num_points=20):
        """
        Gathers the stresses for each element from the element-level calculations.
    
        Parameters
        ----------
        elements : list of Element
            List of elements in the structural system.
        u_nodes_global : numpy.ndarray
            Global displacements of nodes.
        omega : float
            Frequency parameter.
        num_points : int, optional
            Number of points along the element to calculate displacements, default is 20.
    
        Returns
        -------
        element_stresses : dict
            Dictionary with element IDs as keys and their global displacements (as numpy arrays) as values.
        """
        
        # Apply L matrix to get displacements of redundant nodes as well
        u_nodes_global_all = self.L @ u_nodes_global
        # Get the current order
        current_order = self.unique_dofs + self.redundant_dofs
        
        # Create an index map to reorder the elements to match sorted_order
        sorted_order = sorted(current_order)
        index_map = [current_order.index(idx) for idx in sorted_order]
        
        # Reorder u_nodes_global_all using the index map
        u_nodes_global_all_sorted = np.array(u_nodes_global_all)[index_map]
        
        element_stresses = {}
    
        for element in elements:
            # Retrieve global nodal DOF indices for the current element
            global_dof_indices = element.get_node_dof_indices_global()
    
            # Extract the relevant global displacements for this element
            u_element_global = u_nodes_global_all_sorted[global_dof_indices]
    
            # Calculate the displacements for the element using its Displacements method
            s_elem = element.Stresses(u_element_global, omega, num_points)
    
            # Store the displacements in the dictionary using the element's ID
            element_stresses[element.id] = s_elem
    
        return element_stresses
# %% specific methods for this method of analysis

    def find_unique_redundant_dofs(self, B):
        """
        Identifies unique and redundant degrees of freedom (DOFs) based on the constraint matrix B.

        Parameters
        ----------
        B : numpy.ndarray
            The constraint matrix where rows represent constraints and columns represent DOFs.

        Returns
        -------
        unique_dofs : list of int
            Sorted list of unique DOFs.
        redundant_dofs : list of int
            Sorted list of redundant DOFs.
        """
        self.num_dofs = B.shape[1]
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
        zero_dofs = {i for i in range(self.num_dofs ) if np.all(B[:, i] == 0)}
        unique_dofs = (unique_dofs | zero_dofs) - redundant_dofs

        return sorted(unique_dofs), sorted(redundant_dofs)

    def calculate_L(self, B, unique_dofs, redundant_dofs):
        """
        Calculates the transformation matrix L based on the constraint matrix B.

        Parameters
        ----------
        B : numpy.ndarray
            The constraint matrix where rows represent constraints and columns represent DOFs.
        unique_dofs : list of int
            List of unique DOFs.
        redundant_dofs : list of int
            List of redundant DOFs.

        Returns
        -------
        L : numpy.ndarray
            The transformation matrix L.

        Raises
        ------
        ValueError
            If the matrix B_r is singular and cannot be inverted.
        """
        # check if there are any redundant dofs
        if len(redundant_dofs) != 0:
            B_r = B[:, redundant_dofs]
            B_u = B[:, unique_dofs]
    
            try:
                B_r_inv = inv(B_r)
                L_lower = -B_r_inv @ B_u
                L = np.vstack((np.eye(len(unique_dofs)), L_lower))
                return L
            
            except np.linalg.LinAlgError:
                raise ValueError("Matrix B_r is singular and cannot be inverted.")
        # if not, return L solely based on unique dofs
        else:
            L = np.eye(len(unique_dofs))
            return L

    def assign_dof_indices_old(self, nodes, elements):
        """
        Assigns global indices to the degrees of freedom (DOFs) of all elements.
        The element dofs are leading here, nodal dofs will be used to connect element dofs

        Parameters
        ----------
        elements : list of Element
            List of elements in the structural system.

        Returns
        -------
        dof_indices : dict
            Dictionary mapping (node_id, element_id) to a dict of DOFs and their global indices.
        """
        # TODO - currently we loop first over nodes which felt that it makes sense. However the global matrix before constraining by L is broken up. There are no block matrices per element in that case
        dof_indices = defaultdict(dict)
        global_index = 0  # Start a global index counter for all DOFs in the system

        for node in nodes:
            for element in node.connected_elements:
                
                # use the DOF configuration of the element for this specific node
                for dof in element.dofs[node.id].keys():
                    dof_indices[(node.id, element.id)][dof] = global_index
                    # also assign to the element itself as this can be useful
                    element.dof_indices[node.id][dof] = global_index
                    global_index += 1  # Increment global index for each DOF

        return dof_indices
    
    def assign_dof_indices(self, nodes, elements):
        """
        Assigns global indices to the degrees of freedom (DOFs) of all elements.
        The element DOFs are leading here; nodal DOFs will be used to connect element DOFs.
    
        Parameters
        ----------
        nodes : list of Node
            List of nodes in the structural system.
        elements : list of Element
            List of elements in the structural system.
    
        Returns
        -------
        dof_indices : defaultdict
            Dictionary mapping (node_id, element_id) to a dict of DOFs and their global indices.
        """
        dof_indices = defaultdict(dict)
        global_index = 0  # Start a global index counter for all DOFs in the system
        
        for element in elements:
    
            for node in element.nodes:
                node_id = node.id
                element_dof_container = element.dof_containers[node.id]
                node_dof_container = node.dof_container
                
                # Use the DOF configuration of the element for this specific node
                for dof_name, dof in element_dof_container.dofs.items():
                    dof_indices[(node_id, element.id)][dof_name] = global_index
                    
                    # Assign the index also to the dof in the global dof container of the element
                    dof.index = global_index
                    
                    # TODO - SEE HOW TO HANDLE DOF INDICES FOR ELEMENT DOFS AND NODE DOFS IN A NEAT WAY. NOW WE HAVE THREE DIFFERENT DICTIONARIES STORING THEM. NOT THAT EFFICIENT ATM
                    # same for the node, but only if present in node
                    if node_dof_container.has_dof(dof_name):
                        # add index only if it hasn't been assigned yet - # TODO this may not be that useful..
                        if not node_dof_container.dofs[dof_name].index:
                            node_dof_container.dofs[dof_name].index = global_index
                        
                    global_index += 1  # Increment global index for each DOF
        
        num_dof = global_index                    
        return dof_indices, num_dof


    def build_matrix_B(self, nodes, elements, dof_indices, num_dof):
        """
        Builds the constraint matrix B for the structural system.
    
        Parameters
        ----------
        nodes : list of Node
            List of nodes in the structural system.
        elements : list of Element
            List of elements in the structural system.
        dof_indices : dict
            Dictionary mapping (node_id, element_id) to a dict of DOFs and their global indices.
    
        Returns
        -------
        B : numpy.ndarray
            The constraint matrix B.
        """
        # num_dof = sum(len(dofs) for element in elements for dofs in element.dof_container.values())  # Calculate total DOFs based on the elements global dofs at each connection
    
        # Preallocate a larger matrix to avoid expensive np.vstack() operations
        B = np.zeros((num_dof, num_dof))
        constraint_counter = 0  # Keeps track of which row in B we are filling
    
        for node in nodes:
            if len(node.connected_elements) > 1:  # Only consider nodes connected to multiple elements
                # find the element classes and the indices of its dofs that are connected to the current node
                connected_indices = [(element, dof_indices[(node.id, element.id)])  
                                     for element in node.connected_elements
                                     if (node.id, element.id) in dof_indices]
                
                # only continue if connected indices are actually found
                if connected_indices:
                    # 
                    for dof_name, nodal_dof in node.dof_container.dofs.items():
                        # Collect DOFs from elements that are monolithically connected
                        all_dofs = [
                            indices[dof_name] for element, indices in connected_indices
                            if dof_name in indices and element.dof_containers[node.id].dofs[dof_name].value == nodal_dof.value
                            and element.constraint_types[node.id][dof_name] == "monolithic"
                        ]
                        # continue only if we still have connected dofs
                        if len(all_dofs) > 1:
                            # Use the first DOF as the reference for compatibility
                            primary_index = all_dofs[0]
                            for index in all_dofs[1:]:
                                row = np.zeros(num_dof)
                                row[primary_index] = 1  # Primary DOF
                                row[index] = -1  # Enforce compatibility with other DOFs
                                B[constraint_counter, :] = row  # Fill row directly instead of stacking
                                constraint_counter += 1
    
        # Return the trimmed matrix if we overestimated the number of constraints
        return B[:constraint_counter, :]


    def connectivity(self, nodes, elements):
        """
        Computes the connectivity information for the structural system.

        Parameters
        ----------
        nodes : list of Node
            List of nodes in the structural system.
        elements : list of Element
            List of elements in the structural system.

        Returns
        -------
        B : numpy.ndarray
            The constraint matrix B.
        L : numpy.ndarray
            The transformation matrix L.
        dof_indices : dict
            Dictionary mapping (node_id, element_id) to a dict of DOFs and their global indices.
        unique_dofs : list of int
            Sorted list of unique DOFs.
        redundant_dofs : list of int
            Sorted list of redundant DOFs.
        free_dofs: list
            List of indices in `unique_dofs` that are free.
        constrained_dofs: list
            List of lists where each sublist contains a DOF index and its value.
        """
        
        self.dof_indices, num_dof = self.assign_dof_indices(nodes, elements)
        
        self.B = self.build_matrix_B(nodes, elements, self.dof_indices,num_dof)
        
        self.unique_dofs, self.redundant_dofs = self.find_unique_redundant_dofs(self.B)
        
        self.L = self.calculate_L(self.B, self.unique_dofs, self.redundant_dofs)
        
        self.free_dofs, self.constrained_dofs = self.classify_free_constrained_dofs(nodes, elements, np.array(self.unique_dofs), self.dof_indices)

    
    def classify_free_constrained_dofs_old(self, nodes, elements, unique_dofs, dof_indices):
        """
        Classifies unique DOFs as free or constrained (fixed or prescribed) and extracts their indices.
        
        Parameters
        --------
        nodes: list of Node
            List of all nodes.
        elements: list of Element
            List of elements in the structural system.
        unique_dofs: list
            List of indices considered unique.
        dof_indices: dict
            Maps (node_id, element_id) to a dict of DOFs and their indices.
        
        Returns
        -------
        free_dofs: list
            List of indices in `unique_dofs` that are free.
        constrained_dofs: list
            List of lists where each sublist contains a DOF index and its value.
        """
        
        # create dict to look up via O(1) instead of O(n) for larger system        
        unique_dof_dict = {old_index: idx for idx, old_index in enumerate(unique_dofs)}

        # introduce empty lists
        free_dofs = []
        constrained_dofs = {}
        
        # Iterate over all nodes and their connected elements to classify DOFs
        for node in nodes:
            
            # loop over all connected elements
            for element in node.connected_elements:
                # get element dofs
                element_dofs = element.dofs[node.id]
                # get dofs from dof_indices dictionary
                for dof, old_index in dof_indices[(node.id,element.id)].items():
                    
                    # check whether the current dof is in the unique dofs
                    new_index = unique_dof_dict.get(old_index)
                    
                    # only continue if only 1 location is found
                    if new_index is not None:
                        
                        # get value from element
                        dof_value = element_dofs[dof]
                        
                        # set accordingly, i.e. if None its a free node otherwise its prescribed
                        if dof_value is None:
                            free_dofs.append(new_index)
                        else:
                            constrained_dofs[new_index] = dof_value
        
         # TODO - return sorted list, not per se needed I guess. Furthermore I could change the constrained dofs list to a dictionary
        # return sorted(free_dofs), sorted(constrained_dofs, key=lambda x: x[0])
        return np.array(free_dofs), constrained_dofs
    
    def classify_free_constrained_dofs(self, nodes, elements, unique_dofs, dof_indices):
        """
        Classifies unique DOFs as free or constrained (fixed or prescribed) and extracts their indices.
    
        Parameters
        ----------
        nodes : list of Node
            List of all nodes.
        elements : list of Element
            List of elements in the structural system.
        unique_dofs : list
            List of indices considered unique.
        dof_indices : dict
            Maps (node_id, element_id) to a dict of DOFs and their indices.
    
        Returns
        -------
        free_dofs : numpy.ndarray
            Array of indices in `unique_dofs` that are free.
        constrained_dofs : dict
            Dictionary mapping DOF indices to their prescribed values.
        """
        # Create a mapping from old DOF indices to new indices in unique_dofs
        unique_dof_dict = {old_index: idx for idx, old_index in enumerate(unique_dofs)}
    
        free_dofs = []
        constrained_dofs = {}
    
        # Iterate over all elements and their nodes to classify DOFs
        for element in elements:
            for node in element.nodes:
                node_id = node.id
                element_id = element.id
                element_dof_container = element.dof_containers[node_id]
    
                # Get DOFs from dof_indices dictionary
                for dof_name, old_index in dof_indices[(node_id, element_id)].items():
                    # Check whether the current DOF is in the unique DOFs
                    new_index = unique_dof_dict.get(old_index)
    
                    if new_index is not None:
                        # Get value from element's DOFContainer
                        dof_value = element_dof_container.get_dof(dof_name).value
    
                        # Classify DOF as free or constrained
                        if dof_value is None:
                            free_dofs.append(new_index)
                        else:
                            constrained_dofs[new_index] = dof_value
    
        return np.array(free_dofs), constrained_dofs