# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:52:01 2024

@author: rensv
"""

# Import dependencies
import numpy as np
from scipy.linalg import inv
from ...elements import ElementFactory

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
        Ndofs = 3  # Assuming 2D with [x, z, phi] = 3 DOFs per node
        k_global = np.zeros((Ndofs * len(nodes), Ndofs * len(nodes)), complex)
        
        for e in elements:
            elmat = e.Stiffness(omega)
            idofs = e.GlobalDofs()
            k_global[np.ix_(idofs, idofs)] += elmat
    
        return k_global

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
        Ndofs = 3  # Assuming 2D with [x, z, phi] = 3 DOFs per node
        f_global = np.zeros(Ndofs * len(nodes), complex)
        
        for element in elements:
            if not element.element_nodal_loads:
                continue
            
            left_dofs = element.nodes[0].dofs
            right_dofs = element.nodes[1].dofs
            dofs = np.hstack([left_dofs, right_dofs])
            
            for element_nodal_load in element.element_nodal_loads:
                f_global[np.ix_(dofs)] += element.EvaluateDistributedLoad(element_nodal_load, omega)
        
        for n in nodes:
            if not n.nodal_forces:
                continue
            
            for nodal_force in n.nodal_forces:
                force_components = np.array([(force(omega) if callable(force) else force) for force in nodal_force])
                f_global[n.dofs] += force_components
            
        return f_global

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
        
        constrained_indices = [dof[0] for dof in self.constrained_dofs]
        constrained_values = [dof[1] for dof in self.constrained_dofs]
        
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
        constrained_indices = [dof[0] for dof in self.constrained_dofs]
        constrained_values = [dof[1] for dof in self.constrained_dofs]
        
        Kcf = k_global[np.ix_(constrained_indices, self.free_dofs)]
        Kcc = k_global[np.ix_(constrained_indices, constrained_indices)]
        
        return (Kcf @ u_free) + (Kcc @ constrained_values) - f_global[constrained_indices]

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
        constrained_indices = [dof[0] for dof in self.constrained_dofs]
        constrained_values = [dof[1] for dof in self.constrained_dofs]
        
        u_full = np.zeros(len(self.free_dofs) + len(constrained_indices), dtype=complex)
        
        u_full[self.free_dofs] = u_free
        u_full[constrained_indices] = constrained_values
        
        return u_full

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
        B_r = B[:, redundant_dofs]
        B_u = B[:, unique_dofs]

        try:
            B_r_inv = inv(B_r)
            L_lower = -B_r_inv @ B_u
            L = np.vstack((np.eye(len(unique_dofs)), L_lower))
            return L
        except np.linalg.LinAlgError:
            raise ValueError("Matrix B_r is singular and cannot be inverted.")

    def assign_dof_indices(self, nodes, elements):
        """
        Assigns global indices to the degrees of freedom (DOFs) of all elements.

        Parameters
        ----------
        elements : list of Element
            List of elements in the structural system.

        Returns
        -------
        dof_indices : dict
            Dictionary mapping (node_id, element_id) to a dict of DOFs and their global indices.
        """
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

    def build_matrix_B(self, nodes, elements, dof_indices):
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
                    # TODO - CHANGE SUCH THAT THE ASSUMPTION IS NOT TAKEN AND WE CHECK PER ELEMENT AT THE SPECIFIC NODE
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
        
        dof_indices = self.assign_dof_indices(nodes, elements)
        
        self.B = self.build_matrix_B(nodes, elements, dof_indices)
        
        unique_dofs, redundant_dofs = self.find_unique_redundant_dofs(self.B)
        
        self.L = self.calculate_L(self.B, unique_dofs, redundant_dofs)
        
        self.free_dofs, self.constrained_dofs = self.classify_free_constrained_dofs(nodes, unique_dofs, dof_indices)
        
        

    def classify_free_constrained_dofs(self, nodes, unique_dofs, dof_indices):
        """
        Classifies unique DOFs as free or constrained (fixed or prescribed) and extracts their indices.
        
        Parameters
        --------
        nodes: list of Node
            List of all nodes.
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
        reverse_dof_lookup = {index: (node_id, dof_name) for (node_id, element_id), dofs in dof_indices.items() for dof_name, index in dofs.items()}
        
        free_dofs = []
        constrained_dofs = []

        for index, dof_index in enumerate(unique_dofs):
            if dof_index not in reverse_dof_lookup:
                continue  # Skip if no matching node-dof pair is found
            node_id, dof_name = reverse_dof_lookup[dof_index]
            node = next((n for n in nodes if n.id == node_id), None)
            
            if not node:
                continue  # Skip if no node is found
            
            dof_value = node.dofs[dof_name]
            if dof_value is None:
                free_dofs.append(index)
            else:
                constrained_dofs.append([index, dof_value])

        return free_dofs, constrained_dofs