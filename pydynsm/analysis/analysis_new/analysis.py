# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:52:01 2024

@author: rensv
"""

import numpy as np

class Analysis:
    """
    A class to handle structural analysis computations such as stiffness matrices, force vectors, and displacement calculations.
    """

    def __init__(self):
        self.constrained_dofs = []
        self.free_dofs = []

    def UpdateDofsDecorator(func):
        """Decorator to update constrained DOFs before method call."""
        def wrapper(self, nodes, *args, **kwargs):
            self._UpdateConstrainedDofs(nodes)  # Update DOFs before the function call
            return func(self, nodes, *args, **kwargs)
        return wrapper

    def _UpdateConstrainedDofs(self, nodes):
        """
        Update the list of constrained DOFs based on registered nodes.
        """
        self.constrained_dofs.clear()
        for node in nodes:
            self.constrained_dofs.extend(node.constrained_dofs)
        self._UpdateFreeDofs(nodes)

    def _UpdateFreeDofs(self, nodes):
        """
        Update the list of free DOFs.
        """
        all_dofs = set(range(3 * len(nodes)))  # Assuming 2D with [x, z, phi] = 3 DOFs per node
        constrained_indices = {dof[0] for dof in self.constrained_dofs}
        self.free_dofs[:] = list(all_dofs - constrained_indices)

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

    @UpdateDofsDecorator
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

    @UpdateDofsDecorator
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
