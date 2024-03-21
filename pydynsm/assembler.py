# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 19:21:45 2024

@author: rensv
"""

# %% import dependencies
import numpy as np
import matplotlib.pyplot as plt

# %% import others

from .import plotter

# %% main class definition

class Assembler:
    
    # number of dofs, for now we have 2D and thus: [x, z, phi] = 3
    
    Ndofs = 3
    
    def __init__(self, name, plotter=plotter):
        '''
        Initialisation of assembler.
        '''
        
        self.name = name  # Assembler / project name
        
        self.elements = []  # Empty elements list
        
        self.nodes = []  # Empty nodes list
        
        self.constrained_dofs = [] # Stores lists of [DOF index, value]
        self.free_dofs = [] # stores the free dofs
        
        
        # load dependencies that are injected
        self.StructurePlotter = plotter.StructurePlotter()
        
        print(f"Assembler '{self.name}' successfully initialised")
        
    def RegisterNode(self, node):
        '''
        Adds a node to the assembler if it's not already registered.
        '''
        if node not in self.nodes:
            self.nodes.append(node)         # add node to list
            self._UpdateConstrainedDofs()    # update lists of constrained and free dofs

    def RegisterElement(self, element):
        '''
        Adds an element to the assembler if it's not already registered.
        '''
        if element not in self.elements:
            self.elements.append(element)
    
    def PlotStructure(self, plot_elements=False):
        
        self.StructurePlotter.PlotNodes(self.nodes)
        
        if plot_elements:
            self.StructurePlotter.PlotElements(self.elements)
            
        self.StructurePlotter.ShowStructure(f'Structure Configuration: {self.name}')

    def GlobalStiffness(self, omega):
        '''
        Assembles the global unconstrained stiffness matrix
        
        # TODO: should be replaced with the matrix multiplication as in the paper Thanasis found!
        '''
        # introduce global stiffness matrix with size [3*Nnodes, 3*Nnodes] (as we have: x,z,phi)
        k_global = np.zeros( (Assembler.Ndofs*len(self.nodes),Assembler.Ndofs*len(self.nodes) ), complex)
        
        for e in self.elements:
            # get global stiffness matrix of the element and its global dofs
            elmat = e.Stiffness(omega)
            idofs = e.GlobalDofs()
            
            # assign to the global stiffness matrix of the structure
            k_global[np.ix_(idofs,idofs)] += elmat
    
        return k_global            

    def GlobalForce(self, omega):
        '''
        Assembles the global unconstrained force vector
        
        # TODO: should be replaced with the matrix multiplication as in the paper Thanasis found!
        '''
        f_global = np.zeros(Assembler.Ndofs*len(self.nodes),complex)
        
        for n in self.nodes:
            f_global[n.dofs] += n.p
            
        return f_global
    
    def GlobalConstrainedStiffness(self, omega):
        '''
        Constrains the global stiffness matrix with the use of static condensation (I think?)
        
        # TODO: check what exactly happens here and explain it in docstring
        '''
        
        # TODO: Right now we update every time we need it, should make either a decorator or an observer function to dynamically handle this and not have to think of it anymore
        # update constrained and free dofs list
        self._UpdateConstrainedDofs()    
        
        return self.GlobalStiffness(omega)[np.ix_(self.free_dofs,self.free_dofs)]
    
    def GlobalConstrainedForce(self, omega):
        '''
        Constrains the global stiffness matrix with the use of static condensation (I think?)
        
        # TODO: check what exactly happens here and explain it in docstring
        '''
        
        # update constrained and free dofs list
        self._UpdateConstrainedDofs()    
        
        # extract the constrained dofs and their values
        constrained_dofs = [dof[0] for dof in self.constrained_dofs] 
        constrained_values = [value[1] for value in self.constrained_dofs] 
        
        # global free,constrained stiffness matrix
        K_free_constrained = self.GlobalStiffness(omega)[np.ix_(self.free_dofs,constrained_dofs)]
        
        # global free forcing vector
        F_free = self.GlobalForce(omega)[self.free_dofs]
        
        return F_free - K_free_constrained @ constrained_values
    
    def SolveUfree(self, Kc_global, fc_global):
        '''
        Solves the free displacements
        '''

        return np.linalg.inv(Kc_global) @ fc_global
    
    def SupportReactions(self, k_global, u_free, f_global):
        '''
        Gets the support reactions
        '''
        # update constrained and free dofs list
        self._UpdateConstrainedDofs()    
        
        # extract the constrained dofs and their values
        constrained_dofs = [dof[0] for dof in self.constrained_dofs] 
        constrained_values = [value[1] for value in self.constrained_dofs] 
        
        Kcf = k_global[np.ix_(constrained_dofs,self.free_dofs)]
        Kcc = k_global[np.ix_(constrained_dofs,constrained_dofs)]
        
        return (Kcf @ u_free) + (Kcc @ constrained_values) - f_global[constrained_dofs]

    def FullDisplacement(self, u_free):
        '''
        function that returns the full displacement of the whole structure (i.e. calculated free dofs and constrained ones)
        '''
        
        # extract the constrained dofs and their values
        constrained_dofs = [dof[0] for dof in self.constrained_dofs] 
        constrained_values = [value[1] for value in self.constrained_dofs] 
        
        # initialise full displacement vector
        u_full = np.zeros(len(self.free_dofs) + len(constrained_dofs), dtype=complex)
        
        # populate        
        u_full[self.free_dofs] = u_free
        u_full[constrained_dofs] = constrained_values
        
        return u_full
    
    def _UpdateConstrainedDofs(self):
        """
        Update the list of constrained DOFs based on registered nodes.
        """
        self.constrained_dofs.clear()
        for node in self.nodes:
            self.constrained_dofs.extend(node.constrained_dofs)  # Assuming node.constrained_dofs is accessible
        self._UpdateFreeDofs()

    def _UpdateFreeDofs(self):
        """
        Update the list of free DOFs.
        """
        all_dofs = set(range(Assembler.Ndofs * len(self.nodes)))  
        constrained_indices = {dof[0] for dof in self.constrained_dofs}
        self.free_dofs = list(all_dofs - constrained_indices)        
        
    # TODO
    def SaveAssembler(self):
       '''
       Function to be built that can save the current structure you are working on
       '''
       pass
    
    # TODO
    @classmethod
    def LoadAssembler(cls):
        '''
        Class method that will load a specific saved assembler / project to quickly continue with it
        '''
        pass
    
        
            
        
