# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 19:21:45 2024

@author: rensv
"""

# %% import dependencies


# %% import others

from . import plotter
from .analysis import analysis_old
from .analysis import analysis_new
# %% main class definition

class Assembler:
        
    def __init__(self, name, analysis_type = 'old'):
        '''
        Initialisation of assembler.
        '''
        
        self.name = name  # Assembler / project name
        
        self.elements = []  # Empty elements list
        
        self.nodes = []  # Empty nodes list
       
        
        # load dependencies that are injected
        self.StructurePlotter = plotter.StructurePlotter()
        if analysis_type == 'old':
            analysis_type = analysis_old
        elif analysis_type == 'new':
            analysis_type = analysis_new
            
        self.Node = analysis_type.Node
        self.Element = analysis_type.Element
        self.Analysis = analysis_type.Analysis()
            
        print(f"Assembler '{self.name}' successfully initialised")

    
# %% Nodes 
    
    def RegisterNode(self, node):
        '''
        Adds a node to the assembler if it's not already registered.
        '''
        if node not in self.nodes:
            self.nodes.append(node)          # add node to list
    
    def CreateNode(self, *args, **kwargs):
        """Creates a node and registers it automatically with the assembler."""
        new_node = self.Node(*args, **kwargs)
        self.RegisterNode(new_node)
        return new_node        

# %% Elements
  
    def RegisterElement(self, element):
        '''
        Adds an element to the assembler if it's not already registered.
        '''
        if element not in self.elements:
            self.elements.append(element)
    
    def CreateElement(self, nodes, element_type=None, props={}):
        """Creates an element and registers it automatically with the assembler."""
        new_element = self.Element(nodes)
        self.RegisterElement(new_element)
        
        # TODO - if element_type is given, set the section as well - need to test this
        if element_type:
            new_element.SetSection(element_type, props)
        return new_element     
    
# %% Plotting    
    def PlotStructure(self, plot_elements=False):
        
        self.StructurePlotter.PlotNodes(self.nodes)
        
        if plot_elements:
            self.StructurePlotter.PlotElements(self.elements)
            
        self.StructurePlotter.ShowStructure(f'Structure Configuration: {self.name}')

# %% ANALAYSIS 
# i.e. Stiffness, Force, Displacements

    def GlobalStiffness(self, omega):
        """
        Assembles the global unconstrained stiffness matrix.

        Parameters
        ----------
        omega : float
            Frequency parameter.

        Returns
        -------
        k_global : numpy.ndarray
            Global stiffness matrix.
        """
        return self.Analysis.GlobalStiffness(self.nodes, self.elements, omega)

    def GlobalForce(self, omega):
        """
        Assembles the global unconstrained force vector.

        Parameters
        ----------
        omega : float
            Frequency parameter.

        Returns
        -------
        f_global : numpy.ndarray
            Global force vector.
        """
        return self.Analysis.GlobalForce(self.nodes, self.elements, omega)

    def GlobalConstrainedStiffness(self, omega):
        """
        Constrains the global stiffness matrix.

        Parameters
        ----------
        omega : float
            Frequency parameter.

        Returns
        -------
        k_constrained : numpy.ndarray
            Constrained global stiffness matrix.
        """
        return self.Analysis.GlobalConstrainedStiffness(self.nodes, self.elements, omega)

    def GlobalConstrainedForce(self, omega):
        """
        Constrains the global force vector.

        Parameters
        ----------
        omega : float
            Frequency parameter.

        Returns
        -------
        f_constrained : numpy.ndarray
            Constrained global force vector.
        """
        return self.Analysis.GlobalConstrainedForce(self.nodes, self.elements, omega)

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
        return self.Analysis.SolveUfree(Kc_global, fc_global)

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
        return self.Analysis.SupportReactions(k_global, u_free, f_global)

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
        return self.Analysis.FullDisplacement(u_free)
    
    def ElementDisplacements(self, u_full, omega, num_points=20):
        '''
        Returns all element displacements
        '''
        
        return self.Analysis.ElementDisplacements(self.elements, u_full, omega, num_points)
    def ElementForces(self, u_full, omega, num_points=20):
        '''
        Returns all element forces
        '''
        
        return self.Analysis.ElementForces(self.elements, u_full, omega, num_points)
    def ElementStresses(self, u_full, omega, num_points=20):
        '''
        Returns all element stresses
        '''
        
        return self.Analysis.ElementStresses(self.elements, u_full, omega, num_points)
    
# %% POSTPROCESS

    def PlotElementDisplacements(self,Element_displacements: dict, scale: float = 1.0) -> None:
        '''
        Plots the displacements
        '''
        self.StructurePlotter.PlotDisplacements(self.elements, Element_displacements, scale=scale)
        
    
    def PlotMoments(self,Element_forces: dict, scale: float = 1.0) -> None:
         '''
         Plots the moments
         '''
         self.StructurePlotter.Plotmoments(self.elements, Element_forces, scale=scale)
    
    def PlotAxialforces(self,Element_forces: dict, scale: float = 1.0) -> None:
         '''
         Plots the axial forces
         '''
         self.StructurePlotter.Plotaxialforces(self.elements, Element_forces, scale=scale)
    
    def PlotShearforces(self,Element_forces: dict, scale: float = 1.0) -> None:
         '''
         Plots the shear forces
         '''
         self.StructurePlotter.Plotshearforces(self.elements, Element_forces, scale=scale)
    
    def PlotShearstresses(self,Element_stresses: dict, scale: float = 1.0) -> None:
         '''
         Plots the shear forces
         '''
         self.StructurePlotter.Plotshearstresses(self.elements, Element_stresses, scale=scale)
    def PlotAxialstresses(self,Element_stresses: dict, scale: float = 1.0) -> None:
         '''
         Plots the shear forces
         '''
         self.StructurePlotter.Plotaxialstresses(self.elements, Element_stresses, scale=scale)

# %% Testing area
    def get_dofs_elements(self):
        return self.Analysis.assign_dof_indices(self.nodes, self.elements)

    def get_B_matrix(self):
        dof_indices, num_dof = self.get_dofs_elements()
        
        return dof_indices, self.Analysis.build_matrix_B(self.nodes, self.elements, dof_indices,num_dof)
    
    def run_connectivity(self):
        self.Analysis.connectivity(self.nodes, self.elements)
        
# %% TODO's        
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
    
        
            
        
