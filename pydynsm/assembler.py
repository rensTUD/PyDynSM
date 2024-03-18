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
    
    def __init__(self, name, plotter=plotter):
        '''
        Initialisation of assembler.
        '''
        self.elements = []  # Empty elements list
        self.nodes = []  # Empty nodes list
        self.name = name  # Assembler / project name
        
        
        # load dependencies that are injected
        self.StructurePlotter = plotter.StructurePlotter()
        
        print(f"Assembler '{self.name}' successfully initialised")
        
    def RegisterNode(self, node):
        '''
        Adds a node to the assembler if it's not already registered.
        '''
        if node not in self.nodes:
            self.nodes.append(node)

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
        
        # plt.figure(figsize=(10, 6))

        # for node in self.nodes:
        #     self.StructurePlotter.PlotNode(node)
        
        # if plot_elements:
        #     for element in self.elements:
        #         self.StructurePlotter.PlotElement(element)
        
        # self.StructurePlotter.ShowStructure(f'Structure Configuration: {self.name}')
    
    
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
    
        
            
        
