# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 19:21:45 2024

@author: rensv
"""

# %% import dependencies
import numpy as np
import matplotlib.pyplot as plt

# %% main class definition

class Assembler:
    
    def __init__(self, name):
        '''
        Initialisation of assembler.
        '''
        self.elements = []  # Empty elements list
        self.nodes = []  # Empty nodes list
        self.name = name  # Assembler / project name
        
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
    
    # TODO - SHOULD BE MOVED TO A SEPARATE CLASS CALLED PLOTTER OR SOMETHING
    def PlotNodes(self, plot_elements=False):
        '''
        Plots the nodes and elements of the structure on a 2D plot.
        '''
        plt.figure(figsize=(10, 6))
        
        # Plot nodes
        x_coords = [node.x for node in self.nodes]
        z_coords = [node.z for node in self.nodes]
        plt.scatter(x_coords, z_coords, color='blue', marker='o', label='Nodes')
        
        # Label nodes
        for i, node in enumerate(self.nodes):
            plt.text(node.x, node.z, f'Node {i+1}', fontsize=9, ha='right')
        
        if plot_elements:
            # Plot elements as lines between nodes
            for element in self.elements:
                x_values = [element.nodes[0].x, element.nodes[1].x]
                z_values = [element.nodes[0].z, element.nodes[1].z]
                plt.plot(x_values, z_values, 'k-', linewidth=2, label='Element' if element == self.elements[0] else "")
        
        plt.title(f'Structure Configuration: {self.name}')
        plt.xlabel('X Coordinate')
        plt.ylabel('Z Coordinate')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.show()
    
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
    
        
            
        
