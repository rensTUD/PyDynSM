# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 19:21:45 2024

@author: rensv
"""

# %% import dependencies
import numpy as np


# %% main class definition

class Assembler:
    
    def __init__(self, name):
        '''
        Initialisation of assembler
        '''
                        
        # initialise empty elements list
        self.elements = []
        
        # empty nodes list
        self.nodes = []
        
        # assembler / project name
        self.name = name
        
        print(f"Assembler '{self.name}' successfully initialised")
        
    def RegisterNode(self, node):
        if node not in self.nodes:
            self.nodes.append(node)

    def RegisterElement(self, element):
        if element not in self.elements:
            self.elements.append(element)
    
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
    
        
            
        
