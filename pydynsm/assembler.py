# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 19:21:45 2024

@author: rensv
"""

# %% import dependencies
import numpy as np

from .constrainer import Constrainer

# %% main class definition

class Assembler:
    
    def __init__(self):
        '''
        Initialisation of assembler
        '''
        
        # # initialise dependencies
        # self.Constrainer = Constrainer()
        
        # initialise empty elements list
        self.elements = []
    
    def AddElements(self, elements):
        '''
        Function that adds element(s) to itself
        '''
        
        # check if input is a list or not
        if not isinstance(elements,list):
            elements = [elements]
            
        self.elements.extend(elements)
        
