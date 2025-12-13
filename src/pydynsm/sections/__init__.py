# -*- coding: utf-8 -*-
"""
Created on [Date]

@author: rensv
"""

from .section import Section, SectionFactory

# Import all section classes to register them
from .rectangle import Rectangle
from .circle import Circle
from .hollow_circle import HollowCircle
from .i_section import ISection

__all__ = ['Section', 'SectionFactory', 'Rectangle', 'Circle', 'HollowCircle', 'ISection']
