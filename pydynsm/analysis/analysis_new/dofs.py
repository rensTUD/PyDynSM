# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 12:44:33 2024

@author: rvanderleijden
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional
from typing import Dict

# %% DOF class

@dataclass
class DOF:
    name: str
    value: Optional[float] = None # Prescribed value or None if free
    index: Optional[int] = None   # Global index of dof after running connectivity

# %% DOFcontainer class

class DOFContainer:
    def __init__(self):
        self.dofs: Dict[str, DOF] = {}

    def set_dof(self, dof_name: str, value: Optional[float] = None):
        """Set or update the value of a dof DOF."""
        self.dofs[dof_name] = DOF(name=dof_name, value=value)

    def get_dof(self, dof_name: str) -> Optional[DOF]:
        """Retrieve a DOF."""
        return self.dofs.get(dof_name)

    def has_dof(self, dof_name: str) -> bool:
        """Check if a DOF exists."""
        return dof_name in self.dofs

