# Static/Dynamic Analysis Workflow

**Purpose**: Assemble global stiffness matrices and force vectors, apply constraints, and solve for free displacements. This workflow is frequency-domain analysis - all operations require an `omega` (frequency) parameter.

## Entry Point

```python
Kc_global = s1.GlobalConstrainedStiffness(omega)
Fc_global = s1.GlobalConstrainedForce(omega)
u_free = s1.SolveUfree(Kc_global, Fc_global)
```

## Sequence of Execution

### 1. Assemble Global Stiffness Matrix

**Entry**: `s1.GlobalStiffness(omega)`  
**File**: `src/pydynsm/assembler.py:92-106`

**Call Chain**:
1. `Assembler.GlobalStiffness(omega)` → delegates to `self.Analysis.GlobalStiffness(self.nodes, self.elements, omega)` (line 106)
2. `Analysis.GlobalStiffness(nodes, elements, omega)` → `src/pydynsm/analysis/analysis.py:20-59`
   - Initializes global stiffness matrix: `k_global = np.zeros((self.num_dofs, self.num_dofs), complex)` (line 39)
   - Loops through elements (line 42):
     - Calls `e.Stiffness(omega)` → `src/pydynsm/analysis/element.py:300-326` (line 43)
       - Gets global and local DOF indices (lines 306-307)
       - Initializes full 12x12 stiffness matrix (line 310)
       - Loops through element types (line 313):
         - Gets element-specific DOFs (line 315)
         - Calls `element_type.LocalStiffness(omega)` → element implementation (e.g., `rod_1D_exp.py`) (line 318)
         - Accumulates stiffness in full matrix (line 318)
       - Applies rotation: `K_global_full = (self.R.T @ K_full) @ self.R` (line 321)
       - Extracts relevant DOFs: `K_global = K_global_full[np.ix_(global_dof_indices, global_dof_indices)]` (line 324)
     - Gets element DOF indices: `idofs = e.get_element_node_dof_indices_global()` (line 44)
     - Assembles into global matrix: `k_global[np.ix_(idofs, idofs)] = elmat` (line 45)
   - Sorts by unique/redundant DOFs (lines 48-51):
     - Extracts submatrices: `k_uu`, `k_ur`, `k_ru`, `k_rr`
     - Reassembles as block matrix `k_global_ur` (lines 53-54)
   - Applies L transformation: `K_global = self.L.T @ k_global_ur @ self.L` (line 57)
   - Returns constrained global stiffness matrix (line 59)

**Data Passed**:
- Input: `omega` (float) - frequency parameter
- Output: `K_global` (numpy.ndarray, complex) - global unconstrained stiffness matrix
- Intermediate data:
  - Element stiffness matrices (12x12 or smaller, complex)
  - Rotation matrices `R` (12x12, real)
  - Transformation matrix `L` (from connectivity analysis)

**Caching/Reuse**:
- Element rotation matrices `R` cached in `element.R` (computed during element creation)
- Transformation matrix `L` cached in `Analysis.L` (computed during connectivity)
- Element stiffness matrices computed fresh for each `omega` (frequency-dependent)
- No caching of global stiffness matrix (recomputed each call)

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.L`, `self.unique_dofs`, etc.)
- Invalid `omega` → may cause numerical issues in element stiffness computation
- Singular element stiffness → may propagate to global matrix
- No explicit error handling

---

### 2. Assemble Global Force Vector

**Entry**: `s1.GlobalForce(omega)`  
**File**: `src/pydynsm/assembler.py:108-122`

**Call Chain**:
1. `Assembler.GlobalForce(omega)` → delegates to `self.Analysis.GlobalForce(self.nodes, self.elements, omega)` (line 122)
2. `Analysis.GlobalForce(nodes, elements, omega)` → `src/pydynsm/analysis/analysis.py:61-122`
   - Initializes global force vector: `f_global = np.zeros(self.num_dofs, complex)` (line 81)
   - Processes distributed loads (lines 84-93):
     - Loops through elements (line 84)
     - Skips if no loads (line 86)
     - Gets element DOF indices (line 90)
     - Calls `element.EvaluateDistributedLoad(element.element_loads, omega)` → `src/pydynsm/analysis/element.py:362-392` (line 93)
       - Gets local DOF indices (line 367)
       - Initializes local load vector (line 373)
       - Loops through element types (line 376):
         - Gets element-specific DOFs (line 378)
         - Evaluates loads (constant or callable `f(omega)`) (lines 381-384)
         - Calls `element.FullDistributedLoad(element_type, q_evaluated, omega)` → element implementation (line 387)
         - Accumulates in local vector (line 387)
       - Transforms to global: `q_glob = self.R.T @ q_loc` (line 390)
       - Returns global force vector (line 392)
     - Adds to global vector: `f_global[np.ix_(dofs)] += element_force` (line 93)
   - Processes nodal loads (lines 96-111):
     - Loops through nodes (line 96)
     - Skips if no loads (line 97)
     - Gets node DOF container (line 101)
     - Loops through nodal loads (line 104):
       - Checks if DOF exists (line 106)
       - Gets DOF index (line 108)
       - Evaluates load (constant or callable `f(omega)`) (line 110)
       - Adds to global vector: `f_global[index] += value` (line 111)
   - Sorts by unique/redundant DOFs (lines 114-117):
     - Extracts subvectors: `f_global_u`, `f_global_r`
     - Concatenates: `f_global_ur = np.hstack((f_global_u, f_global_r))`
   - Applies L transformation: `F_global = self.L.T @ f_global_ur` (line 120)
   - Returns global force vector (line 122)

**Data Passed**:
- Input: `omega` (float) - frequency parameter
- Output: `F_global` (numpy.ndarray, complex) - global unconstrained force vector
- Intermediate data:
  - Element distributed loads (evaluated at `omega`)
  - Nodal loads (evaluated at `omega`)
  - Local force vectors (transformed to global)

**Caching/Reuse**:
- Load functions stored in `node.nodal_loads` and `element.element_loads`
- Loads evaluated fresh for each `omega` (frequency-dependent)
- No caching of global force vector (recomputed each call)

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.L`, `self.unique_dofs`, etc.)
- Invalid `omega` → may cause issues in load evaluation
- Load function errors → exception propagates
- No explicit error handling

---

### 3. Apply Constraints to Stiffness Matrix

**Entry**: `s1.GlobalConstrainedStiffness(omega)`  
**File**: `src/pydynsm/assembler.py:124-138`

**Call Chain**:
1. `Assembler.GlobalConstrainedStiffness(omega)` → delegates to `self.Analysis.GlobalConstrainedStiffness(self.nodes, self.elements, omega)` (line 138)
2. `Analysis.GlobalConstrainedStiffness(nodes, elements, omega)` → `src/pydynsm/analysis/analysis.py:124-143`
   - Calls `self.GlobalStiffness(nodes, elements, omega)` (line 142)
   - Extracts free DOF submatrix: `k_global[np.ix_(self.free_dofs, self.free_dofs)]` (line 143)
   - Returns constrained stiffness matrix (line 143)

**Data Passed**:
- Input: `omega` (float)
- Output: `K_ff` (numpy.ndarray, complex) - constrained stiffness matrix (free DOFs only)
- Uses: `self.free_dofs` from connectivity analysis

**Caching/Reuse**:
- Calls `GlobalStiffness()` which has no caching
- `free_dofs` cached from connectivity analysis

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.free_dofs`)
- Singular constrained matrix → will cause issues in solving
- No explicit error handling

---

### 4. Apply Constraints to Force Vector

**Entry**: `s1.GlobalConstrainedForce(omega)`  
**File**: `src/pydynsm/assembler.py:140-154`

**Call Chain**:
1. `Assembler.GlobalConstrainedForce(omega)` → delegates to `self.Analysis.GlobalConstrainedForce(self.nodes, self.elements, omega)` (line 154)
2. `Analysis.GlobalConstrainedForce(nodes, elements, omega)` → `src/pydynsm/analysis/analysis.py:145-172`
   - Calls `self.GlobalStiffness(nodes, elements, omega)` (line 163)
   - Calls `self.GlobalForce(nodes, elements, omega)` (line 164)
   - Gets constrained DOF indices and values (lines 166-167)
   - Extracts coupling submatrix: `K_free_constrained = k_global[np.ix_(self.free_dofs, constrained_indices)]` (line 169)
   - Extracts free DOF force: `F_free = f_global[self.free_dofs]` (line 170)
   - Applies constraint correction: `F_constrained = F_free - K_free_constrained @ constrained_values` (line 172)
   - Returns constrained force vector (line 172)

**Data Passed**:
- Input: `omega` (float)
- Output: `F_constrained` (numpy.ndarray, complex) - constrained force vector (free DOFs only)
- Uses: `self.free_dofs`, `self.constrained_dofs` from connectivity analysis

**Caching/Reuse**:
- Calls `GlobalStiffness()` and `GlobalForce()` which have no caching
- `free_dofs` and `constrained_dofs` cached from connectivity analysis

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.free_dofs`, `self.constrained_dofs`)
- No explicit error handling

---

### 5. Solve for Free Displacements

**Entry**: `s1.SolveUfree(Kc_global, Fc_global)`  
**File**: `src/pydynsm/assembler.py:174-190`

**Call Chain**:
1. `Assembler.SolveUfree(Kc_global, Fc_global)` → delegates to `self.Analysis.SolveUfree(Kc_global, Fc_global)` (line 190)
2. `Analysis.SolveUfree(Kc_global, Fc_global)` → `src/pydynsm/analysis/analysis.py:208-224`
   - Solves linear system: `np.linalg.solve(Kc_global, Fc_global)` (line 224)
   - Returns free displacements (line 224)

**Data Passed**:
- Input: 
  - `Kc_global` (numpy.ndarray, complex) - constrained stiffness matrix
  - `Fc_global` (numpy.ndarray, complex) - constrained force vector
- Output: `u_free` (numpy.ndarray, complex) - free displacements

**Caching/Reuse**: None - solves fresh each time

**Failure Modes**:
- Singular matrix → `LinAlgError` from `np.linalg.solve()`
- Ill-conditioned matrix → numerical errors
- Dimension mismatch → `ValueError`
- No explicit error handling - exception propagates

---

### 6. Compute Full Displacement Vector (Optional)

**Entry**: `s1.FullDisplacement(u_free)`  
**File**: `src/pydynsm/assembler.py:212-226`

**Call Chain**:
1. `Assembler.FullDisplacement(u_free)` → delegates to `self.Analysis.FullDisplacement(u_free)` (line 226)
2. `Analysis.FullDisplacement(u_free)` → `src/pydynsm/analysis/analysis.py:252-274`
   - Gets constrained DOF indices and values (lines 266-267)
   - Initializes full displacement vector (line 269)
   - Sets free DOFs: `u_full[np.ix_(self.free_dofs)] = u_free` (line 271)
   - Sets constrained DOFs: `u_full[np.ix_(constrained_indices)] = constrained_values` (line 272)
   - Returns full displacement vector (line 274)

**Data Passed**:
- Input: `u_free` (numpy.ndarray, complex) - free displacements
- Output: `u_full` (numpy.ndarray, complex) - full displacement vector (free + constrained)
- Uses: `self.free_dofs`, `self.constrained_dofs` from connectivity analysis

**Caching/Reuse**: None

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.free_dofs`, `self.constrained_dofs`)
- Dimension mismatch → `IndexError`
- No explicit error handling

---

### 7. Compute Support Reactions (Optional)

**Entry**: `s1.SupportReactions(k_global, u_free, f_global)`  
**File**: `src/pydynsm/assembler.py:192-210`

**Call Chain**:
1. `Assembler.SupportReactions(k_global, u_free, f_global)` → delegates to `self.Analysis.SupportReactions(k_global, u_free, f_global)` (line 210)
2. `Analysis.SupportReactions(k_global, u_free, f_global)` → `src/pydynsm/analysis/analysis.py:226-250`
   - Gets constrained DOF indices and values (lines 244-245)
   - Extracts submatrices:
     - `Kcf = k_global[np.ix_(constrained_indices, self.free_dofs)]` (line 247)
     - `Kcc = k_global[np.ix_(constrained_indices, constrained_indices)]` (line 248)
   - Computes reactions: `(Kcf @ u_free) + (Kcc @ constrained_values) - f_global[constrained_indices]` (line 250)
   - Returns support reactions (line 250)

**Data Passed**:
- Input:
  - `k_global` (numpy.ndarray, complex) - unconstrained global stiffness matrix
  - `u_free` (numpy.ndarray, complex) - free displacements
  - `f_global` (numpy.ndarray, complex) - unconstrained global force vector
- Output: `reactions` (numpy.ndarray, complex) - support reactions at constrained DOFs
- Uses: `self.free_dofs`, `self.constrained_dofs` from connectivity analysis

**Caching/Reuse**: None

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.free_dofs`, `self.constrained_dofs`)
- Dimension mismatch → `IndexError`
- No explicit error handling

---

## Summary

**Complete Analysis Sequence**:
1. `s1.GlobalStiffness(omega)` - Assemble unconstrained stiffness matrix
2. `s1.GlobalForce(omega)` - Assemble unconstrained force vector
3. `s1.GlobalConstrainedStiffness(omega)` - Extract free DOF stiffness submatrix
4. `s1.GlobalConstrainedForce(omega)` - Compute constrained force vector
5. `s1.SolveUfree(Kc_global, Fc_global)` - Solve for free displacements
6. `s1.FullDisplacement(u_free)` - Reconstruct full displacement vector (optional)
7. `s1.SupportReactions(k_global, u_free, f_global)` - Compute support reactions (optional)

**Critical Requirements**:
- `run_connectivity()` must be called first
- All operations require `omega` (frequency) parameter
- Stiffness and force matrices must be computed for same `omega`

**Data Flow**:
- Input: `omega` (frequency parameter)
- Intermediate: Element stiffness matrices, force vectors
- Output: 
  - `K_global`, `F_global`: Unconstrained matrices
  - `Kc_global`, `Fc_global`: Constrained matrices
  - `u_free`: Free displacements
  - `u_full`: Full displacements (optional)
  - `reactions`: Support reactions (optional)

**Caching/Reuse**:
- **No caching** of global matrices (recomputed each call)
- Element rotation matrices `R` cached (computed once)
- Transformation matrix `L` cached (from connectivity)
- DOF classifications cached (from connectivity)
- Element stiffness matrices computed fresh (frequency-dependent)

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing L, free_dofs, etc.)
- Singular stiffness matrix → `LinAlgError` in `SolveUfree()`
- Ill-conditioned matrix → numerical errors
- Invalid `omega` → may cause numerical issues
- No explicit error handling in most methods
