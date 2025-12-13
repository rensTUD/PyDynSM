# Post-Processing Workflow

**Purpose**: Compute element-level results (displacements, forces, stresses) from global displacement solution. Results are computed at multiple points along each element for visualization.

## Entry Point

```python
u_elem = s1.FullDisplacement(u_free)
displacements = s1.ElementDisplacements(u_elem, omega, num_points=100)
forces = s1.ElementForces(u_elem, omega, num_points=100)
stresses = s1.ElementStresses(u_elem, omega, num_points=100)
```

## Sequence of Execution

### 1. Get Full Displacement Vector

**Entry**: `s1.FullDisplacement(u_free)`  
**File**: `src/pydynsm/assembler.py:212-226`

**Call Chain**: (Same as in static_dynamic_analysis.md)
1. `Assembler.FullDisplacement(u_free)` → `self.Analysis.FullDisplacement(u_free)` (line 226)
2. `Analysis.FullDisplacement(u_free)` → `src/pydynsm/analysis/analysis.py:252-274`
   - Reconstructs full displacement vector including constrained DOFs
   - Returns `u_full` in unique DOF space

**Data Passed**:
- Input: `u_free` (numpy.ndarray, complex) - free displacements
- Output: `u_full` (numpy.ndarray, complex) - full displacement vector in unique DOF space

**Caching/Reuse**: None

**Failure Modes**: See static_dynamic_analysis.md

---

### 2. Compute Element Displacements

**Entry**: `s1.ElementDisplacements(u_full, omega, num_points=20, local_axes=False)`  
**File**: `src/pydynsm/assembler.py:228-233`

**Call Chain**:
1. `Assembler.ElementDisplacements(u_full, omega, num_points, local_axes)` → delegates to `self.Analysis.ElementDisplacements(self.elements, u_full, omega, num_points, local_axes)` (line 233)
2. `Analysis.ElementDisplacements(elements, u_nodes_global, omega, num_points, local_axes)` → `src/pydynsm/analysis/analysis.py:326-374`
   - Applies L matrix to get all DOF displacements: `u_nodes_global_all = self.L @ u_nodes_global` (line 348)
   - Gets current DOF order: `current_order = self.unique_dofs + self.redundant_dofs` (line 350)
   - Creates index map to reorder to sorted order (lines 352-354)
   - Reorders displacement vector: `u_nodes_global_all_sorted = np.array(u_nodes_global_all)[index_map]` (line 357)
   - Initializes result dict: `element_displacements = {}` (line 359)
   - Loops through elements (line 361):
     - Gets element DOF indices: `global_dof_indices = element.get_element_node_dof_indices_global()` (line 363)
     - Extracts element displacements: `u_element_global = u_nodes_global_all_sorted[global_dof_indices]` (line 366)
     - Calls `element.Displacements(u_element_global, omega, num_points, local_axes)` → `src/pydynsm/analysis/element.py:799-892` (line 369)
       - Transforms to local coordinates: `u_nodes_local = self.R @ u_element_global` (line 808)
       - Extracts local DOF indices (line 811)
       - Loops through element types (line 814):
         - Gets element-specific DOFs (line 816)
         - Extracts local nodal displacements for element type (line 819)
         - Calls `element_type.LocalDisplacements(u_nodes_local_element_type, omega, num_points)` → element implementation (line 822)
         - Accumulates displacements (line 825)
       - If `local_axes=False`: transforms back to global (lines 828-891)
         - Loops through points (line 830)
         - Transforms displacement vector at each point (line 832)
       - Returns displacement array: shape `(6, num_points)` for global, `(num_local_dofs, num_points)` for local (line 892)
     - Stores in dict: `element_displacements[element.id] = u_elem` (line 372)
   - Returns dict of element displacements (line 374)

**Data Passed**:
- Input:
  - `u_full` (numpy.ndarray, complex) - full displacement vector in unique DOF space
  - `omega` (float) - frequency parameter
  - `num_points` (int) - number of evaluation points along element
  - `local_axes` (bool) - whether to return in local coordinates
- Output: `element_displacements` (dict) - maps `element.id` to displacement array
  - Shape: `(6, num_points)` for global axes, `(num_local_dofs, num_points)` for local axes
  - Each column is displacement at a point along element length
- Intermediate data:
  - `u_nodes_global_all`: Full displacement vector (unique + redundant DOFs)
  - `u_element_global`: Element nodal displacements in global coordinates
  - `u_nodes_local`: Element nodal displacements in local coordinates
  - Element type local displacements

**Caching/Reuse**:
- L matrix cached from connectivity analysis
- Element rotation matrices `R` cached from element creation
- Element displacement computation done fresh for each call

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.L`, `self.unique_dofs`, etc.)
- Invalid `num_points` → may cause issues
- Element without element types → may cause issues in displacement computation
- No explicit error handling

---

### 3. Compute Element Forces

**Entry**: `s1.ElementForces(u_full, omega, num_points=20)`  
**File**: `src/pydynsm/assembler.py:234-239`

**Call Chain**:
1. `Assembler.ElementForces(u_full, omega, num_points)` → delegates to `self.Analysis.ElementForces(self.elements, u_full, omega, num_points)` (line 239)
2. `Analysis.ElementForces(elements, u_nodes_global, omega, num_points)` → `src/pydynsm/analysis/analysis.py:376-424`
   - Same structure as `ElementDisplacements()`:
     - Applies L matrix: `u_nodes_global_all = self.L @ u_nodes_global` (line 398)
     - Reorders to sorted DOF order (lines 400-407)
     - Loops through elements (line 411):
       - Gets element DOF indices (line 413)
       - Extracts element displacements (line 416)
       - Calls `element.Forces(u_element_global, omega, num_points)` → `src/pydynsm/analysis/element.py:895-947` (line 419)
         - Transforms to local coordinates: `u_nodes_local = self.R @ u_element_global` (line 904)
         - Extracts local DOF indices (line 907)
         - Loops through element types (line 910):
           - Gets element-specific DOFs (line 912)
           - Extracts local nodal displacements (line 915)
           - Calls `element_type.LocalForces(u_nodes_local_element_type, omega, num_points)` → element implementation (line 918)
           - Accumulates forces (line 921)
         - Transforms back to global (lines 924-945)
         - Returns force array: shape `(6, num_points)` (line 947)
       - Stores in dict: `element_forces[element.id] = f_elem` (line 422)
   - Returns dict of element forces (line 424)

**Data Passed**:
- Input: Same as `ElementDisplacements()`
- Output: `element_forces` (dict) - maps `element.id` to force array
  - Shape: `(6, num_points)` - forces at each point along element
  - Rows: [Fx, Fz, Fy, Mx, Mz, My] in global coordinates

**Caching/Reuse**: Same as `ElementDisplacements()`

**Failure Modes**: Same as `ElementDisplacements()`

---

### 4. Compute Element Stresses

**Entry**: `s1.ElementStresses(u_full, omega, num_points=20)`  
**File**: `src/pydynsm/assembler.py:240-245`

**Call Chain**:
1. `Assembler.ElementStresses(u_full, omega, num_points)` → delegates to `self.Analysis.ElementStresses(self.elements, u_full, omega, num_points)` (line 245)
2. `Analysis.ElementStresses(elements, u_nodes_global, omega, num_points)` → `src/pydynsm/analysis/analysis.py:426-474`
   - Same structure as `ElementDisplacements()` and `ElementForces()`:
     - Applies L matrix: `u_nodes_global_all = self.L @ u_nodes_global` (line 448)
     - Reorders to sorted DOF order (lines 450-457)
     - Loops through elements (line 461):
       - Gets element DOF indices (line 463)
       - Extracts element displacements (line 466)
       - Calls `element.Stresses(u_element_global, omega, num_points)` → `src/pydynsm/analysis/element.py:949-1007` (line 469)
         - Transforms to local coordinates: `u_nodes_local = self.R @ u_element_global` (line 958)
         - Extracts local DOF indices (line 961)
         - Loops through element types (line 964):
           - Gets element-specific DOFs (line 966)
           - Extracts local nodal displacements (line 969)
           - Calls `element_type.LocalStresses(u_nodes_local_element_type, omega, num_points)` → element implementation (line 972)
           - Accumulates stresses (line 975)
         - Returns stress array: shape `(num_stress_components, num_points)` (line 1007)
       - Stores in dict: `element_stresses[element.id] = s_elem` (line 472)
   - Returns dict of element stresses (line 474)

**Data Passed**:
- Input: Same as `ElementDisplacements()`
- Output: `element_stresses` (dict) - maps `element.id` to stress array
  - Shape: `(num_stress_components, num_points)` - stresses at each point
  - Components depend on element type (e.g., axial, bending, shear)

**Caching/Reuse**: Same as `ElementDisplacements()`

**Failure Modes**: Same as `ElementDisplacements()`

---

## Summary

**Complete Post-Processing Sequence**:
1. `s1.FullDisplacement(u_free)` - Reconstruct full displacement vector
2. `s1.ElementDisplacements(u_full, omega, num_points)` - Compute element displacements
3. `s1.ElementForces(u_full, omega, num_points)` - Compute element forces (optional)
4. `s1.ElementStresses(u_full, omega, num_points)` - Compute element stresses (optional)

**Critical Requirements**:
- `run_connectivity()` must be called first
- `u_free` must be from solved system (from `SolveUfree()`)
- All operations require `omega` parameter (must match analysis `omega`)

**Data Flow**:
- Input: `u_free` (from analysis), `omega` (frequency), `num_points` (resolution)
- Processing:
  - Apply L matrix to get all DOF displacements
  - Extract element-level displacements
  - Transform to local coordinates
  - Compute element-level results (displacements/forces/stresses)
  - Transform back to global coordinates (if needed)
- Output: Dicts mapping `element.id` to result arrays
  - Displacements: `(6, num_points)` or `(num_local_dofs, num_points)`
  - Forces: `(6, num_points)`
  - Stresses: `(num_stress_components, num_points)`

**Caching/Reuse**:
- L matrix cached from connectivity analysis
- Element rotation matrices `R` cached from element creation
- Element-level computations done fresh for each call
- No caching of post-processing results

**Failure Modes**:
- Connectivity not run → `AttributeError` (missing `self.L`, etc.)
- Invalid `u_free` → dimension mismatch errors
- `omega` mismatch → incorrect results (no validation)
- Element without element types → may cause issues
- Invalid `num_points` → may cause issues
- No explicit error handling in most methods
