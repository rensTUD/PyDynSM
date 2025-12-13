# Connectivity Analysis Workflow

**Purpose**: Analyze the connectivity of the structural system, assign global DOF indices, build constraint matrices, and classify free/constrained DOFs. This is a **critical step** that must be executed before any analysis operations.

## Entry Point

```python
s1.run_connectivity()
```

## Sequence of Execution

### 1. Entry Point

**Entry**: `s1.run_connectivity()`  
**File**: `src/pydynsm/assembler.py:307-308`

**Call Chain**:
1. `Assembler.run_connectivity()` → delegates to `self.Analysis.connectivity(self.nodes, self.elements)` (line 308)

**Data Passed**:
- Input: `self.nodes` (list of Node), `self.elements` (list of Element)
- Output: None (results stored in `Analysis` instance variables)

**Caching/Reuse**: Results stored in `Analysis` instance and reused for all subsequent analysis operations

**Failure Modes**:
- No explicit error handling - will fail if nodes/elements not properly configured

---

### 2. Connectivity Analysis

**Entry**: `Analysis.connectivity(nodes, elements)`  
**File**: `src/pydynsm/analysis/analysis.py:688-725`

**Call Chain**:
1. `Analysis.connectivity()` orchestrates the full connectivity analysis:
   - `self.assign_dof_indices(nodes, elements)` → line 717
   - `self.build_matrix_B(nodes, elements, self.dof_indices, num_dof)` → line 719
   - `self.find_unique_redundant_dofs(self.B)` → line 721
   - `self.calculate_L(self.B, self.unique_dofs, self.redundant_dofs)` → line 723
   - `self.classify_free_constrained_dofs(nodes, elements, np.array(self.unique_dofs), self.dof_indices)` → line 725

**Data Passed**:
- Input: `nodes` (list), `elements` (list)
- Output: Results stored in `Analysis` instance:
  - `self.dof_indices`: dict mapping `(node_id, element_id)` to dict of DOF names to global indices
  - `self.B`: constraint matrix (numpy.ndarray)
  - `self.unique_dofs`: list of unique DOF indices
  - `self.redundant_dofs`: list of redundant DOF indices
  - `self.L`: transformation matrix (numpy.ndarray)
  - `self.free_dofs`: numpy array of free DOF indices
  - `self.constrained_dofs`: dict mapping constrained DOF indices to their values
  - `self.num_dofs`: total number of DOFs

**Caching/Reuse**: All results cached in `Analysis` instance for use in all subsequent operations

**Failure Modes**:
- Elements without sections/element types → may cause issues in DOF assignment
- No explicit error handling

---

### 3. Assign DOF Indices

**Entry**: `Analysis.assign_dof_indices(nodes, elements)`  
**File**: `src/pydynsm/analysis/analysis.py:585-629`

**Call Chain**:
1. `Analysis.assign_dof_indices()` → line 585
   - Loops through all elements (line 605)
   - For each element, loops through nodes (line 607)
   - For each node, loops through DOFs in element's DOF container (line 613)
   - Assigns sequential global index to each DOF (lines 614-625)
     - Stores in `dof_indices[(node_id, element.id)][dof_name] = global_index`
     - Assigns index to element's DOF: `dof.index = global_index` (line 617)
     - Assigns index to node's DOF if present: `node_dof_container.dofs[dof_name].index = global_index` (lines 621-624)
   - Increments `global_index` for each DOF (line 626)
   - Returns `dof_indices` dict and `num_dof` count (line 629)

**Data Passed**:
- Input: `nodes` (list), `elements` (list)
- Output: 
  - `dof_indices`: defaultdict mapping `(node_id, element_id)` to dict of `{dof_name: global_index}`
  - `num_dof`: total number of DOFs assigned
- Side effects: Updates DOF indices in element and node DOF containers

**Caching/Reuse**: DOF indices stored in multiple places:
- `Analysis.dof_indices` dict
- `element.dof_containers[node.id].dofs[dof_name].index`
- `node.dof_container.dofs[dof_name].index`

**Failure Modes**:
- Elements without proper DOF containers → AttributeError
- No explicit error handling

---

### 4. Build Constraint Matrix B

**Entry**: `Analysis.build_matrix_B(nodes, elements, dof_indices, num_dof)`  
**File**: `src/pydynsm/analysis/analysis.py:632-685`

**Call Chain**:
1. `Analysis.build_matrix_B()` → line 632
   - Estimates number of constraints (line 640)
   - Initializes B matrix: `np.zeros((num_constraints, num_dof))` (line 641)
   - Loops through nodes (line 644)
   - For each node, loops through connected elements (line 645)
   - For each element pair at same node, builds constraint rows (lines 648-682)
     - Checks if DOFs are 'monolithic' (line 655)
     - Creates constraint row: `+1` for first element's DOF, `-1` for other elements' DOFs (lines 656-680)
     - Fills B matrix row (line 681)
   - Returns trimmed B matrix (line 685)

**Data Passed**:
- Input: `nodes`, `elements`, `dof_indices`, `num_dof`
- Output: `B` matrix (numpy.ndarray) where:
  - Rows = constraints (one per DOF connection)
  - Columns = DOFs
  - `B @ u_all = 0` enforces compatibility

**Caching/Reuse**: B matrix stored in `Analysis.B` and used for L matrix calculation

**Failure Modes**:
- Singular or ill-conditioned B matrix → may cause issues in L calculation
- No explicit error handling

---

### 5. Find Unique and Redundant DOFs

**Entry**: `Analysis.find_unique_redundant_dofs(B)`  
**File**: `src/pydynsm/analysis/analysis.py:478-511`

**Call Chain**:
1. `Analysis.find_unique_redundant_dofs()` → line 478
   - Sets `self.num_dofs = B.shape[1]` (line 494)
   - Initializes sets: `redundant_dofs`, `unique_dofs` (lines 495-496)
   - Loops through B matrix rows (line 498)
   - For each row:
     - Finds positive entries → adds to `unique_dofs` (lines 499-503)
     - Finds negative entries → adds to `redundant_dofs` (lines 500-505)
   - Finds zero columns (DOFs with no constraints) → adds to `unique_dofs` (line 508)
   - Removes redundant DOFs from unique set (line 509)
   - Returns sorted lists (line 511)

**Data Passed**:
- Input: `B` matrix (numpy.ndarray)
- Output: 
  - `unique_dofs`: sorted list of unique DOF indices
  - `redundant_dofs`: sorted list of redundant DOF indices
- Side effect: Sets `self.num_dofs`

**Caching/Reuse**: Results stored in `Analysis.unique_dofs` and `Analysis.redundant_dofs`

**Failure Modes**:
- Empty B matrix → may cause issues
- No explicit error handling

---

### 6. Calculate Transformation Matrix L

**Entry**: `Analysis.calculate_L(B, unique_dofs, redundant_dofs)`  
**File**: `src/pydynsm/analysis/analysis.py:513-552`

**Call Chain**:
1. `Analysis.calculate_L()` → line 513
   - If no redundant DOFs: returns identity matrix (lines 549-552)
   - Otherwise:
     - Extracts submatrix `B_ru` (rows for redundant DOFs, columns for unique DOFs) (line 520)
     - Extracts submatrix `B_rr` (rows for redundant DOFs, columns for redundant DOFs) (line 521)
     - Computes `L_ru = -inv(B_rr) @ B_ru` (line 522)
     - Builds L matrix as block matrix (lines 524-548):
       - Identity for unique DOFs
       - `L_ru` for redundant DOFs
     - Returns L matrix (line 548)

**Data Passed**:
- Input: `B` matrix, `unique_dofs` (list), `redundant_dofs` (list)
- Output: `L` matrix (numpy.ndarray) where:
  - `u_all = L @ u_unique` transforms from unique DOF space to full DOF space
  - Shape: `(num_dofs, len(unique_dofs))`

**Caching/Reuse**: L matrix stored in `Analysis.L` and used in all matrix assembly operations

**Failure Modes**:
- Singular `B_rr` matrix → `LinAlgError` from `inv()`
- No explicit error handling - exception propagates

---

### 7. Classify Free and Constrained DOFs

**Entry**: `Analysis.classify_free_constrained_dofs(nodes, elements, unique_dofs, dof_indices)`  
**File**: `src/pydynsm/analysis/analysis.py:787-837`

**Call Chain**:
1. `Analysis.classify_free_constrained_dofs()` → line 787
   - Creates mapping from old DOF indices to new indices in unique_dofs (line 810)
   - Loops through all elements (line 816)
   - For each element, loops through nodes (line 817)
   - For each node, loops through DOFs in `dof_indices` (line 823)
   - Checks if DOF is in unique DOFs (line 825)
   - Gets DOF value from element's DOF container (line 829)
   - Classifies DOF:
     - If `dof_value is None` → adds to `free_dofs` (lines 832-833)
     - Otherwise → adds to `constrained_dofs` dict with value (lines 834-835)
   - Returns numpy array of free DOF indices and dict of constrained DOFs (line 837)

**Data Passed**:
- Input: `nodes`, `elements`, `unique_dofs` (numpy array), `dof_indices` (dict)
- Output:
  - `free_dofs`: numpy array of free DOF indices (subset of unique_dofs)
  - `constrained_dofs`: dict mapping DOF index to prescribed value

**Caching/Reuse**: Results stored in `Analysis.free_dofs` and `Analysis.constrained_dofs`

**Failure Modes**:
- DOF values not properly set → may misclassify DOFs
- No explicit error handling

---

## Summary

**Complete Connectivity Analysis Sequence**:
1. `s1.run_connectivity()` - Entry point
2. `Analysis.connectivity()` - Orchestrates analysis
3. `assign_dof_indices()` - Assigns global DOF indices
4. `build_matrix_B()` - Builds constraint matrix
5. `find_unique_redundant_dofs()` - Identifies DOF types
6. `calculate_L()` - Computes transformation matrix
7. `classify_free_constrained_dofs()` - Separates free/constrained DOFs

**Critical Requirements**:
- Must be called **after** model creation (nodes, elements, sections, element types)
- Must be called **before** any analysis operations (GlobalStiffness, GlobalForce, etc.)
- Only needs to be called **once** per model

**Data Flow**:
- Input: `nodes`, `elements` from Assembler
- Output: All results stored in `Analysis` instance:
  - `dof_indices`: DOF index mapping
  - `B`: Constraint matrix
  - `unique_dofs`, `redundant_dofs`: DOF classification
  - `L`: Transformation matrix
  - `free_dofs`, `constrained_dofs`: DOF constraint classification
  - `num_dofs`: Total DOF count

**Caching/Reuse**:
- All connectivity results cached in `Analysis` instance
- Used by all subsequent analysis operations:
  - `GlobalStiffness()` uses `L`, `unique_dofs`, `redundant_dofs`
  - `GlobalForce()` uses `L`, `unique_dofs`, `redundant_dofs`
  - `GlobalConstrainedStiffness()` uses `free_dofs`
  - `GlobalConstrainedForce()` uses `free_dofs`, `constrained_dofs`
  - `ElementDisplacements()` uses `L`, `unique_dofs`, `redundant_dofs`

**Failure Modes**:
- Singular B matrix → LinAlgError in `calculate_L()`
- Missing DOF indices → AttributeError
- No explicit error handling in most methods - exceptions propagate
