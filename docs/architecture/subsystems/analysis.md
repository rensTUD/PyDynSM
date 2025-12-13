# Analysis Subsystem

**Location**: `src/pydynsm/analysis/`

## Responsibilities

The Analysis subsystem is the core computational engine of PyDynSM. It handles:

- DOF (Degree of Freedom) assignment and indexing
- Constraint matrix (B) construction for multi-element connections
- Transformation matrix (L) computation for redundant DOF elimination
- Global stiffness matrix assembly
- Global force vector assembly
- Constraint application (free vs constrained DOFs)
- Linear system solving
- Post-processing (element displacements, forces, stresses)

## Public API

### Analysis Class (`analysis.py`)

#### Connectivity and DOF Management

- `connectivity(nodes, elements)`: Main entry point - computes all connectivity information
  - Calls `assign_dof_indices()`, `build_matrix_B()`, `find_unique_redundant_dofs()`, `calculate_L()`, `classify_free_constrained_dofs()`
  - Stores results in instance variables: `dof_indices`, `B`, `unique_dofs`, `redundant_dofs`, `L`, `free_dofs`, `constrained_dofs`

- `assign_dof_indices(nodes, elements) -> tuple[dict, int]`: Assigns global DOF indices to all element DOFs
- `build_matrix_B(nodes, elements, dof_indices, num_dof) -> np.ndarray`: Builds constraint matrix B for monolithic connections
- `find_unique_redundant_dofs(B) -> tuple[list, list]`: Identifies unique and redundant DOFs from constraint matrix
- `calculate_L(B, unique_dofs, redundant_dofs) -> np.ndarray`: Computes transformation matrix L
- `classify_free_constrained_dofs(nodes, elements, unique_dofs, dof_indices) -> tuple[np.ndarray, dict]`: Separates free vs constrained DOFs

#### Matrix Assembly

- `GlobalStiffness(nodes, elements, omega) -> np.ndarray`: Assembles global unconstrained stiffness matrix
- `GlobalForce(nodes, elements, omega) -> np.ndarray`: Assembles global unconstrained force vector
- `GlobalConstrainedStiffness(nodes, elements, omega) -> np.ndarray`: Returns constrained stiffness matrix (free DOFs only)
- `GlobalConstrainedForce(nodes, elements, omega) -> np.ndarray`: Returns constrained force vector (free DOFs only)
- `GlobalConstrainedSystem(nodes, elements, omega) -> tuple[np.ndarray, np.ndarray]`: Returns (K_ff, F_constrained)

#### Solving

- `SolveUfree(Kc_global, fc_global) -> np.ndarray`: Solves linear system for free displacements
- `SupportReactions(k_global, u_free, f_global) -> np.ndarray`: Calculates support reactions at constrained DOFs
- `FullDisplacement(u_free) -> np.ndarray`: Reconstructs full displacement vector (free + constrained)
- `SolveEigenvector(nodes, elements, omega, fixed_index=0) -> np.ndarray`: Computes eigenvector at specific frequency

#### Post-processing

- `ElementDisplacements(elements, u_nodes_global, omega, num_points=20, local_axes=False) -> dict`: Returns displacements for each element
- `ElementForces(elements, u_nodes_global, omega, num_points=20) -> dict`: Returns forces for each element
- `ElementStresses(elements, u_nodes_global, omega, num_points=20) -> dict`: Returns stresses for each element

### Element Class (`element.py`)

- `Element.__init__(nodes)`: Creates element, connects to nodes, initializes DOF containers
- `Element.Stiffness(omega) -> np.ndarray`: Assembles element stiffness matrix from all element types
- `Element.SetSection(section_type, props)`: Sets cross-section via SectionFactory
- `Element.SetElementType(element_type, **props)`: Adds structural element type via ElementFactory
- `Element.AddDistributedLoad(**loads)`: Adds distributed loads
- `Element.EvaluateDistributedLoad(element_loads, omega) -> np.ndarray`: Evaluates distributed loads to nodal loads
- `Element.Displacements(u_nodes_global, omega, num_points=20, local_axes=False) -> np.ndarray`: Computes element displacements
- `Element.Forces(u_nodes_global, omega, num_points=20) -> np.ndarray`: Computes element forces
- `Element.Stresses(u_nodes_global, omega, num_points=20) -> np.ndarray`: Computes element stresses
- `Element.prescribe_dof(node, **dofs)`: Prescribes DOF values
- `Element.fix_dof(node, *dofs)`: Fixes DOFs (sets to 0)
- `Element.free_dof(node, *dofs)`: Frees DOFs (sets to None)
- `Element.couple_dof(node, *dofs)`: Ensures monolithic connection
- `Element.decouple_dof(node, *dofs)`: Allows independent movement

### Node Class (`node.py`)

- `Node.__init__(x, z, y=0, config='2D', dof_config=None)`: Creates node with coordinates and DOF configuration
- `Node.prescribe_node(**dofs)`: Prescribes DOF values, propagates to connected elements
- `Node.fix_node(*dofs)`: Fixes DOFs (sets to 0)
- `Node.free_node(*dofs)`: Frees DOFs (sets to None)
- `Node.add_load(**loads)`: Adds nodal loads
- `Node.remove_load(dof)`: Removes nodal load
- `Node.connect_element(element)`: Registers connected element
- `Node.get_coords() -> tuple`: Returns (x, z, y) coordinates
- `Node.get_dof_indices() -> list[int]`: Returns list of global DOF indices

### DOF Classes (`dofs.py`)

- `DOF(name: str, value: Optional[float] = None, index: Optional[int] = None)`: Dataclass representing a degree of freedom
- `DOFContainer`: Container for DOF objects
  - `set_dof(dof_name, value=None)`: Creates/updates DOF
  - `get_dof(dof_name) -> Optional[DOF]`: Retrieves DOF
  - `has_dof(dof_name) -> bool`: Checks existence
  - `set_dof_value(dof_name, value)`: Updates DOF value

## Key Data Structures

### Analysis Instance Variables

- `self.dof_indices: dict`: Maps `(node_id, element_id)` to dict of DOF names to global indices
- `self.B: np.ndarray`: Constraint matrix (rows = constraints, cols = DOFs)
- `self.unique_dofs: list[int]`: List of unique DOF indices
- `self.redundant_dofs: list[int]`: List of redundant DOF indices
- `self.L: np.ndarray`: Transformation matrix (maps unique DOFs to all DOFs)
- `self.free_dofs: np.ndarray`: Array of free DOF indices (subset of unique_dofs)
- `self.constrained_dofs: dict[int, float]`: Maps constrained DOF indices to their prescribed values
- `self.num_dofs: int`: Total number of DOFs

### Element Data Structures

- `element.nodes: list[Node]`: Connected nodes
- `element.dof_containers: dict[int, DOFContainer]`: Global DOF containers per node
- `element.local_dof_container: dict[int, DOFContainer]`: Local DOF containers per node
- `element.constraint_types: dict[int, dict[str, str]]`: Constraint type per node/DOF ('monolithic' or 'independent')
- `element.element_types: dict[str, StructuralElement]`: Dictionary of element type instances
- `element.section: Section`: Cross-section instance
- `element.element_loads: defaultdict(dict)`: Distributed loads per DOF
- `element.R: np.ndarray`: 12x12 rotation matrix (global to local transformation)
- `element.L: float`: Element length

### Node Data Structures

- `node.dof_container: DOFContainer`: Container for node's DOFs
- `node.nodal_loads: defaultdict(dict)`: Nodal loads per DOF
- `node.connected_elements: list[Element]`: List of connected elements
- `node.x, node.z, node.y: float`: Coordinates
- `node.id: int`: Unique node identifier (class variable `Node.nn`)

## Lifecycle

### Analysis Lifecycle

1. **Initialization**: `Analysis()` instance created (no parameters)
2. **Connectivity**: `connectivity(nodes, elements)` called
   - DOF indices assigned
   - Constraint matrix B built
   - Unique/redundant DOFs identified
   - Transformation matrix L computed
   - Free/constrained DOFs classified
3. **Assembly**: `GlobalStiffness()` / `GlobalForce()` called
   - Element stiffness matrices assembled
   - Transformation matrix L applied
4. **Constraint Application**: `GlobalConstrainedStiffness()` / `GlobalConstrainedForce()` called
   - Free DOF submatrices extracted
   - Constrained DOF contributions subtracted from force vector
5. **Solving**: `SolveUfree()` called
   - Linear system solved
6. **Post-processing**: `ElementDisplacements()` / `ElementForces()` / `ElementStresses()` called
   - L matrix applied to get all DOF displacements
   - Element-level computations performed

### Element Lifecycle

1. **Creation**: `Element(nodes)` called
   - Nodes connected bidirectionally
   - DOF containers initialized from nodes
   - Geometrical properties computed (length, rotation matrix)
2. **Section Setup**: `SetSection(section_type, props)` called
   - Section instance created via SectionFactory
3. **Element Type Setup**: `SetElementType(element_type, **props)` called
   - StructuralElement instance created via ElementFactory
   - Local DOFs added to local_dof_container
   - Global constraints applied to local DOFs
4. **Load/Constraint Setup**: `AddDistributedLoad()`, `prescribe_dof()`, etc. called
5. **Analysis**: `Stiffness(omega)` called during global assembly
6. **Post-processing**: `Displacements()`, `Forces()`, `Stresses()` called

### Node Lifecycle

1. **Creation**: `Node(x, z, y, config, dof_config)` called
   - DOF container initialized based on config
   - Unique ID assigned (class variable incremented)
2. **Element Connection**: `connect_element(element)` called automatically when element created
3. **DOF Prescription**: `prescribe_node()`, `fix_node()`, `free_node()` called
   - Changes propagated to connected elements
4. **Load Application**: `add_load()` called
5. **Analysis**: Node data accessed during DOF assignment and force assembly

## Invariants

- **Connectivity before analysis**: `connectivity()` must be called before any analysis operations
- **DOF index consistency**: All DOF indices must be assigned and consistent across nodes and elements
- **L matrix validity**: L matrix must satisfy `u_all = L @ u_unique` where u_all includes both unique and redundant DOFs
- **B matrix structure**: B matrix rows represent constraints: `B @ u_all = 0` for compatible DOFs
- **Element-node connection**: Elements and nodes maintain bidirectional references
- **DOF value semantics**: `None` = free, `0` = fixed, `float` = prescribed displacement
- **Frequency parameter**: All stiffness/force computations require `omega` parameter
- **Complex values**: All matrices and results are complex-valued (for damping via `E * (1+2j*ksi)`)
- **Section before element type**: Section must be set before element type (element type requires section)
- **Local DOF consistency**: Local DOFs must be subset of global DOFs influenced by element types

## Dependencies

- `numpy`
- `scipy.linalg.inv`
- `collections.defaultdict`
- `pydynsm.elements.ElementFactory`
- `pydynsm.sections.SectionFactory` (via Element)

## File Sizes

- `analysis.py`: ~837 lines (largest file)
- `element.py`: ~1199 lines (largest class file)
- `node.py`: ~253 lines
- `dofs.py`: ~48 lines
