# Model Creation Workflow

**Purpose**: Build the structural model by creating nodes, elements, sections, element types, constraints, and loads.

## Entry Point

```python
s1 = Assembler('project_name')
```

## Sequence of Execution

### 1. Assembler Initialization

**Entry**: `Assembler.__init__(name)`  
**File**: `src/pydynsm/assembler.py:20-39`

**Call Chain**:
- Creates empty lists: `self.nodes = []`, `self.elements = []`
- Injects dependencies:
  - `self.StructurePlotter = plotter.StructurePlotter()` (line 33)
  - `self.Node = analysis.Node` (line 34)
  - `self.Element = analysis.Element` (line 35)
  - `self.Analysis = analysis.Analysis()` (line 36)
  - `self.ElementFactory = ElementFactory` (line 37)

**Data Created**:
- `Assembler` instance with empty collections

**Caching/Reuse**: None

**Failure Modes**:
- No explicit error handling - Python will raise exceptions for import failures

---

### 2. Create Node

**Entry**: `s1.CreateNode(x, z, y=0, config='2D', dof_config=None)`  
**File**: `src/pydynsm/assembler.py:51-55`

**Call Chain**:
1. `Assembler.CreateNode()` → `self.Node(*args, **kwargs)` (line 53)
2. `Node.__init__(x, z, y, config, dof_config)` → `src/pydynsm/analysis/node.py:55-95`
   - Sets coordinates: `self.x, self.z, self.y` (lines 85-87)
   - Initializes DOF container: `self.dof_container = DOFContainer()` (line 78)
   - Sets up DOFs based on config: loops through `self.dof_config` (lines 81-82)
   - Assigns unique ID: `self.id = Node.nn; Node.nn += 1` (lines 90-91)
   - Initializes empty collections: `self.nodal_loads`, `self.connected_elements` (lines 94-95)
3. `Assembler.RegisterNode(new_node)` → `src/pydynsm/assembler.py:44-49`
   - Adds node to `self.nodes` list if not already present (line 49)

**Data Passed**:
- Input: `x, z, y, config, dof_config`
- Output: `Node` instance with:
  - `node.id`: unique integer ID
  - `node.dof_container`: DOFContainer with DOFs initialized
  - `node.x, node.z, node.y`: coordinates
  - `node.nodal_loads`: defaultdict(dict) for loads
  - `node.connected_elements`: list for connected elements

**Caching/Reuse**: 
- Node IDs are sequential class variables (`Node.nn`)
- Nodes stored in `Assembler.nodes` list

**Failure Modes**:
- Invalid `config` string → KeyError if not in `dof_configurations`
- Invalid `dof_config` → may cause issues in DOF setup
- No explicit error handling for coordinate validation

---

### 3. Create Element

**Entry**: `s1.CreateElement(nodes, element_type=None, props={})`  
**File**: `src/pydynsm/assembler.py:66-74`

**Call Chain**:
1. `Assembler.CreateElement()` → `self.Element(nodes)` (line 68)
2. `Element.__init__(nodes)` → `src/pydynsm/analysis/element.py:71-129`
   - Assigns element ID: `self.id = Element.ne; Element.ne += 1` (lines 82-83)
   - Stores nodes: `self.nodes = nodes` (line 86)
   - Connects element to nodes: `node.connect_element(self)` for each node (lines 89-90)
   - Initializes DOF containers: copies from nodes (lines 92-114)
     - Creates `DOFContainer` for each node (line 103)
     - Copies DOFs from node's `dof_container` to element's `dof_containers[node.id]` (lines 106-109)
     - Sets constraint types to 'monolithic' initially (line 112)
   - Initializes local DOF containers: empty dict per node (line 117)
   - Computes geometrical properties: `self.geometrical_properties()` (line 120)
     - Calculates length `self.L` (line 412-493 in element.py)
     - Computes rotation matrix `self.R` (12x12 matrix, line 493)
   - Initializes collections: `self.element_types = {}`, `self.section = None`, `self.element_loads = defaultdict(dict)` (lines 123-129)
3. `Assembler.RegisterElement(new_element)` → `src/pydynsm/assembler.py:59-64`
   - Adds element to `self.elements` list if not already present (line 64)

**Data Passed**:
- Input: `nodes` (list of Node objects)
- Output: `Element` instance with:
  - `element.id`: unique integer ID
  - `element.nodes`: list of connected nodes
  - `element.dof_containers`: dict mapping node_id to DOFContainer
  - `element.local_dof_container`: dict mapping node_id to DOFContainer (initially empty)
  - `element.constraint_types`: dict mapping node_id to dict of DOF constraint types
  - `element.R`: 12x12 rotation matrix (global to local)
  - `element.L`: element length
  - `element.element_types`: empty dict (to be populated by SetElementType)
  - `element.section`: None (to be set by SetSection)
  - `element.element_loads`: defaultdict(dict) for distributed loads

**Caching/Reuse**:
- Element IDs are sequential class variables (`Element.ne`)
- Elements stored in `Assembler.elements` list
- Rotation matrix `R` computed once and cached
- Length `L` computed once and cached

**Failure Modes**:
- Invalid nodes (not Node instances) → AttributeError when accessing node properties
- Less than 2 nodes → IndexError when accessing nodes[0], nodes[1]
- No explicit error handling for invalid node connections

---

### 4. Set Section

**Entry**: `element.SetSection(section_type, props)`  
**File**: `src/pydynsm/analysis/element.py:496-531`

**Call Chain**:
1. `Element.SetSection(section_type, props)` → line 496
2. `SectionFactory.CreateSection(section_type, **props)` → `src/pydynsm/sections/section.py` (line 523)
   - Factory looks up registered section class
   - Creates section instance with provided dimensions
   - Section computes geometric properties: `A, I_y, I_z, W_y, W_z`
3. Stores section: `self.section = section` (line 526)

**Data Passed**:
- Input: `section_type` (str), `props` (dict of dimensions)
- Output: `Section` instance stored in `element.section` with:
  - `section.A`: cross-sectional area
  - `section.I_y`, `section.I_z`: moments of inertia
  - `section.W_y`, `section.W_z`: section moduli

**Caching/Reuse**:
- Section instance cached in `element.section`
- Section properties computed once during creation

**Failure Modes**:
- Unknown `section_type` → KeyError from SectionFactory
- Missing required props → ValueError from SectionFactory validation
- Exception caught and re-raised with message (lines 529-531)

---

### 5. Set Element Type

**Entry**: `element.SetElementType(element_type, **material_and_element_props)`  
**File**: `src/pydynsm/analysis/element.py:533-587`

**Call Chain**:
1. `Element.SetElementType()` → line 533
2. Validates section exists: checks `self.section is None` (lines 556-557)
   - Raises `ValueError` if section not set
3. `ElementFactory.CreateElement(element_type, section=self.section, L=self.L, **material_and_element_props)` → `src/pydynsm/elements/structuralelement.py` (lines 564-569)
   - Factory looks up registered element class
   - Creates StructuralElement instance with:
     - Section object (extracts A, I_y, I_z, etc.)
     - Element length L
     - Material properties (E, rho, ksi, G, nu)
     - Element-specific properties (kd, cd, T, k, etc.)
4. Adds local DOFs: loops through nodes and element's DOFs (lines 572-576)
   - Adds DOFs to `self.local_dof_container[node.id]` for each DOF in `element.dofs`
5. Applies global constraints to local DOFs: `self.apply_global_constraints_to_local_dofs()` (line 579)
6. Stores element type: `self.element_types[element_type] = element` (line 582)

**Data Passed**:
- Input: `element_type` (str), material properties (E, rho, ksi, etc.), element-specific props
- Output: `StructuralElement` instance stored in `element.element_types[element_type]` with:
  - Access to section properties (A, I_y, I_z, etc.)
  - Material properties (E, rho, ksi, etc.)
  - Element-specific properties
  - `element.dofs`: list of local DOF names (e.g., ['x'] for Rod, ['z', 'phi_y'] for Beam)

**Caching/Reuse**:
- Element type instances cached in `element.element_types` dict
- Multiple element types can be added to same element (composition pattern)

**Failure Modes**:
- Section not set → `ValueError: "Section must be set first using SetSection()"`
- Unknown `element_type` → KeyError from ElementFactory
- Missing required material props → ValueError from ElementFactory validation
- Exception caught and re-raised with message (lines 585-587)

---

### 6. Apply Constraints

Constraints can be applied at two levels: **node-level** (affects all connected elements) or **element-level** (affects only specific element). Additionally, DOFs can be **coupled** (monolithic connection) or **decoupled** (independent movement, creating hinges).

**Note**: Prescribed displacements can be constant values or frequency-dependent functions (lambda functions). For details on using lambda functions, see section 8.

**Constraint Interaction Model**

The constraint system uses a **hierarchical model** with one-way propagation:

```
Node.dof_container (source of truth for node)
    ↓ (propagates via apply_dof_change_to_elements)
Element.dof_containers[node.id] (per-element copy, can be overridden)
    ↓ (used during connectivity analysis)
classify_free_constrained_dofs() reads from element.dof_containers
```

**Key Interaction Rules**:

1. **Node → Element (automatic propagation)**:
   - Node-level constraints automatically propagate to all connected elements via `apply_dof_change_to_elements()` → `element.apply_global_dof_change()` (see `src/pydynsm/analysis/node.py:114-140`)
   - When you set `node.fix_node('x')`, it updates `node.dof_container['x'] = 0` and then calls `apply_dof_change_to_elements()` which updates `element.dof_containers[node.id]['x'] = 0` for all connected elements

2. **Element → Node (no propagation)**:
   - Element-level constraints only affect that specific element and do NOT update the node
   - When you set `element.fix_dof(node, 'x')`, it only updates `element.dof_containers[node.id]['x'] = 0`, leaving `node.dof_container['x']` unchanged

3. **Precedence during analysis**:
   - During connectivity analysis, `classify_free_constrained_dofs()` reads from `element.dof_containers[node_id]` (see `src/pydynsm/analysis/analysis.py:820, 829`)
   - Element-level constraints take precedence for that specific element
   - This allows different constraints per element at the same node

4. **Initial state**:
   - When elements are created, they copy DOFs from nodes (see `src/pydynsm/analysis/element.py:106-109`)
   - Initially, `node.dof_container` and `element.dof_containers[node.id]` match

5. **Divergence allowed**:
   - After creation, node and element DOF values can diverge if element-level constraints are set
   - Node-level constraints will overwrite element-level constraints if set after them

**Example**:
```python
# Node-level constraint (affects all elements at node1)
node1.fix_node('x')  # Sets node.dof_container['x'] = 0, propagates to all elements

# Element-level constraint (overrides for this element only)
element0.fix_dof(node1, 'x')  # Sets element0.dof_containers[node1.id]['x'] = 0
                               # Does NOT change node.dof_container['x']
                               # During analysis, element0 uses its own value
```

**Use Cases**:
- **Node-level**: When all elements at a node should have the same constraint (e.g., all fixed)
- **Element-level**: When different elements at the same node need different constraints (e.g., hinge at one element, fixed at another)

---

#### 6a. Node-Level Constraints

**Entry**: `node.fix_node(*dofs)`, `node.free_node(*dofs)`, `node.prescribe_node(**dofs)`  
**File**: `src/pydynsm/analysis/node.py:154-210`

**Call Chain** (Fix Node):
1. `Node.fix_node(*dofs)` → line 154
   - Calls `self.prescribe_node(**{dof: 0 for dof in dofs})` (line 164)
   - Sets DOF values to `0` in `node.dof_container`
   - Propagates to all connected elements via `apply_dof_change_to_elements()` (line 138)

**Call Chain** (Free Node):
1. `Node.free_node(*dofs)` → line 166
   - Calls `self.prescribe_node(**{dof: None for dof in dofs})` (line 176)
   - Sets DOF values to `None` in `node.dof_container`
   - Propagates to all connected elements

**Call Chain** (Prescribe Node):
1. `Node.prescribe_node(**dofs)` → line 178
   - Validates DOFs exist in `node.dof_container` (line 195)
   - Sets DOF values in `node.dof_container` (line 199)
   - Propagates to all connected elements via `apply_dof_change_to_elements()` (line 138)

**Data Passed**:
- Input: DOF names (e.g., 'x', 'z', 'phi_y') and optional values
- Output: DOF values set in `node.dof_container`:
  - `0` = fixed
  - `None` = free
  - `float` = prescribed displacement (constant)
  - `callable` = prescribed displacement function `f(omega)` (see section 8 for lambda functions)
- Side effects: Changes propagated to all `node.connected_elements` via:
  - `apply_dof_change_to_elements(dof_name, value)` → `src/pydynsm/analysis/node.py:114-140` (line 203)
  - For each connected element: `element.apply_global_dof_change(node, dof_name, value)` → `src/pydynsm/analysis/element.py:1133-1162` (line 138)
  - Updates `element.dof_containers[node.id]` and local DOFs in all connected elements

**Examples**:

```python
# Constant prescribed displacement
node1.prescribe_node(z=0.001)  # 1 mm displacement in z-direction

# Lambda function for frequency-dependent prescribed displacement (see section 8)
omega_target = 100
u_z = lambda omega: 1e-3 if omega == omega_target else 0
node1.prescribe_node(z=u_z)  # 1 mm only at omega=100, 0 otherwise

# More complex frequency-dependent function
def U_omega(omega):
    """Prescribed displacement varies with frequency"""
    if omega < 50:
        return 0.001 * omega / 50
    elif omega < 200:
        return 0.001
    else:
        return 0.001 * (300 - omega) / 100
node2.prescribe_node(x=U_omega)

# Multiple DOFs (mix of constant and callable)
node3.prescribe_node(x=0.0005, z=lambda omega: 0.002 if omega > 150 else 0)
```

**Caching/Reuse**: DOF values stored in `node.dof_container`, used during connectivity analysis

**Interaction with Element-Level Constraints**:
- Node-level constraints **propagate down** to all connected elements
- If an element already has an element-level constraint set, the node-level constraint will **overwrite** it for that element
- This ensures consistency: node-level constraints affect all elements at that node
- See "Constraint Interaction Model" at the top of section 6 for complete details

**Failure Modes**:
- Invalid DOF name → `ValueError` if DOF not in node's configuration (line 196)
- Error handling: Exception caught and printed (line 210)

---

#### 6b. Element-Level Constraints

**Entry**: `element.fix_dof(node, *dofs)`, `element.free_dof(node, *dofs)`, `element.prescribe_dof(node, **dofs)`  
**File**: `src/pydynsm/analysis/element.py:625-714`

**Call Chain** (Fix DOF):
1. `Element.fix_dof(node, *dofs)` → line 625
   - Calls `self.prescribe_dof(node, **{dof_name: 0 for dof_name in dofs})` (line 637)
   - Sets DOF values to `0` in `element.dof_containers[node.id]`

**Call Chain** (Free DOF):
1. `Element.free_dof(node, *dofs)` → line 639
   - Calls `self.prescribe_dof(node, **{dof_name: None for dof_name in dofs})` (line 650)
   - Sets DOF values to `None` in `element.dof_containers[node.id]`

**Call Chain** (Prescribe DOF):
1. `Element.prescribe_dof(node, **dofs)` → line 653
   - Gets element's DOF container: `dof_container = self.dof_containers[node.id]` (line 670)
   - If DOF exists: calls `self._update_dof(node, dof_name, value)` (line 675)
     - Updates global DOF in `element.dof_containers[node.id]` (line 706)
     - Updates local DOFs via `map_global_to_local_dofs()` (lines 709-713)
   - If DOF doesn't exist: checks if influenced by local DOFs (lines 677-686)
   - Updates DOF value in element's DOF container

**Data Passed**:
- Input: `node` (Node object), DOF names and optional values
- Output: DOF values set in `element.dof_containers[node.id]`:
  - `0` = fixed
  - `None` = free
  - `float` = prescribed displacement (constant)
  - `callable` = prescribed displacement function `f(omega)` (see section 8 for lambda functions)
- Side effects: 
  - Local DOFs updated if global DOF changes (via `_update_dof()` → `map_global_to_local_dofs()`)
  - **Does NOT update** `node.dof_container` (element-level constraints are local to element)

**Examples**:

```python
# Constant prescribed displacement
element1.prescribe_dof(node1, z=0.001)  # 1 mm displacement in z-direction

# Lambda function for frequency-dependent prescribed displacement (see section 8)
omega_target = 100
u_z = lambda omega: 1e-3 if omega == omega_target else 0
element1.prescribe_dof(node1, z=u_z)  # 1 mm only at omega=100, 0 otherwise

# Multiple DOFs
element2.prescribe_dof(node2, x=0.0005, z=lambda omega: 0.002 if omega > 150 else 0)
```

**Caching/Reuse**: DOF values stored in `element.dof_containers[node.id]`, used during connectivity analysis

**Interaction with Node-Level Constraints**:
- Element-level constraints **do NOT propagate up** to the node
- Element-level constraints **override** node-level constraints for that specific element
- During connectivity analysis, `classify_free_constrained_dofs()` reads from `element.dof_containers[node_id]` (line 820, 829), so element-level values take precedence
- **Use case**: Allows different constraints per element at the same node (e.g., one element fixed, another free)
- See "Constraint Interaction Model" at the top of section 6 for complete details

**Failure Modes**:
- Invalid DOF name → `ValueError` if DOF not influenced by local DOFs (line 686)
- Error handling: Exception caught and printed (lines 687-688)

---

#### 6c. Couple/Decouple DOFs (Hinges)

**Entry**: `element.couple_dof(node, *dofs)`, `element.decouple_dof(node, *dofs)`  
**File**: `src/pydynsm/analysis/element.py:715-745`

**Call Chain** (Couple DOF):
1. `Element.couple_dof(node, *dofs)` → line 715
   - Sets constraint type to 'monolithic': `self.constraint_types[node.id][dof_name] = "monolithic"` (line 729)
   - Ensures DOF is fully constrained across all elements at this node

**Call Chain** (Decouple DOF):
1. `Element.decouple_dof(node, *dofs)` → line 732
   - Sets constraint type to 'independent': `self.constraint_types[node.id][dof_name] = "independent"` (line 745)
   - Allows DOF to move independently in this element (creates hinge)

**Data Passed**:
- Input: `node` (Node object), DOF names (e.g., 'phi_y' for rotational hinge)
- Output: Constraint type set in `element.constraint_types[node.id][dof_name]`:
  - `'monolithic'` = fully coupled (default)
  - `'independent'` = decoupled (hinge)
- Side effects: Affects constraint matrix B construction during connectivity analysis

**Caching/Reuse**: Constraint types stored in `element.constraint_types`, used by `build_matrix_B()` during connectivity analysis

**Failure Modes**:
- DOF not in element's DOF container → silently ignored (line 728, 744)
- No explicit error handling

**Important Notes**:
- **Decoupling creates hinges**: When a DOF is decoupled, it allows independent movement at that node, creating a hinge connection
- **Constraint matrix B**: Decoupled DOFs are NOT included in constraint equations in `build_matrix_B()` (see `connectivity_analysis.md`)
- **Common use case**: Decouple `'phi_y'` (rotation) to create moment-release hinges in beams
- **After decoupling**: You may need to `free_dof()` the decoupled DOF to avoid singular matrices

---

### 7. Apply Loads

Loads can be applied at two levels: **nodal loads** (point loads at nodes) and **distributed loads** (loads distributed along elements). Both support constant values or frequency-dependent functions (lambda functions).

**Note**: For details on using lambda functions for frequency-dependent loads, see section 8.

**Key Concepts**:
- **Nodal loads**: Applied in **global coordinates** at specific nodes
- **Distributed loads**: Applied in **local coordinates** along elements
- **Load types**: Can be constants or callable functions `f(omega)` for frequency-dependent loading (see section 8 for lambda functions)
- **Evaluation**: Loads are evaluated during analysis when `GlobalForce()` is called with a specific `omega` value

---

#### 7a. Nodal Loads

**Entry**: `node.add_load(**loads)`  
**File**: `src/pydynsm/analysis/node.py:208-224`

**Call Chain**:
1. `Node.add_load(**loads)` → line 208
   - Validates DOF exists in node's configuration (line 219)
   - Stores load in `node.nodal_loads[dof_name] = load` (line 222)
   - Warns if overwriting existing load (lines 220-221)
   - Load can be:
     - **Constant**: A numeric value (e.g., `1000.0`) - applied at all frequencies
     - **Callable**: A function `f(omega)` - evaluated at each frequency during analysis (see section 8 for lambda functions)

**Data Passed**:
- Input: Keyword arguments where key is DOF name (e.g., `'x'`, `'z'`, `'phi_y'`) and value is load (constant or callable)
- Output: Loads stored in `node.nodal_loads[dof_name]` dictionary
- Coordinate system: **Global coordinates** (loads are applied directly to global DOFs)

**Load Evaluation** (during analysis):
- When `Analysis.GlobalForce(nodes, elements, omega)` is called → `src/pydynsm/analysis/analysis.py:61-122`
- For each nodal load (lines 104-111):
  - Gets DOF index from `node.dof_container` (line 108)
  - Evaluates load: `value = load(omega) if callable(load) else load` (line 110)
  - Adds to global force vector: `f_global[index] += value` (line 111)

**Examples**:

```python
# Constant load (applied at all frequencies)
node1.add_load(z=1000.0)  # 1000 N in z-direction

# Lambda function for frequency-dependent load (see section 8 for details)
omega_p = 100  # Target frequency
p_z = lambda omega: 1e6 if omega == omega_p else 0
node1.add_load(z=p_z)  # 1e6 N only at omega=100, 0 otherwise

# More complex frequency-dependent function
def F_omega(omega):
    """Load magnitude varies with frequency"""
    if omega < 50:
        return 1000 * omega / 50
    elif omega < 200:
        return 1000
    else:
        return 1000 * (300 - omega) / 100
node2.add_load(x=F_omega)

# Multiple DOFs (mix of constant and callable)
node3.add_load(x=500.0, z=lambda omega: 2000 if omega > 150 else 0)
```

**Caching/Reuse**: 
- Loads stored in `node.nodal_loads` dictionary
- Evaluated fresh at each `omega` value during analysis
- Same load function can be reused across multiple nodes

**Failure Modes**:
- Invalid DOF name → Warning printed, load not added (line 224)
- DOF not in node's configuration → Load ignored with warning message

---

#### 7b. Distributed Loads

**Note**: Distributed loads support lambda functions for frequency-dependent loading. For details, see section 8.

**Entry**: `element.AddDistributedLoad(**loads)`  
**File**: `src/pydynsm/analysis/element.py:331-344`

**Call Chain**:
1. `Element.AddDistributedLoad(**loads)` → line 331
   - Stores loads in `self.element_loads[dof] = load` (line 344)
   - Warns if overwriting existing load (lines 341-343)
   - Load can be:
     - **Constant**: A numeric value (e.g., `1e3`) - applied at all frequencies
     - **Callable**: A function `f(omega)` - evaluated at each frequency during analysis (see section 8 for lambda functions)

**Data Passed**:
- Input: Keyword arguments where key is DOF name in **local coordinates** (e.g., `'x'` for axial, `'z'` for transverse) and value is load (constant or callable)
- Output: Loads stored in `element.element_loads[dof]` dictionary
- Coordinate system: **Local coordinates** (loads are applied in element's local coordinate system)

**Load Evaluation** (during analysis):
- When `Analysis.GlobalForce(nodes, elements, omega)` is called → `src/pydynsm/analysis/analysis.py:61-122`
- For each element with loads (lines 84-93):
  - Calls `element.EvaluateDistributedLoad(element.element_loads, omega)` → `src/pydynsm/analysis/element.py:362-392`
  - **Step 1 - Load evaluation** (lines 376-384):
    - Loops through element types (line 376)
    - Gets element-specific DOFs (line 378)
    - Evaluates loads for each DOF:
      ```python
      q_evaluated = [
          element_loads.get(dof,0)(omega) if callable(element_loads.get(dof,0)) 
          else element_loads.get(dof,0)
          for dof in dofs
      ]
      ```
    - If load is callable: calls `load(omega)`; if constant: uses value directly
  - **Step 2 - Local force vector** (line 387):
    - Calls `element.FullDistributedLoad(element_type, q_evaluated, omega)` → line 394
    - Which calls `element_type.LocalDistributedLoad(q, omega)` → element-specific implementation (line 405)
    - Element type converts distributed load to equivalent nodal forces in local coordinates
  - **Step 3 - Global transformation** (line 390):
    - Transforms local forces to global: `q_glob = self.R.T @ q_loc`
    - Returns global force vector for element
  - **Step 4 - Assembly** (line 93):
    - Adds element forces to global force vector: `f_global[np.ix_(dofs)] += element_force`

**Examples**:

```python
# Constant distributed load (applied at all frequencies)
element1.AddDistributedLoad(z=1e3)  # 1000 N/m in local z-direction

# Lambda function for frequency-dependent distributed load (see section 8 for details)
omega_p = 100
q_z = lambda omega: 1e6 if omega == omega_p else 0
element2.AddDistributedLoad(z=q_z)  # 1e6 N/m only at omega=100

# Multiple DOFs with different load types
element3.AddDistributedLoad(
    x=500.0,  # Constant axial load
    z=lambda omega: 2000 if omega > 150 else 0  # Frequency-dependent transverse load
)

# Complex frequency-dependent function
def Q_omega(omega):
    """Distributed load varies with frequency"""
    frequency = omega / (2 * np.pi)
    if 10 <= frequency <= 50:
        return 1000 * np.sin(frequency * np.pi / 50)
    else:
        return 0
element4.AddDistributedLoad(x=Q_omega)
```

**Important Notes**:
- **Local coordinates**: Distributed loads are defined in the element's local coordinate system
- **Element-specific DOFs**: Only DOFs relevant to the element type are processed (e.g., `'x'` for Rod, `'z'` and `'phi_y'` for Beam)
- **Multiple element types**: If an element has multiple element types, loads are processed separately for each type
- **Transformation**: Local distributed loads are automatically transformed to global coordinates via rotation matrix `R`

**Caching/Reuse**: 
- Loads stored in `element.element_loads` dictionary
- Evaluated fresh at each `omega` value during analysis
- Same load function can be reused across multiple elements

**Failure Modes**:
- Invalid DOF name → Load stored but may be ignored if DOF not in element type
- No explicit error handling for invalid DOF names
- Missing element type → Load evaluation may fail if element type not set

---

**Summary**:

**Complete Load Application Sequence**:
1. Define load (constant or callable function)
2. Apply to node: `node.add_load(**loads)` (global coordinates)
3. Apply to element: `element.AddDistributedLoad(**loads)` (local coordinates)
4. During analysis: Loads evaluated at each `omega` value
5. Forces assembled into global force vector

**Key Differences**:
- **Nodal loads**: Global coordinates, direct application to DOF
- **Distributed loads**: Local coordinates, converted to equivalent nodal forces by element type
- **Constant loads**: Same value at all frequencies
- **Callable loads**: Evaluated at each frequency during analysis

**Data Flow**:
- Loads stored in: `node.nodal_loads[dof]` or `element.element_loads[dof]`
- Evaluated during: `Analysis.GlobalForce(nodes, elements, omega)`
- For callable loads: `load(omega)` called with current frequency
- Result: Global force vector `F_global` assembled from all loads

---

## 8. Lambda Functions and Frequency-Dependent Values

**Purpose**: Lambda functions (or any callable) allow values to vary with frequency `omega`, enabling frequency-domain analysis with frequency-dependent prescribed displacements and loads.

**How It Works**:

1. **Definition**: Value is defined as a callable function that takes `omega` as parameter:
   ```python
   # Simple lambda
   u_z = lambda omega: 1e-3 if omega == omega_target else 0
   
   # Named function
   def F_omega(omega):
       return 1000 * np.sin(omega / 10)
   ```

2. **Storage**: Function is stored directly (not evaluated):
   ```python
   node.prescribe_node(z=u_z)  # Stores the function, not the result
   # node.dof_container['z'].value = <function u_z>
   
   node.add_load(z=F_omega)  # Stores the function, not the result
   # node.nodal_loads['z'] = <function F_omega>
   ```

3. **Evaluation**: During analysis, function is called with current `omega`:
   ```python
   # For prescribed displacements (in GlobalConstrainedForce)
   value = dof.value(omega) if callable(dof.value) else dof.value
   
   # For loads (in GlobalForce or EvaluateDistributedLoad)
   value = load(omega) if callable(load) else load
   ```

4. **Frequency Sweep**: For each frequency in analysis:
   - `omega = 0` → `u_z(0)` → returns `0`
   - `omega = 100` → `u_z(100)` → returns `1e-3`
   - `omega = 200` → `u_z(200)` → returns `0`

**Common Use Cases**:

```python
# 1. Impulse at specific frequency
omega_target = 100
prescribed_disp = lambda omega: 1e-3 if omega == omega_target else 0
load = lambda omega: 1e6 if omega == omega_target else 0

# 2. Frequency band
def band_value(omega):
    if 50 <= omega <= 150:
        return 1000
    return 0

# 3. Harmonic with phase
def harmonic_value(omega):
    return 1000 * np.exp(1j * omega * 0.1)  # Complex-valued

# 4. Frequency-dependent magnitude
def magnitude_varying(omega):
    return 1000 * (omega / 100) ** 2  # Quadratic variation

# 5. FFT-based (from time-domain signal)
def FFT_value(omega):
    frequency = omega / (2 * np.pi)
    idx = np.where(freq == frequency)[0]
    return P_omega[idx[0]] if len(idx) > 0 else 0
```

**Best Practices**:
- Use lambda for simple conditional values
- Use named functions for complex logic (better debugging)
- Ensure function handles all possible `omega` values (no undefined cases)
- For FFT-based values, handle cases where frequency may not be in FFT bins
- Consider using `np.isclose()` for floating-point comparisons instead of `==`

**Performance Considerations**:
- Lambda functions are evaluated once per frequency per value
- For large frequency sweeps, ensure lambda functions are efficient
- Consider caching results if same function is called multiple times with same `omega`

**Where Lambda Functions Are Used**:
- **Prescribed displacements**: `node.prescribe_node(**dofs)` and `element.prescribe_dof(node, **dofs)` - see section 6
- **Nodal loads**: `node.add_load(**loads)` - see section 7a
- **Distributed loads**: `element.AddDistributedLoad(**loads)` - see section 7b

**Note**: The evaluation mechanism is the same for all uses: if the value is callable, it's called with `omega`; otherwise, the constant value is used directly.

---

## Summary

**Complete Model Creation Sequence**:
1. `Assembler(name)` - Initialize assembler
2. `CreateNode(x, z, ...)` - Create nodes (repeat as needed)
3. `CreateElement([node1, node2])` - Create elements (repeat as needed)
4. `element.SetSection(section_type, props)` - Set cross-section geometry
5. `element.SetElementType(element_type, **props)` - Set element type with material properties
6. Apply constraints:
   - `node.fix_node(*dofs)` - Fix node DOFs (affects all connected elements)
   - `node.free_node(*dofs)` - Free node DOFs
   - `node.prescribe_node(**dofs)` - Prescribe node DOF values
   - `element.fix_dof(node, *dofs)` - Fix element DOFs at specific node
   - `element.free_dof(node, *dofs)` - Free element DOFs at specific node
   - `element.prescribe_dof(node, **dofs)` - Prescribe element DOF values
   - `element.decouple_dof(node, *dofs)` - Decouple DOFs (create hinges)
   - `element.couple_dof(node, *dofs)` - Couple DOFs (ensure monolithic connection)
7. `node.add_load(**loads)` / `element.AddDistributedLoad(**loads)` - Apply loads (repeat as needed)

**Critical Dependencies**:
- Section must be set before element type
- Nodes must exist before creating elements
- Elements must have sections and element types before analysis
- Decouple DOFs before applying other constraints if creating hinges
- After decoupling, may need to `free_dof()` to avoid singular matrices

**Constraint Interaction**: See "Constraint Interaction Model" at the top of section 6 for complete details on how node-level and element-level constraints interact.

**Data Flow**:
- Nodes → stored in `Assembler.nodes`
- Elements → stored in `Assembler.elements`
- Sections → stored in `element.section`
- Element types → stored in `element.element_types`
- Constraints → stored in:
  - DOFContainer values (None/0/float) in `node.dof_container` and `element.dof_containers[node.id]`
  - Constraint types ('monolithic'/'independent') in `element.constraint_types[node.id][dof_name]`
- Loads → stored in `node.nodal_loads` and `element.element_loads`
