# How to Read This Codebase

This guide introduces the 5-10 key concepts a new developer must understand before contributing to PyDynSM.

## 1. Dynamic Stiffness Matrix (DSM) Method

**What it is**: PyDynSM implements the Dynamic Stiffness Matrix method for frequency-domain structural dynamics analysis. Unlike finite element methods, DSM uses exact solutions of the differential equations of motion.

**Key insight**: Stiffness matrices are **frequency-dependent** and **complex-valued** (for damping). Every analysis operation requires an `omega` (frequency) parameter.

**Where to see it**: 
- Element implementations in `src/pydynsm/elements/` compute `LocalStiffness(omega)`
- All analysis methods take `omega` as a parameter
- Results are complex-valued NumPy arrays

**Read first**: `src/pydynsm/elements/rod_1D_exp.py` - simplest element implementation

---

## 2. Element Composition Pattern

**What it is**: An `Element` (from `analysis/element.py`) is a **container** that can hold multiple `StructuralElement` instances. This allows combining behaviors (e.g., rod + beam in one element).

**Key insight**: 
- `Element` = geometric entity (nodes, length, rotation, DOFs)
- `StructuralElement` = physical behavior (stiffness computation)
- One `Element` can contain multiple `StructuralElement` instances via `element.element_types` dict

**Where to see it**:
- `Element.SetElementType()` adds a `StructuralElement` to `element.element_types`
- `Element.Stiffness(omega)` accumulates stiffness from all element types
- `src/pydynsm/analysis/element.py` lines 300-326

**Read first**: `src/pydynsm/analysis/element.py` `Stiffness()` method

---

## 3. DOF (Degree of Freedom) System

**What it is**: The DOF system manages degrees of freedom at multiple levels:
- **Node DOFs**: Defined by node configuration (`'2D'`, `'3D'`, etc.)
- **Element DOFs**: Copied from nodes, can be modified per element
- **Local DOFs**: DOFs specific to element types (e.g., rod has only `'x'`, beam has `['z', 'phi_y']`)
- **Global DOF indices**: Assigned during connectivity analysis

**Key insight**: 
- DOF values: `None` = free, `0` = fixed, `float` = prescribed displacement
- Global DOF indices are assigned after `run_connectivity()` is called
- Local DOFs map to global DOFs via rotation matrix `R`

**Where to see it**:
- `src/pydynsm/analysis/dofs.py` - DOF and DOFContainer classes
- `src/pydynsm/analysis/node.py` - Node DOF management
- `src/pydynsm/analysis/element.py` - Element DOF containers
- `src/pydynsm/analysis/analysis.py` - DOF index assignment

**Read first**: `src/pydynsm/analysis/dofs.py`, then `Node.__init__()` in `node.py`

---

## 4. Constraint Matrix (B) and Transformation Matrix (L)

**What it is**: When multiple elements connect at a node, their DOFs must be compatible. The constraint matrix `B` enforces this, and the transformation matrix `L` eliminates redundant DOFs.

**Key insight**:
- **B matrix**: Rows = constraints, Columns = DOFs. `B @ u_all = 0` for compatible DOFs
- **Unique DOFs**: DOFs that appear with positive sign in B (or are unconstrained)
- **Redundant DOFs**: DOFs that appear with negative sign in B (dependent on unique DOFs)
- **L matrix**: Transforms from unique DOF space to full DOF space: `u_all = L @ u_unique`

**Where to see it**:
- `src/pydynsm/analysis/analysis.py`:
  - `build_matrix_B()` - constructs constraint matrix
  - `find_unique_redundant_dofs()` - identifies DOF types
  - `calculate_L()` - computes transformation matrix
  - `GlobalStiffness()` - applies L matrix: `K = L.T @ k_global @ L`

**Read first**: `Analysis.connectivity()` method, then `calculate_L()`

---

## 5. Factory Pattern for Elements and Sections

**What it is**: Both elements and sections use a factory pattern for registration and creation by name.

**Key insight**:
- Classes are registered via decorators: `@ElementFactory.ElementType('Rod')` or `@SectionFactory.SectionType('rectangle')`
- Registration happens at **import time** (side effect of importing the class)
- Creation happens at **runtime** via factory methods: `ElementFactory.CreateElement(name, **kwargs)`
- Factory validates required parameters automatically

**Where to see it**:
- `src/pydynsm/elements/structuralelement.py` - ElementFactory class
- `src/pydynsm/sections/section.py` - SectionFactory class
- `src/pydynsm/elements/rod_1D_exp.py` - Example element registration
- `src/pydynsm/sections/rectangle.py` - Example section registration

**Read first**: `ElementFactory.RegisterElement()` and `ElementFactory.CreateElement()` methods

---

## 6. Coordinate Transformations (Global ↔ Local)

**What it is**: Elements have local coordinate systems (along element axis), but the global system uses arbitrary orientations. A 12x12 rotation matrix `R` transforms between them.

**Key insight**:
- **Local coordinates**: Element lies along x-axis, z is perpendicular
- **Global coordinates**: Element can be at any angle in 3D space
- **Rotation matrix R**: 12x12 matrix (6 DOFs per node × 2 nodes) transforms global to local
- Stiffness: `K_global = R.T @ K_local @ R`
- Displacements: `u_local = R @ u_global`

**Where to see it**:
- `src/pydynsm/analysis/element.py`:
  - `geometrical_properties()` - computes R matrix (lines 412-493)
  - `Stiffness()` - applies rotation (line 321)
  - `Displacements()` - transforms to local, then back to global (lines 799-892)

**Read first**: `Element.geometrical_properties()` method

---

## 7. Frequency-Domain Analysis

**What it is**: All analysis is performed in the frequency domain. Results are complex-valued (real part = in-phase, imaginary part = out-of-phase with excitation).

**Key insight**:
- Every computation requires `omega` (angular frequency) parameter
- Material damping included via complex modulus: `E_complex = E * (1 + 2j * ksi)`
- Results are complex: `u = u_real + j * u_imag`
- Time-domain response requires inverse Fourier transform (not implemented in core)

**Where to see it**:
- All `LocalStiffness(omega)` methods in element implementations
- `Analysis.GlobalStiffness(nodes, elements, omega)`
- Element wavenumber computation depends on omega

**Read first**: Any element's `ElementWaveNumbers(omega)` method

---

## 8. Section vs Element Type Separation

**What it is**: Cross-sectional geometry (sections) is separated from material properties and element behavior (element types).

**Key insight**:
- **Section**: Geometry only (A, I_y, I_z, W_y, W_z) - no material properties
- **Element Type**: Material properties (E, rho, ksi) + behavior (DSM formulation)
- **Two-step setup**: First `SetSection()`, then `SetElementType()`
- Element type receives section object and extracts geometric properties

**Where to see it**:
- `src/pydynsm/analysis/element.py`:
  - `SetSection()` - creates section via SectionFactory
  - `SetElementType()` - creates element type via ElementFactory, passes section
- `src/pydynsm/elements/rod_1D_exp.py` - element type extracts `section.A` in `__init__`

**Read first**: `Element.SetSection()` and `Element.SetElementType()` methods

---

## 9. Assembler as Orchestrator

**What it is**: The `Assembler` class is the main user-facing API that orchestrates all subsystems. It uses dependency injection to avoid tight coupling.

**Key insight**:
- Assembler doesn't create classes directly - it receives them as dependencies
- All user interactions go through Assembler
- Assembler delegates to Analysis, Plotter, ElementFactory
- Assembler manages collections (nodes, elements lists)

**Where to see it**:
- `src/pydynsm/assembler.py` - entire file
- `__init__()` - dependency injection (lines 20-39)
- All methods delegate to injected dependencies

**Read first**: `Assembler.__init__()` to see dependency injection pattern

---

## 10. Analysis Workflow

**What it is**: The typical analysis workflow follows a specific sequence that must be understood.

**Key insight - Required sequence**:
1. Create Assembler
2. Create nodes and elements
3. Set sections and element types
4. Apply constraints and loads
5. **Call `run_connectivity()`** - **CRITICAL STEP**
6. Assemble global matrices
7. Solve system
8. Post-process results

**Where to see it**:
- `examples/Examples_functionality.py` - complete workflow example
- `src/pydynsm/assembler.py` - all workflow methods
- `src/pydynsm/analysis/analysis.py` - `connectivity()` method

**Read first**: Any example in `examples/` directory to see full workflow

---

## Recommended Reading Order

1. **Start here**: `examples/main_notebook3.1.py` - simple example showing workflow
2. **Understand DOFs**: `src/pydynsm/analysis/dofs.py` and `node.py`
3. **Understand Elements**: `src/pydynsm/elements/rod_1D_exp.py` (simplest element)
4. **Understand Analysis**: `src/pydynsm/analysis/analysis.py` `connectivity()` method
5. **Understand Assembler**: `src/pydynsm/assembler.py` - see how it orchestrates
6. **Dive deeper**: Read subsystem documentation in `docs/architecture/subsystems/`

## Common Pitfalls

- **Forgetting `run_connectivity()`**: Analysis will fail if connectivity not computed first
- **Mixing local/global coordinates**: Be careful which coordinate system you're working in
- **Frequency parameter**: All analysis requires `omega` - don't forget it
- **Complex values**: Results are complex - use `np.real()` and `np.imag()` for visualization
- **Section before element type**: Must set section before element type (element type needs section properties)

## Key Files to Understand

- `src/pydynsm/assembler.py` - Main API
- `src/pydynsm/analysis/analysis.py` - Core computations
- `src/pydynsm/analysis/element.py` - Element container
- `src/pydynsm/elements/structuralelement.py` - Element base class
- `src/pydynsm/elements/rod_1D_exp.py` - Simplest element implementation
- `src/pydynsm/sections/section.py` - Section base class
