# Elements Subsystem

**Location**: `src/pydynsm/elements/`

## Responsibilities

The Elements subsystem provides:

- Abstract base class (`StructuralElement`) defining the interface for all element types
- Factory pattern (`ElementFactory`) for element registration and creation by name
- Concrete implementations of structural elements implementing DSM (Dynamic Stiffness Matrix) formulations:
  - Rod elements (1D, foundation, Rayleigh-Bishop, Rayleigh-Love)
  - Beam elements (Euler-Bernoulli, with foundation, tensioned, shear)

## Public API

### StructuralElement (Abstract Base Class)

**Location**: `src/pydynsm/elements/structuralelement.py`

#### Abstract Methods (must be implemented by subclasses)

- `ElementWaveNumbers(omega) -> np.ndarray`: Returns wavenumbers for the element at frequency omega
- `LocalStiffness(omega) -> np.ndarray`: Returns local stiffness matrix (2x2 for rods, 4x4 for beams typically)
- `LocalElementDisplacements(u_nodes_local, omega, num_points) -> list[np.ndarray]`: Returns displacement field along element
- `LocalElementForces(u_nodes_local, omega, num_points) -> list[np.ndarray]`: Returns force field along element
- `LocalElementStresses(u_nodes_local, omega, num_points) -> list[np.ndarray]`: Returns stress field along element
- `Coefficients(u_nodes_local, omega) -> np.ndarray`: Returns solution coefficients from boundary conditions
- `LocalDistributedLoad(q, omega) -> np.ndarray`: Converts distributed load q to nodal load vector

#### Helper Methods (implemented in base class)

- `GetNodalDofs(dofs, Ndof) -> list[int]`: Maps element DOFs to nodal DOF indices
- `FullStiffness(omega) -> np.ndarray`: Expands local stiffness to full 6x6 or 12x12 matrix
- `FullDistributedLoad(q, omega) -> np.ndarray`: Expands local load to full vector
- `FullDisplacement(u_node, omega) -> np.ndarray`: Expands local displacement to full vector
- `FullElementDisplacements(u_nodes_global, omega, num_points) -> list`: Expands element displacements to full DOF space

### ElementFactory

**Location**: `src/pydynsm/elements/structuralelement.py`

- `RegisterElement(element_class)`: Registers an element class (called automatically via decorator)
- `CreateElement(element_name, **kwargs) -> StructuralElement`: Creates element instance by name
- `ListElementTypes()`: Prints list of available element types
- `@ElementFactory.ElementType(name)`: Decorator to register element classes

### Element Implementations

All elements are registered via `@ElementFactory.ElementType(name)` decorator and inherit from `StructuralElement`.

#### Rod Elements

- **`Rod1D`** (`rod_1D_exp.py`, name: `'Rod'`): Basic 1D rod element
  - DOFs: `['x']` (axial)
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`

- **`Rod1D_foundation`** (`rod_1D_foundation.py`, name: `'Rod Foundation'`): Rod on elastic foundation
  - DOFs: `['x']`
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`, `k` (foundation stiffness)

- **`RayleighBishopRod`** (`rb_rod.py`): Rayleigh-Bishop rod (includes rotational inertia)
  - DOFs: `['x']`
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`

- **`RayleighLoveRod`** (`rl_rod.py`): Rayleigh-Love rod (includes lateral inertia)
  - DOFs: `['x']`
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`

#### Beam Elements

- **`EulerBernoulliBeam`** (`eb_beam_exp.py`, name: `'EulerBernoulli Beam'`): Classic Euler-Bernoulli beam
  - DOFs: `['z', 'phi_y']` (transverse displacement and rotation)
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`

- **`EulerBernoulliBeamFoundation`** (`eb_beam_foundation.py`, name: `'EulerBernoulli Beam Foundation'`): EB beam on foundation
  - DOFs: `['z', 'phi_y']`
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`, `kd` (distributed foundation stiffness), `cd` (distributed damping)

- **`TensionedEulerBernoulliBeam`** (`eb_beam_tensioned.py`, name: `'Tensioned EulerBernoulli Beam'`): EB beam with axial tension
  - DOFs: `['z', 'phi_y']`
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`, `T` (tension force)

- **`EulerBernoulliBeamFoundationEndAttachment`** (`eb_beam_foundation_end_attachment.py`): EB beam with foundation and end attachments
  - DOFs: `['z', 'phi_y']`
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`, `kd`, `cd`, plus end attachment parameters

- **`RayleighBeam`** (`rayleigh_beam.py`): Rayleigh beam (includes rotational inertia)
  - DOFs: `['z', 'phi_y']`
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`

- **`Shear1D`** (`shear_beam_1D_exp.py`, name: `'Shear Beam'`): Shear beam element
  - DOFs: `['z']` (transverse only, no rotation)
  - Parameters: `section`, `L`, `E`, `rho`, `ksi`, `G` (shear modulus)

- **Timoshenko Beam** (`timoshenko_beam.py`): Noted as incomplete in code comments

## Key Data Structures

### StructuralElement Instance Variables

- `self.dofs: list[str]`: List of DOF names this element contributes to (e.g., `['x']` for rod, `['z', 'phi_y']` for beam)
- `self.element_name: str`: Name of element type (class variable, set by decorator)
- `self.section: Section`: Reference to section object (provides A, I_y, I_z, etc.)
- `self.L: float`: Element length
- `self.E: complex`: Complex modulus (includes damping: `E * (1+2j*ksi)`)
- `self.rho: float`: Material density
- `self.ksi: float`: Damping ratio
- Element-specific properties (e.g., `k` for foundation, `T` for tension, `G` for shear)

### ElementFactory Class Variables

- `ElementFactory.elements: dict`: Maps element name to dict containing:
  - `'class'`: Element class
  - `'required_params'`: List of required parameter names
  - `'all_params'`: List of all parameter names

## Lifecycle

### Element Registration (at import time)

1. Element class defined with `@ElementFactory.ElementType(name)` decorator
2. Decorator sets `element.element_name = name`
3. `ElementFactory.RegisterElement(element_class)` called
4. Factory inspects `__init__` signature to extract required/optional parameters
5. Element class stored in `ElementFactory.elements` dict

### Element Creation (at runtime)

1. User calls `Element.SetElementType(element_type, **props)`
2. `ElementFactory.CreateElement(element_type, section=section, L=L, **props)` called
3. Factory validates required parameters are provided
4. Element class instantiated: `element_class(section=section, L=L, **props)`
5. Element instance stored in `element.element_types[element_type]` dict

### Element Usage (during analysis)

1. **Stiffness computation**: `element.Stiffness(omega)` calls `element_type.LocalStiffness(omega)` for each element type
2. **Load evaluation**: `element.EvaluateDistributedLoad()` calls `element_type.LocalDistributedLoad(q, omega)`
3. **Post-processing**: `element.Displacements()` calls `element_type.LocalElementDisplacements(u_local, omega, num_points)`

## Invariants

- **Element name required**: All element classes must define `element_name` class variable (set by decorator)
- **Abstract methods**: All abstract methods must be implemented by concrete element classes
- **Section dependency**: Element creation requires `section` parameter (Section object)
- **Length dependency**: Element creation requires `L` parameter (element length)
- **Material properties**: Elements require `E`, `rho`, `ksi` at minimum
- **DOF consistency**: Element `dofs` list must match DOFs used in LocalStiffness matrix dimensions
- **Frequency parameter**: All computations require `omega` (frequency) parameter
- **Complex stiffness**: LocalStiffness returns complex-valued matrix (for damping)
- **Wavenumber computation**: ElementWaveNumbers must return valid wavenumbers for given omega
- **Registration uniqueness**: Element names must be unique (factory enforces this)

## Common Implementation Pattern

All element implementations follow this pattern:

1. **Initialization**: Extract section properties, store material properties, compute derived quantities
2. **ElementWaveNumbers**: Solve characteristic equation to get wavenumbers (complex, frequency-dependent)
3. **LocalStiffness**: Build 2x2 or 4x4 stiffness matrix using wavenumbers and boundary conditions
4. **Coefficients**: Solve for solution coefficients from nodal displacements
5. **LocalElementDisplacements**: Evaluate displacement field using coefficients and wavenumbers
6. **LocalElementForces**: Compute forces from displacement derivatives
7. **LocalElementStresses**: Compute stresses from forces and section properties
8. **LocalDistributedLoad**: Convert distributed load to equivalent nodal loads

## Dependencies

- `numpy`
- `abc.ABC`, `abc.abstractmethod`
- `inspect` (for parameter validation in factory)
- `pydynsm.sections.Section` (elements require section objects)

## File Sizes

- `structuralelement.py`: ~265 lines
- `eb_beam_foundation_end_attachment.py`: ~544 lines (largest element)
- `eb_beam_foundation.py`: ~492 lines
- `shear_beam_1D_exp.py`: ~342 lines
- `rod_1D_exp.py`: ~342 lines
- `eb_beam_tensioned.py`: ~284 lines
- Other elements: ~200-350 lines each
