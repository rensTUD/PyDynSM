# Sections Subsystem

**Location**: `src/pydynsm/sections/`

## Responsibilities

The Sections subsystem provides:

- Abstract base class (`Section`) defining the interface for cross-sectional geometry
- Factory pattern (`SectionFactory`) for section registration and creation by name
- Concrete implementations computing geometric properties from dimensions:
  - Rectangle, Circle, HollowCircle, ISection

Sections compute **only geometric properties** (area, moments of inertia, section moduli). Material properties (E, rho, ksi) are provided separately when setting element types.

## Public API

### Section (Abstract Base Class)

**Location**: `src/pydynsm/sections/section.py`

#### Abstract Methods

- `compute_properties()`: Must compute and set `_A`, `_I_y`, `_I_z`, `_W_y`, `_W_z`

#### Properties (read-only)

- `A: float`: Cross-sectional area [m²]
- `I_y: float`: Second moment of area about y-axis (bending in x-z plane) [m⁴]
- `I_z: float`: Second moment of area about z-axis (bending in x-y plane) [m⁴]
- `W_y: float`: Section modulus about y-axis [m³]
- `W_z: float`: Section modulus about z-axis [m³]

### SectionFactory

**Location**: `src/pydynsm/sections/section.py`

- `RegisterSection(section_class)`: Registers a section class (called automatically via decorator)
- `CreateSection(section_type, **props) -> Section`: Creates section instance by type name
- `ListSectionTypes()`: Prints list of available section types
- `@SectionFactory.SectionType(name)`: Decorator to register section classes

### Section Implementations

All sections are registered via `@SectionFactory.SectionType(name)` decorator and inherit from `Section`.

#### Rectangle (`rectangle.py`, name: `'rectangle'`)

- **Parameters**:
  - `width: float`: Width in y-direction [m]
  - `height: float`: Height in z-direction [m]
- **Properties**:
  - `A = width * height`
  - `I_y = width * height³ / 12`
  - `I_z = height * width³ / 12`
  - `W_y = I_y / (height / 2) = width * height² / 6`
  - `W_z = I_z / (width / 2) = height * width² / 6`

#### Circle (`circle.py`, name: `'circle'`)

- **Parameters**:
  - `diameter: float`: Diameter [m]
- **Properties**:
  - `A = π * (diameter/2)²`
  - `I_y = I_z = π * diameter⁴ / 64` (circular symmetry)
  - `W_y = W_z = I_y / (diameter/2) = π * diameter³ / 32`

#### HollowCircle (`hollow_circle.py`, name: `'hollow_circle'`)

- **Parameters**:
  - `outer_diameter: float`: Outer diameter [m]
  - `inner_diameter: float`: Inner diameter [m]
- **Properties**:
  - `A = π * (outer_diameter² - inner_diameter²) / 4`
  - `I_y = I_z = π * (outer_diameter⁴ - inner_diameter⁴) / 64`
  - `W_y = W_z = I_y / (outer_diameter/2)`

#### ISection (`i_section.py`, name: `'i_section'`)

- **Parameters**:
  - `flange_width: float`: Flange width [m]
  - `flange_thickness: float`: Flange thickness [m]
  - `web_height: float`: Web height [m]
  - `web_thickness: float`: Web thickness [m]
- **Properties**:
  - `A = 2 * flange_width * flange_thickness + web_height * web_thickness`
  - `I_y`, `I_z` computed using parallel axis theorem
  - `W_y`, `W_z` computed from I and extreme fiber distances

## Key Data Structures

### Section Instance Variables

- `self._A: float`: Cross-sectional area (private, accessed via property)
- `self._I_y: float`: Second moment about y-axis (private)
- `self._I_z: float`: Second moment about z-axis (private)
- `self._W_y: float`: Section modulus about y-axis (private)
- `self._W_z: float`: Section modulus about z-axis (private)
- Section-specific dimensions (e.g., `self.width`, `self.height` for Rectangle)

### SectionFactory Class Variables

- `SectionFactory.sections: dict`: Maps section name to dict containing:
  - `'class'`: Section class
  - `'required_params'`: List of required parameter names
  - `'all_params'`: List of all parameter names

## Lifecycle

### Section Registration (at import time)

1. Section class defined with `@SectionFactory.SectionType(name)` decorator
2. Decorator sets `section.section_name = name`
3. `SectionFactory.RegisterSection(section_class)` called
4. Factory inspects `__init__` signature to extract required/optional parameters
5. Section class stored in `SectionFactory.sections` dict

### Section Creation (at runtime)

1. User calls `element.SetSection(section_type, props)`
2. `SectionFactory.CreateSection(section_type, **props)` called
3. Factory validates required parameters are provided
4. Section class instantiated: `section_class(**props)`
5. `__init__` calls `compute_properties()` (via `super().__init__()`)
6. Section instance stored in `element.section`

### Section Usage (during element operations)

1. **Element creation**: Element types access `section.A`, `section.I_y`, etc. during initialization
2. **Stiffness computation**: Element types use section properties in stiffness matrix calculations
3. **Stress computation**: Element types use `section.W_y`, `section.W_z` for stress calculations

## Invariants

- **Section name required**: All section classes must define `section_name` class variable (set by decorator)
- **Abstract method**: `compute_properties()` must be implemented by all concrete sections
- **Property computation**: All properties (`_A`, `_I_y`, `_I_z`, `_W_y`, `_W_z`) must be set in `compute_properties()`
- **Positive dimensions**: All dimension parameters must be positive (sections validate this)
- **Property access**: Properties are read-only (accessed via `@property` decorators)
- **Coordinate system**: 
  - x: horizontal/axial direction
  - z: vertical direction (height)
  - y: in-plane horizontal direction (width)
- **I_y vs I_z**: 
  - `I_y`: Bending in x-z plane (about y-axis)
  - `I_z`: Bending in x-y plane (about z-axis)
- **Section modulus**: `W = I / c` where `c` is distance to extreme fiber
- **Registration uniqueness**: Section names must be unique (factory enforces this)
- **No material properties**: Sections store only geometry, not material properties

## Coordinate System

The coordinate system is consistent across all sections:

- **x-axis**: Horizontal/axial direction (along element length)
- **z-axis**: Vertical direction (typically height dimension)
- **y-axis**: In-plane horizontal direction (typically width dimension)

For bending:
- Bending in **x-z plane** uses `I_y` (rotation about y-axis)
- Bending in **x-y plane** uses `I_z` (rotation about z-axis)

## Dependencies

- `numpy` (for mathematical computations)
- `abc.ABC`, `abc.abstractmethod`
- `inspect` (for parameter validation in factory)

## File Sizes

- `section.py`: ~194 lines
- `i_section.py`: ~130 lines (most complex)
- `rectangle.py`: ~81 lines
- `circle.py`: ~80 lines
- `hollow_circle.py`: ~80 lines
