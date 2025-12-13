# PyDynSM Architecture Overview

**Purpose**: PyDynSM is a Dynamic Stiffness Matrix (DSM) method solver for structural dynamics analysis. It provides a Python framework for modeling 1D structural elements (rods, beams) and performing frequency-domain dynamic analysis.

## System Overview

PyDynSM follows a layered architecture where:

1. **Assembler** (`src/pydynsm/assembler.py`) - Main user-facing API that orchestrates all subsystems
2. **Analysis** (`src/pydynsm/analysis/`) - Core computational engine for DOF management, matrix assembly, and solving
3. **Elements** (`src/pydynsm/elements/`) - Library of structural element types implementing DSM formulations
4. **Sections** (`src/pydynsm/sections/`) - Cross-sectional geometry definitions
5. **Plotter** (`src/pydynsm/plotter/`) - Visualization of structures and results
6. **Logging** (`src/pydynsm/logging/`) - Logging infrastructure (currently minimal)

## Subsystem Documentation

Detailed documentation for each subsystem:

- [Assembler](subsystems/assembler.md) - Main API and orchestration
- [Analysis](subsystems/analysis.md) - Core computational engine
- [Elements](subsystems/elements.md) - Structural element implementations
- [Sections](subsystems/sections.md) - Cross-sectional geometry
- [Plotter](subsystems/plotter.md) - Visualization
- [Logging](subsystems/logging.md) - Logging infrastructure

## Getting Started

New to the codebase? Start with the [Onboarding Guide](onboarding.md) to learn the key concepts.

## Data Flow

```
User Code
    ↓
Assembler (orchestration)
    ↓
Analysis (DOF assignment, matrix assembly)
    ↓
Elements (stiffness computation)
    ↓
Sections (geometry properties)
```

Results flow back through the same path, with visualization handled by Plotter.

## Key Design Patterns

- **Factory Pattern**: Used for both elements (`ElementFactory`) and sections (`SectionFactory`) to enable registration and creation by name
- **Dependency Injection**: Assembler injects dependencies (Analysis, Node, Element classes) at initialization
- **Frequency-Domain Analysis**: All computations are in frequency domain (omega parameter), results are complex-valued

## Execution Entry Points

- **Library Usage**: Import `pydynsm` and use `Assembler` class (primary entry point)
- **Example Scripts**: Located in `examples/` directory
- **Jupyter Notebooks**: Located in `notebooks/` directory
- **Tests**: Located in `tests/` directory

## Data & Persistence

- **No file I/O**: System is fully programmatic - no configuration files or data persistence
- **In-memory only**: All data structures (nodes, elements, results) exist only in memory
- **NumPy arrays**: All numerical data uses `numpy.ndarray` (often `dtype=complex`)

## Code Quality Notes

- **Largest files**: `analysis/element.py` (~1199 lines), `analysis/analysis.py` (~837 lines)
- **Highest coupling**: Assembler, Analysis, and Element classes are central to the system
- See individual subsystem docs for detailed metrics

---

**Document Version**: 2.0  
**Last Updated**: 2025-01-XX  
**Analyzed Codebase**: PyDynSM v0.2.6
