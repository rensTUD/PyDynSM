# Assembler Subsystem

**Location**: `src/pydynsm/assembler.py`

## Responsibilities

The Assembler is the main user-facing API and central orchestrator of PyDynSM. It:

- Manages collections of nodes and elements
- Provides factory methods for creating nodes and elements
- Delegates analysis operations to the Analysis subsystem
- Coordinates visualization through the Plotter subsystem
- Acts as the single entry point for all user interactions

## Public API

### Initialization

```python
Assembler(name: str)
```
Creates a new assembler instance with the given project name. Injects dependencies (Analysis, Node, Element classes, StructurePlotter, ElementFactory).

### Node Management

- `CreateNode(x, z, y=0, config='2D', dof_config=None) -> Node`: Creates and registers a node
- `RegisterNode(node)`: Manually registers an existing node

### Element Management

- `CreateElement(nodes, element_type=None, props={}) -> Element`: Creates and registers an element
- `RegisterElement(element)`: Manually registers an existing element
- `ListElementTypes()`: Lists available element types from ElementFactory

### Analysis Operations

- `GlobalStiffness(omega) -> np.ndarray`: Assembles global unconstrained stiffness matrix
- `GlobalForce(omega) -> np.ndarray`: Assembles global unconstrained force vector
- `GlobalConstrainedStiffness(omega) -> np.ndarray`: Returns constrained stiffness matrix
- `GlobalConstrainedForce(omega) -> np.ndarray`: Returns constrained force vector
- `GlobalConstrainedSystem(omega) -> tuple[np.ndarray, np.ndarray]`: Returns (K_constrained, F_constrained)
- `SolveUfree(Kc_global, fc_global) -> np.ndarray`: Solves for free displacements
- `SupportReactions(k_global, u_free, f_global) -> np.ndarray`: Calculates support reactions
- `FullDisplacement(u_free) -> np.ndarray`: Returns full displacement vector
- `ElementDisplacements(u_full, omega, num_points=20, local_axes=False) -> dict`: Returns element displacements
- `ElementForces(u_full, omega, num_points=20) -> dict`: Returns element forces
- `ElementStresses(u_full, omega, num_points=20) -> dict`: Returns element stresses
- `SolveEigenvector(omega) -> np.ndarray`: Returns eigenvector at specific omega

### Connectivity

- `run_connectivity()`: Computes DOF indices, constraint matrix B, transformation matrix L, and classifies free/constrained DOFs
- `get_dofs_elements() -> tuple`: Returns DOF indices and number of DOFs
- `get_B_matrix() -> tuple`: Returns DOF indices and constraint matrix B

### Visualization

- `PlotStructure(plot_elements=False)`: Plots structure (nodes and optionally elements)
- `PlotElementDisplacements(Element_displacements, scale=1.0)`: Plots element displacements
- `PlotMoments(Element_forces, scale=1.0)`: Plots bending moments
- `PlotAxialforces(Element_forces, scale=1.0)`: Plots axial forces
- `PlotShearforces(Element_forces, scale=1.0)`: Plots shear forces
- `PlotShearstresses(Element_stresses, scale=1.0)`: Plots shear stresses
- `PlotAxialstresses(Element_stresses, scale=1.0)`: Plots axial stresses
- `PlotBendingstresses(Element_stresses, scale=1.0)`: Plots bending stresses

### Persistence (Not Implemented)

- `SaveAssembler()`: Placeholder - not implemented
- `LoadAssembler()`: Placeholder classmethod - not implemented

## Key Data Structures

- `self.name: str`: Project/assembler name
- `self.nodes: list[Node]`: List of registered nodes
- `self.elements: list[Element]`: List of registered elements
- `self.Analysis: Analysis`: Injected Analysis instance
- `self.Node: type`: Injected Node class
- `self.Element: type`: Injected Element class
- `self.StructurePlotter: StructurePlotter`: Injected plotter instance
- `self.ElementFactory: ElementFactory`: Injected element factory

## Lifecycle

1. **Initialization**: User creates `Assembler(name)`
   - Dependencies are injected
   - Empty lists for nodes and elements are created

2. **Model Building**: User creates nodes and elements
   - `CreateNode()` / `CreateElement()` add to internal lists
   - Elements connect to nodes automatically

3. **Configuration**: User sets sections, element types, constraints, loads
   - Done through Element and Node objects directly

4. **Connectivity Analysis**: User calls `run_connectivity()`
   - Must be called before analysis operations
   - Computes DOF indices, constraint matrices, transformation matrices

5. **Analysis**: User calls analysis methods (GlobalStiffness, GlobalForce, etc.)
   - Delegated to `self.Analysis` instance
   - Results returned as NumPy arrays

6. **Post-processing**: User calls ElementDisplacements/Forces/Stresses
   - Delegated to `self.Analysis` instance
   - Results returned as dictionaries

7. **Visualization**: User calls plotting methods
   - Delegated to `self.StructurePlotter` instance

## Invariants

- **Node uniqueness**: Nodes are only added if not already in `self.nodes` list
- **Element uniqueness**: Elements are only added if not already in `self.elements` list
- **Connectivity required**: Analysis operations require `run_connectivity()` to be called first
- **Dependency injection**: All dependencies are injected at initialization, not created internally
- **Frequency parameter**: All analysis methods require `omega` (frequency) parameter
- **Complex results**: All stiffness matrices and results are complex-valued (for damping)

## Dependencies

- `pydynsm.plotter.StructurePlotter`
- `pydynsm.analysis.Node`
- `pydynsm.analysis.Element`
- `pydynsm.analysis.Analysis`
- `pydynsm.elements.ElementFactory`

## File Size

~329 lines
