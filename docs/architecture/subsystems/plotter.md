# Plotter Subsystem

**Location**: `src/pydynsm/plotter/`

## Responsibilities

The Plotter subsystem provides visualization capabilities for:

- Structure visualization (nodes and elements)
- Displacement plots (deformed shape)
- Force diagrams (moments, axial forces, shear forces)
- Stress diagrams (axial, bending, shear stresses)

All plotting is done via matplotlib and displays results (real and imaginary parts) in the frequency domain.

## Public API

### StructurePlotter Class

**Location**: `src/pydynsm/plotter/structureplotter.py`

#### Structure Visualization

- `PlotNodes(nodes)`: Plots node locations as scatter points with labels
- `PlotElements(elements)`: Plots element lines connecting nodes
- `ShowStructure(title='Structure Configuration')`: Displays plot with title, labels, grid, equal aspect ratio

#### Result Visualization (Global)

- `PlotDisplacements(elements, displacements, scale=1.0)`: Plots deformed shape
  - `displacements`: dict mapping `element.id` to `np.ndarray` of shape `(Ndof, num_points)`
  - Shows both real (blue) and imaginary (green) parts

- `Plotmoments(elements, forces, scale=1e-8)`: Plots bending moment diagrams
- `Plotaxialforces(elements, forces, scale=1e-8)`: Plots axial force diagrams
- `Plotshearforces(elements, forces, scale=1e-8)`: Plots shear force diagrams
- `Plotshearstresses(elements, stresses, scale=1e-9)`: Plots shear stress diagrams
- `Plotaxialstresses(elements, stresses, scale=1e-9)`: Plots axial stress diagrams
- `Plotbendingstresses(elements, stresses, scale=1e-9)`: Plots bending stress diagrams

#### Result Visualization (Local/Element-level)

- `plot_element_moment(element, forces)`: Plots moment for single element in local coordinates
- `plot_element_axial_force(element, forces)`: Plots axial force for single element
- `plot_element_shear_force(element, forces)`: Plots shear force for single element
- `plot_element_shear_stress(element, stresses)`: Plots shear stress for single element
- `plot_element_axial_stress(element, stresses)`: Plots axial stress for single element
- `plot_element_bending_stress(element, stresses)`: Plots bending stress for single element

## Key Data Structures

### Input Data Structures

- `nodes: list[Node]`: List of node objects with `x`, `z`, `y` coordinates and `id`
- `elements: list[Element]`: List of element objects with `nodes` and `id`
- `displacements: dict[int, np.ndarray]`: Maps element ID to displacement array of shape `(Ndof, num_points)`
- `forces: dict[int, np.ndarray]`: Maps element ID to force array of shape `(Ndof, num_points)`
- `stresses: dict[int, np.ndarray]`: Maps element ID to stress array of shape `(Ndof, num_points)`

### Plotting Conventions

- **Real part**: Blue solid line (`'b-'`)
- **Imaginary part**: Green dashed line (`'g--'`)
- **Original structure**: Black dashed line (`'k--'`)
- **Nodes**: Red circles (`'red', 'o'`)
- **Displaced nodes**: Blue/green crosses (`'x'`)

## Lifecycle

1. **Initialization**: `StructurePlotter()` created (no parameters, no state)
2. **Structure Plotting**: User calls `PlotNodes()` and/or `PlotElements()`, then `ShowStructure()`
3. **Result Plotting**: User calls result plotting methods with element data dictionaries
4. **Display**: `plt.show()` called automatically in plotting methods

## Invariants

- **No state**: StructurePlotter is stateless (no instance variables)
- **Matplotlib dependency**: All plotting uses `matplotlib.pyplot`
- **Complex data**: All result data is complex-valued (real and imaginary parts plotted separately)
- **Element data format**: Displacement/force/stress data must be dicts mapping element IDs to arrays
- **Array shape**: Result arrays must have shape `(Ndof, num_points)` where Ndof is number of DOFs and num_points is evaluation points
- **Scale factors**: Scale parameters control visualization magnitude (defaults: 1.0 for displacements, 1e-8 for forces, 1e-9 for stresses)
- **Coordinate system**: Plots use global coordinates (x, z) from node positions
- **Perpendicular offset**: Force/stress diagrams are offset perpendicular to element axis for visualization

## Plotting Patterns

### Common Pattern for Force/Stress Plots

1. Extract real and imaginary parts from element data
2. Compute element geometry (start/end coordinates, length, perpendicular direction)
3. Create evaluation points along element (using `np.linspace(0, 1, num_points)`)
4. Offset force/stress values perpendicular to element axis
5. Plot real and imaginary parts
6. Connect diagram back to structure with dashed lines

### Displacement Plot Pattern

1. Extract displacement components (u_x, u_z) from element data
2. Create evaluation points along element
3. Add scaled displacements to original coordinates
4. Plot deformed shape (real and imaginary parts)
5. Mark displaced node positions

## Dependencies

- `matplotlib.pyplot`
- `numpy`

## File Size

- `structureplotter.py`: ~738 lines

## Notes

- All plots are displayed (not saved) - `plt.show()` is called in each method
- Scale factors are provided to control visualization magnitude
- Perpendicular offset direction is computed from element geometry
- Legend entries are added conditionally to avoid duplicates
- Element-level plotting methods use local coordinate system (element along x-axis)
