# Visualization Workflow

**Purpose**: Visualize the structural model and analysis results (displacements, forces, stresses). Visualization uses matplotlib for plotting.

## Entry Point

```python
s1.PlotStructure(plot_elements=False)
s1.PlotElementDisplacements(displacements, scale=1.0)
s1.PlotMoments(forces, scale=1.0)
```

## Sequence of Execution

### 1. Plot Structure

**Entry**: `s1.PlotStructure(plot_elements=False)`  
**File**: `src/pydynsm/assembler.py:80-87`

**Call Chain**:
1. `Assembler.PlotStructure(plot_elements)` → line 80
2. Calls `self.StructurePlotter.PlotNodes(self.nodes)` → `src/pydynsm/plotter/structureplotter.py:18-25` (line 82)
   - Creates figure: `plt.figure(figsize=(10, 6))` (line 19)
   - Loops through nodes (line 20):
     - Plots node: `plt.scatter(node.x, node.z, color='red', marker='o')` (line 24)
     - Adds label: `plt.text(node.x+0.02, node.z+0.02, f'{node.id}')` (line 25)
3. If `plot_elements=True`: calls `self.StructurePlotter.PlotElements(self.elements)` → `src/pydynsm/plotter/structureplotter.py:27-34` (line 85)
   - Loops through elements (line 28):
     - Gets node coordinates: `x_values = [element.nodes[0].x, element.nodes[1].x]` (line 29)
     - Plots element: `plt.plot(x_values, z_values, 'k-', linewidth=2)` (line 31)
     - Adds element ID label at midpoint (lines 32-34)
4. Calls `self.StructurePlotter.ShowStructure(f'Structure Configuration: {self.name}')` → `src/pydynsm/plotter/structureplotter.py` (line 87)
   - Sets title, labels, grid, legend
   - Shows plot: `plt.show()`

**Data Passed**:
- Input: `plot_elements` (bool) - whether to plot elements
- Uses: `self.nodes` (list of Node), `self.elements` (list of Element)
- Output: matplotlib figure displayed

**Caching/Reuse**: None - creates new figure each call

**Failure Modes**:
- Empty nodes/elements lists → empty plot
- Invalid node coordinates → may cause plotting issues
- No explicit error handling

---

### 2. Plot Element Displacements

**Entry**: `s1.PlotElementDisplacements(Element_displacements, scale=1.0)`  
**File**: `src/pydynsm/assembler.py:256-260`

**Call Chain**:
1. `Assembler.PlotElementDisplacements(Element_displacements, scale)` → delegates to `self.StructurePlotter.PlotDisplacements(self.elements, Element_displacements, scale=scale)` (line 260)
2. `StructurePlotter.PlotDisplacements(elements, displacements, scale)` → `src/pydynsm/plotter/structureplotter.py:37-115`
   - Creates figure: `plt.figure(figsize=(10, 8))` (line 52)
   - Sets title, labels, grid (lines 53-56)
   - Plots original structure (nodes and elements) (lines 58-65)
   - Initializes displaced node dictionaries (lines 67-68)
   - Loops through elements (line 71):
     - Gets original node coordinates (lines 73-76)
     - Gets displacement array: `u_elem = displacements[element.id]` (line 79)
     - Extracts real and imaginary parts:
       - `u_x_real = np.real(u_elem[0, :])` (line 80)
       - `u_z_real = np.real(u_elem[1, :])` (line 81)
       - `u_x_imag = np.imag(u_elem[0, :])` (line 82)
       - `u_z_imag = np.imag(u_elem[1, :])` (line 83)
     - Creates points along element: `s = np.linspace(0, 1, num_points)` (line 87)
     - Calculates original positions: `x = x0 + (x1 - x0) * s` (line 90)
     - Adds scaled displacements:
       - `x_disp_real = x + scale * u_x_real` (line 94)
       - `z_disp_real = z + scale * u_z_real` (line 95)
       - `x_disp_imag = x + scale * u_x_imag` (line 96)
       - `z_disp_imag = z + scale * u_z_imag` (line 97)
     - Plots deformed shape (real and imaginary parts) (lines 100-115)
   - Shows plot: `plt.show()`

**Data Passed**:
- Input:
  - `Element_displacements` (dict) - from `ElementDisplacements()` call
  - `scale` (float) - scaling factor for displacements
- Uses: `self.elements` (list of Element)
- Output: matplotlib figure displayed

**Caching/Reuse**: None - creates new figure each call

**Failure Modes**:
- Missing element in displacements dict → `KeyError`
- Invalid displacement array shape → `IndexError`
- No explicit error handling

---

### 3. Plot Moments

**Entry**: `s1.PlotMoments(Element_forces, scale=1.0)`  
**File**: `src/pydynsm/assembler.py:263-267`

**Call Chain**:
1. `Assembler.PlotMoments(Element_forces, scale)` → delegates to `self.StructurePlotter.Plotmoments(self.elements, Element_forces, scale=scale)` (line 267)
2. `StructurePlotter.Plotmoments(elements, Element_forces, scale)` → `src/pydynsm/plotter/structureplotter.py`
   - Similar structure to `PlotDisplacements()`
   - Extracts moment component from force array
   - Plots moment diagram along elements

**Data Passed**:
- Input:
  - `Element_forces` (dict) - from `ElementForces()` call
  - `scale` (float) - scaling factor
- Output: matplotlib figure displayed

**Caching/Reuse**: None

**Failure Modes**: Similar to `PlotDisplacements()`

---

### 4. Plot Axial Forces

**Entry**: `s1.PlotAxialforces(Element_forces, scale=1.0)`  
**File**: `src/pydynsm/assembler.py:269-273`

**Call Chain**:
1. `Assembler.PlotAxialforces(Element_forces, scale)` → delegates to `self.StructurePlotter.Plotaxialforces(self.elements, Element_forces, scale=scale)` (line 273)
2. `StructurePlotter.Plotaxialforces(elements, Element_forces, scale)` → `src/pydynsm/plotter/structureplotter.py`
   - Similar structure to `PlotMoments()`
   - Extracts axial force component (typically Fx)

**Data Passed**: Same as `PlotMoments()`

**Caching/Reuse**: None

**Failure Modes**: Similar to `PlotMoments()`

---

### 5. Plot Shear Forces

**Entry**: `s1.PlotShearforces(Element_forces, scale=1.0)`  
**File**: `src/pydynsm/assembler.py:275-279`

**Call Chain**:
1. `Assembler.PlotShearforces(Element_forces, scale)` → delegates to `self.StructurePlotter.Plotshearforces(self.elements, Element_forces, scale=scale)` (line 279)
2. `StructurePlotter.Plotshearforces(elements, Element_forces, scale)` → `src/pydynsm/plotter/structureplotter.py`
   - Similar structure to `PlotMoments()`
   - Extracts shear force component (typically Fz)

**Data Passed**: Same as `PlotMoments()`

**Caching/Reuse**: None

**Failure Modes**: Similar to `PlotMoments()`

---

### 6. Plot Stresses (Local Element Plots)

**Entry**: `s1.StructurePlotter.plot_element_bending_stress(element, stresses)`  
**File**: `src/pydynsm/plotter/structureplotter.py`

**Call Chain**:
- Similar structure to force plots
- Extracts stress components from stress array
- Plots stress distribution along element

**Data Passed**:
- Input:
  - `element` (Element) - specific element to plot
  - `stresses` (dict) - from `ElementStresses()` call
- Output: matplotlib figure displayed

**Caching/Reuse**: None

**Failure Modes**: Similar to other plotting methods

---

## Summary

**Complete Visualization Sequence**:
1. `s1.PlotStructure(plot_elements=False)` - Plot undeformed structure
2. `s1.PlotElementDisplacements(displacements, scale)` - Plot deformed structure
3. `s1.PlotMoments(forces, scale)` - Plot moment diagrams (optional)
4. `s1.PlotAxialforces(forces, scale)` - Plot axial force diagrams (optional)
5. `s1.PlotShearforces(forces, scale)` - Plot shear force diagrams (optional)
6. `s1.PlotBendingstresses(stresses, scale)` - Plot bending stress (optional)
7. `s1.PlotAxialstresses(stresses, scale)` - Plot axial stress (optional)
8. `s1.PlotShearstresses(stresses, scale)` - Plot shear stress (optional)

**Critical Requirements**:
- For result plots: must have computed displacements/forces/stresses first
- `plot_elements=True` requires elements to be created

**Data Flow**:
- Input: 
  - Structure data: `nodes`, `elements` from Assembler
  - Result data: dicts from `ElementDisplacements()`, `ElementForces()`, `ElementStresses()`
- Processing:
  - Extracts coordinates from nodes/elements
  - Extracts result components from arrays
  - Scales results for visualization
  - Creates matplotlib plots
- Output: matplotlib figures displayed

**Caching/Reuse**:
- No caching of plots - creates new figure each call
- Structure data (nodes, elements) accessed from Assembler
- Result data passed as arguments

**Failure Modes**:
- Missing data in result dicts → `KeyError`
- Invalid array shapes → `IndexError`
- Empty structure → empty plots
- No explicit error handling in most methods
