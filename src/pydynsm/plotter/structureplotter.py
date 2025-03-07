# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 19:52:53 2024

@author: rensv
"""

# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np

# %% Class definition
class StructurePlotter:
    def __init__(self):
        pass
    
    # %% PLOT STRUCTURE GLOBAL
    def PlotNodes(self, nodes):
        plt.figure(figsize=(10, 6))
        for node in nodes:
            # marker = 'o' if len(node.constrained_dofs) == 0 else '^'
            # color = 'blue' if len(node.constrained_dofs) == 0 else 'red'
            # plt.scatter(node.x, node.z, color=color, marker=marker, label=f'{"Free" if len(node.constrained_dofs) == 0 else "Fixed"} Node' if f'{"Free" if len(node.constrained_dofs) == 0 else "Fixed"} Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.scatter(node.x, node.z, color='red', marker='o', label=f'Node: {node.id}')
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')

    def PlotElements(self, elements):
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k-', linewidth=2, label=f'Element: {element.id}' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")
            mid_x = sum(x_values) / 2
            mid_z = sum(z_values) / 2
            plt.text(mid_x, mid_z+0.02, element.id, fontsize=9, color='black')
    
    # %% PLOT RESULTS GLOBAL
    def PlotDisplacements(self, elements, displacements, scale=1.0):
        """
        Plots the displaced shape of the structure based on the provided displacements.

        Parameters
        ----------
        elements : list of Element
            The list of elements in the structural system.
        displacements : list of numpy.ndarray
            List of displacements for each element as returned by `Analysis.ElementDisplacements`.
        scale : float
            Scale factor for displacements to visualize the deformations clearly.
        num_points : int
            Number of points along each element length for plotting the displaced shape.
        """
        plt.figure(figsize=(10, 8))
        plt.title('Real and Imaginary Parts of Displaced Structure')
        plt.xlabel('X Coordinate')
        plt.ylabel('Z Coordinate')
        plt.grid(True)

        # Plot the original structure (nodes and elements)
        for node in [node for element in elements for node in element.nodes]:
            plt.scatter(node.x, node.z, color='red', marker='o', label='Node' if 'Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k--', linewidth=2, label='Element' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")

        displaced_nodes_real = {}
        displaced_nodes_imag = {}

        # Plot the displaced elements (real and imaginary parts)
        for element in elements:
            # Extract the original coordinates of the nodes
            x0 = element.nodes[0].x
            z0 = element.nodes[0].z
            x1 = element.nodes[1].x
            z1 = element.nodes[1].z

            # Get the displacement components (assuming u_elem has shape (6, num_points))
            u_elem = displacements[element.id]
            u_x_real = np.real(u_elem[0, :])
            u_z_real = np.real(u_elem[1, :])
            u_x_imag = np.imag(u_elem[0, :])
            u_z_imag = np.imag(u_elem[1, :])
            
            num_points = len(u_x_real)
            # Create an array of points along the length of the element
            s = np.linspace(0, 1, num_points)

            # Calculate the original positions along the element's length
            x = x0 + (x1 - x0) * s
            z = z0 + (z1 - z0) * s

            # Add the displacements scaled by the specified factor
            x_disp_real = x + scale * u_x_real
            z_disp_real = z + scale * u_z_real
            x_disp_imag = x + scale * u_x_imag
            z_disp_imag = z + scale * u_z_imag

            # Plot the deformed shape of the element (real and imaginary parts)
            plt.plot(x_disp_real, z_disp_real, 'b-', linewidth=1.5, label='Real Part of Displaced Element' if 'Real Part of Displaced Element' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.plot(x_disp_imag, z_disp_imag, 'g--', linewidth=1.5, label='Imaginary Part of Displaced Element' if 'Imaginary Part of Displaced Element' not in plt.gca().get_legend_handles_labels()[1] else "")

            # Save the displaced coordinates of the nodes (start and end points) for both real and imaginary parts
            if element.nodes[0].id not in displaced_nodes_real:
                displaced_nodes_real[element.nodes[0].id] = (x_disp_real[0], z_disp_real[0])
            if element.nodes[1].id not in displaced_nodes_real:
                displaced_nodes_real[element.nodes[1].id] = (x_disp_real[-1], z_disp_real[-1])

            if element.nodes[0].id not in displaced_nodes_imag:
                displaced_nodes_imag[element.nodes[0].id] = (x_disp_imag[0], z_disp_imag[0])
            if element.nodes[1].id not in displaced_nodes_imag:
                displaced_nodes_imag[element.nodes[1].id] = (x_disp_imag[-1], z_disp_imag[-1])

        # Plot the displaced nodes (real and imaginary parts)
        for node_id, (x_disp_real, z_disp_real) in displaced_nodes_real.items():
            plt.scatter(x_disp_real, z_disp_real, color='blue', marker='x', label='Real Part of Displaced Node' if 'Real Part of Displaced Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(x_disp_real+0.02, z_disp_real+0.02, f'{node_id}', fontsize=9, ha='right', color='blue')

        for node_id, (x_disp_imag, z_disp_imag) in displaced_nodes_imag.items():
            plt.scatter(x_disp_imag, z_disp_imag, color='green', marker='x', label='Imaginary Part of Displaced Node' if 'Imaginary Part of Displaced Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(x_disp_imag+0.02, z_disp_imag+0.02, f'{node_id}', fontsize=9, ha='right', color='green')

        # Add legend and show the plot
        plt.legend()
        plt.show()

    def Plotmoments(self, elements, forces, scale=1e-8): # TODO: check the direction of moment, shear force lines
        
        plt.figure(figsize=(10, 8))
        plt.title('Bending moments plot')
        plt.grid(True)

        # Plot the original structure (nodes and elements)
        for node in [node for element in elements for node in element.nodes]:
            plt.scatter(node.x, node.z, color='red', marker='o', label='Node' if 'Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k--', linewidth=2, label='Element' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")

        mmt_nodes_real = {}
        mmt_nodes_imag = {}

        # Plot the displaced elements (real and imaginary parts)
        for element in elements:
            # Extract the original coordinates of the nodes
            x0 = element.nodes[0].x
            z0 = element.nodes[0].z
            x1 = element.nodes[1].x
            z1 = element.nodes[1].z
            
            # Get the displacement components (assuming u_elem has shape (6, num_points))
            f_elem = forces[element.id]
            m_real = np.real(f_elem[-1, :])
            m_imag = np.imag(f_elem[-1, :])
            
            # m_real = m_real/np.max(np.abs(m_real)) # need to be removed after the development of scaling factor
            # m_imag = m_imag
            
            num_points = len(m_real)
            # Create an array of points along the length of the element
            s = np.linspace(0, 1, num_points)

            # Calculate the original positions along the element's length
            x = x0 + (x1 - x0) * s
            z = z0 + (z1 - z0) * s
          # Compute perpendicular direction for offsetting the moment line
            dx = x1 - x0
            dz = z1 - z0
            length = np.sqrt(dx**2 + dz**2)
            perp_x = -dz / length  # Unit normal in X
            perp_z = dx / length   # Unit normal in Z          
            
            moment_x = x + scale * m_real * perp_x
            moment_z = z + scale * m_real * perp_z
            
            moment_x_imag = x + scale * m_imag * perp_x
            moment_z_imag = z + scale * m_imag * perp_z
            
            # Plot the moment line
            plt.plot(moment_x, moment_z, 'b-', linewidth=1, label='Real Part of Bending Moment' if element.id == elements[0].id else "")
            plt.plot(moment_x_imag, moment_z_imag, 'g--', linewidth=1, label='Imaginary Part of Bending Moment' if element.id == elements[0].id else "")
            
            # Connect moment diagram back to the structure
            plt.plot([x[0], moment_x[0]], [z[0], moment_z[0]], 'b--', linewidth=1)
            plt.plot([x[-1], moment_x[-1]], [z[-1], moment_z[-1]], 'b--', linewidth=1)
            plt.plot([x[0], moment_x_imag[0]], [z[0], moment_z_imag[0]], 'g--', linewidth=1)
            plt.plot([x[-1], moment_x_imag[-1]], [z[-1], moment_z_imag[-1]], 'g--', linewidth=1)
        # Add legend and show the plot
        plt.legend()
        plt.show()    

    def Plotaxialforces(self, elements, forces, scale=1e-8):
        
        plt.figure(figsize=(10, 8))
        plt.title('Axial Force Plot')
        plt.grid(True)
  
        # Plot the original structure (nodes and elements)
        for node in [node for element in elements for node in element.nodes]:
            plt.scatter(node.x, node.z, color='red', marker='o', label='Node' if 'Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k--', linewidth=2, label='Element' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")
  
        mmt_nodes_real = {}
        mmt_nodes_imag = {}
  
        # Plot the displaced elements (real and imaginary parts)
        for element in elements:
            # Extract the original coordinates of the nodes
            x0 = element.nodes[0].x
            z0 = element.nodes[0].z
            x1 = element.nodes[1].x
            z1 = element.nodes[1].z
            
            # Get the displacement components (assuming u_elem has shape (6, num_points))
            f_elem = forces[element.id]
            N_real = np.real(f_elem[0, :])
            N_imag = np.imag(f_elem[0, :])
            
            # N_real = N_real/np.max(np.abs(N_real))
            
            num_points = len(N_real)
            # Create an array of points along the length of the element
            s = np.linspace(0, 1, num_points)
  
            # Calculate the original positions along the element's length
            x = x0 + (x1 - x0) * s
            z = z0 + (z1 - z0) * s
          # Compute perpendicular direction for offsetting the moment line
            dx = x1 - x0
            dz = z1 - z0
            length = np.sqrt(dx**2 + dz**2)
            perp_x = dz / length  # Unit normal in X
            perp_z = dx / length   # Unit normal in Z          
            
            axialforce_x = x + scale * N_real * perp_x
            axialforce_z = z - scale * N_real * perp_z
            
            axialforce_x_imag = x + scale * N_imag * perp_x
            axialforce_z_imag = z - scale * N_imag * perp_x
            
            # Plot the moment line
            plt.plot(axialforce_x, axialforce_z, 'b-', linewidth=1, label='Real Part of Axial Force' if element.id == elements[0].id else "")
            plt.plot(axialforce_x_imag, axialforce_z_imag, 'g--', linewidth=1, label='Imaginary Part of Axial Force' if element.id == elements[0].id else "")
        
            # Connect moment diagram back to the structure
            plt.plot([x[0], axialforce_x[0]], [z[0], axialforce_z[0]], 'b--', linewidth=1)
            plt.plot([x[-1], axialforce_x[-1]], [z[-1], axialforce_z[-1]], 'b--', linewidth=1)
            plt.plot([x[0], axialforce_x_imag[0]], [z[0], axialforce_z_imag[0]], 'g--', linewidth=1)
            plt.plot([x[-1], axialforce_x_imag[-1]], [z[-1], axialforce_z_imag[-1]], 'g--', linewidth=1)
  
        # Add legend and show the plot
        plt.legend()
        plt.show()
        
    def Plotshearforces(self, elements, forces, scale=1e-8):
        
        plt.figure(figsize=(10, 8))
        plt.title('Shear Force Plot')
        plt.grid(True)
  
        # Plot the original structure (nodes and elements)
        for node in [node for element in elements for node in element.nodes]:
            plt.scatter(node.x, node.z, color='red', marker='o', label='Node' if 'Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k--', linewidth=2, label='Element' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")
  
        sf_nodes_real = {}
        sf_nodes_imag = {}
  
        # Plot the displaced elements (real and imaginary parts)
        for element in elements:
            # Extract the original coordinates of the nodes
            x0 = element.nodes[0].x
            z0 = element.nodes[0].z
            x1 = element.nodes[1].x
            z1 = element.nodes[1].z
            
            # Get the displacement components (assuming u_elem has shape (6, num_points))
            f_elem = forces[element.id]
            V_real = np.real(f_elem[1, :])
            V_imag = np.imag(f_elem[1, :])
            
            num_points = len(V_real)
            # Create an array of points along the length of the element
            s = np.linspace(0, 1, num_points)
      
            # Calculate the original positions along the element's length
            x = x0 + (x1 - x0) * s
            z = z0 + (z1 - z0) * s
          # Compute perpendicular direction for offsetting the moment line
            dx = x1 - x0
            dz = z1 - z0
            length = np.sqrt(dx**2 + dz**2)
            perp_x = dz / length  # Unit normal in X
            perp_z = dx / length   # Unit normal in Z          
            
            shearforce_x = x - scale * V_real * perp_x
            shearforce_z = z + scale * V_real * perp_z
            
            shearforce_x_imag = x - scale * V_imag * perp_x
            shearforce_z_imag = z + scale * V_imag * perp_z
            
            # Plot the moment line
            plt.plot(shearforce_x, shearforce_z, 'b-', linewidth=1, label='Real Part of Shear Force' if element.id == elements[0].id else "")
            plt.plot(shearforce_x_imag, shearforce_z_imag, 'g--', linewidth=1, label='Imaginary Part of Shear Force' if element.id == elements[0].id else "")
            
            # Connect moment diagram back to the structure
            plt.plot([x[0], shearforce_x[0]], [z[0], shearforce_z[0]], 'b--', linewidth=1)
            plt.plot([x[-1], shearforce_x[-1]], [z[-1], shearforce_z[-1]], 'b--', linewidth=1)
            plt.plot([x[0], shearforce_x_imag[0]], [z[0], shearforce_z_imag[0]], 'g--', linewidth=1)
            plt.plot([x[-1], shearforce_z_imag[-1]], [z[-1], shearforce_z_imag[-1]], 'g--', linewidth=1)  
        # Add legend and show the plot
        plt.legend()
        plt.show()
 
    def Plotshearstresses(self, elements, stresses, scale=1e-9):
     
        plt.figure(figsize=(10, 8))
        plt.title('Shear Stress Plot')
        plt.grid(True)
    
        # Plot the original structure (nodes and elements)
        for node in [node for element in elements for node in element.nodes]:
            plt.scatter(node.x, node.z, color='red', marker='o', label='Node' if 'Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k--', linewidth=2, label='Element' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")
    
        ss_nodes_real = {}
        ss_nodes_imag = {}
    
        # Plot the displaced elements (real and imaginary parts)
        for element in elements:
            # Extract the original coordinates of the nodes
            x0 = element.nodes[0].x
            z0 = element.nodes[0].z
            x1 = element.nodes[1].x
            z1 = element.nodes[1].z
            
            # Get the displacement components (assuming u_elem has shape (6, num_points))
            sigma_elem = stresses[element.id]
            tau_real = np.real(sigma_elem[1, :])
            tau_imag = np.imag(sigma_elem[1, :])
            
            num_points = len(tau_real)
            # Create an array of points along the length of the element
            s = np.linspace(0, 1, num_points)
      
            # Calculate the original positions along the element's length
            x = x0 + (x1 - x0) * s
            z = z0 + (z1 - z0) * s
          # Compute perpendicular direction for offsetting the moment line
            dx = x1 - x0
            dz = z1 - z0
            length = np.sqrt(dx**2 + dz**2)
            perp_x = dz / length  # Unit normal in X
            perp_z = dx / length   # Unit normal in Z          
            
            shearstress_x = x - scale * tau_real * perp_x
            shearstress_z = z + scale * tau_real * perp_z
            
            shearstress_x_imag = x - scale * tau_imag * perp_x
            shearstress_z_imag = z + scale * tau_imag * perp_z
            
            # Plot the shear stress
            plt.plot(shearstress_x, shearstress_z, 'b-', linewidth=1, label='Real Part of Shear Stress' if element.id == elements[0].id else "")
            plt.plot(shearstress_x_imag, shearstress_z_imag, 'g--', linewidth=1, label='Imaginary Part of Shear Stress' if element.id == elements[0].id else "")
            
            # Connect moment diagram back to the structure
            plt.plot([x[0], shearstress_x[0]], [z[0], shearstress_z[0]], 'b--', linewidth=1)
            plt.plot([x[-1], shearstress_x[-1]], [z[-1], shearstress_z[-1]], 'b--', linewidth=1)
            plt.plot([x[0], shearstress_x_imag[0]], [z[0], shearstress_z_imag[0]], 'g--', linewidth=1)
            plt.plot([x[-1], shearstress_z_imag[-1]], [z[-1], shearstress_z_imag[-1]], 'g--', linewidth=1)  
        # Add legend and show the plot
        plt.legend()
        plt.show()

    def Plotaxialstresses(self, elements, stresses, scale=1e-9):
     
        plt.figure(figsize=(10, 8))
        plt.title('Axial Stress Plot')
        plt.grid(True)
    
        # Plot the original structure (nodes and elements)
        for node in [node for element in elements for node in element.nodes]:
            plt.scatter(node.x, node.z, color='red', marker='o', label='Node' if 'Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k--', linewidth=2, label='Element' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")
    
        as_nodes_real = {}
        as_nodes_imag = {}
    
        # Plot the displaced elements (real and imaginary parts)
        for element in elements:
            # Extract the original coordinates of the nodes
            x0 = element.nodes[0].x
            z0 = element.nodes[0].z
            x1 = element.nodes[1].x
            z1 = element.nodes[1].z
            
            # Get the displacement components (assuming u_elem has shape (6, num_points))
            sigma_elem = stresses[element.id]
            sigma_xx_real = np.real(sigma_elem[0, :])
            sigma_xx_imag = np.imag(sigma_elem[0, :])
            
            num_points = len(sigma_xx_real)
            # Create an array of points along the length of the element
            s = np.linspace(0, 1, num_points)
      
            # Calculate the original positions along the element's length
            x = x0 + (x1 - x0) * s
            z = z0 + (z1 - z0) * s
          # Compute perpendicular direction for offsetting the moment line
            dx = x1 - x0
            dz = z1 - z0
            length = np.sqrt(dx**2 + dz**2)
            perp_x = dz / length  # Unit normal in X
            perp_z = dx / length   # Unit normal in Z          
            
            axialstress_x = x + scale * sigma_xx_real * perp_x
            axialstress_z = z - scale * sigma_xx_real * perp_z
            
            axialstress_x_imag = x + scale * sigma_xx_imag * perp_x
            axialstress_z_imag = z + scale * sigma_xx_imag * perp_z
            
            # Plot the shear stress
            plt.plot(axialstress_x, axialstress_z, 'b-', linewidth=1, label='Real Part of Axial Stress' if element.id == elements[0].id else "")
            plt.plot(axialstress_x_imag, axialstress_z_imag, 'g--', linewidth=1, label='Imaginary Part of Axial Stress' if element.id == elements[0].id else "")
            
            # Connect moment diagram back to the structure
            plt.plot([x[0], axialstress_x[0]], [z[0], axialstress_z[0]], 'b--', linewidth=1)
            plt.plot([x[-1], axialstress_x[-1]], [z[-1], axialstress_z[-1]], 'b--', linewidth=1)
            plt.plot([x[0], axialstress_x_imag[0]], [z[0], axialstress_z_imag[0]], 'g--', linewidth=1)
            plt.plot([x[-1], axialstress_z_imag[-1]], [z[-1], axialstress_z_imag[-1]], 'g--', linewidth=1)  
        # Add legend and show the plot
        plt.legend()
        plt.show()

    def Plotbendingstresses(self, elements, stresses, scale=1e-9):
     
        plt.figure(figsize=(10, 8))
        plt.title('Bending Stress Plot')
        plt.xlabel('X Coordinate')
        plt.ylabel('Z Coordinate')
        plt.grid(True)
    
        # Plot the original structure (nodes and elements)
        for node in [node for element in elements for node in element.nodes]:
            plt.scatter(node.x, node.z, color='red', marker='o', label='Node' if 'Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(node.x+0.02, node.z+0.02, f'{node.id}', fontsize=9, ha='right')
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k--', linewidth=2, label='Element' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")
    
        as_nodes_real = {}
        as_nodes_imag = {}
    
        # Plot the displaced elements (real and imaginary parts)
        for element in elements:
            # Extract the original coordinates of the nodes
            x0 = element.nodes[0].x
            z0 = element.nodes[0].z
            x1 = element.nodes[1].x
            z1 = element.nodes[1].z
            
            # Get the displacement components (assuming u_elem has shape (6, num_points))
            sigma_elem = stresses[element.id]
            sigma_yy_real = np.real(sigma_elem[-1, :])
            sigma_yy_imag = np.imag(sigma_elem[-1, :])
            
            num_points = len(sigma_yy_real)
            # Create an array of points along the length of the element
            s = np.linspace(0, 1, num_points)
      
            # Calculate the original positions along the element's length
            x = x0 + (x1 - x0) * s
            z = z0 + (z1 - z0) * s
          # Compute perpendicular direction for offsetting the moment line
            dx = x1 - x0
            dz = z1 - z0
            length = np.sqrt(dx**2 + dz**2)
            perp_x = dz / length  # Unit normal in X
            perp_z = dx / length   # Unit normal in Z          
            
            bendingstress_x = x - scale * sigma_yy_real * perp_x
            bendingstress_z = z + scale * sigma_yy_real * perp_z
            
            bendingstress_x_imag = x - scale * sigma_yy_imag * perp_x
            bendingstress_z_imag = z + scale * sigma_yy_imag * perp_z
            
            # Plot the shear stress
            plt.plot(bendingstress_x, bendingstress_z, 'b-', linewidth=1, label='Real Part of Bending Stress' if element.id == elements[0].id else "")
            plt.plot(bendingstress_x_imag, bendingstress_z_imag, 'g--', linewidth=1, label='Imaginary Part of Bending Stress' if element.id == elements[0].id else "")
            
            # Connect moment diagram back to the structure
            plt.plot([x[0], bendingstress_x[0]], [z[0], bendingstress_z[0]], 'b--', linewidth=1)
            plt.plot([x[-1], bendingstress_x[-1]], [z[-1], bendingstress_z[-1]], 'b--', linewidth=1)
            plt.plot([x[0], bendingstress_x_imag[0]], [z[0], bendingstress_z_imag[0]], 'g--', linewidth=1)
            plt.plot([x[-1], bendingstress_z_imag[-1]], [z[-1], bendingstress_z_imag[-1]], 'g--', linewidth=1)  
        # Add legend and show the plot
        plt.legend()
        plt.show()
    
    # %% PLOT RESULTS LOCAL 
    def plot_element_moment(self, element, forces):
        """Plots the bending moment for a single element."""
        plt.figure(figsize=(8, 6))
        plt.title(f'Bending Moments for Element {element.id}')
        plt.grid(True)

        # Extract the original coordinates of the nodes
        x0, z0 = 0, 0
        x1, z1 = element.L, 0

        # Get moment components
        f_elem = forces[element.id]
        m_real = np.real(f_elem[-1, :])
        m_imag = np.imag(f_elem[-1, :])

        num_points = len(m_real)
        s = np.linspace(0, 1, num_points)

        x = x0 + (x1 - x0) * s
        z = z0 + (z1 - z0) * s

        dx, dz = x1 - x0, z1 - z0
        length = np.sqrt(dx**2 + dz**2)
        perp_x, perp_z = -dz / length, dx / length

        moment_x = x + m_real * perp_x
        moment_z = z + m_real * perp_z
        moment_x_imag = x + m_imag * perp_x
        moment_z_imag = z + m_imag * perp_z

        plt.plot(moment_x, moment_z, 'b-', linewidth=1, label='Real Part')
        plt.plot(moment_x_imag, moment_z_imag, 'g--', linewidth=1, label='Imaginary Part')

        plt.legend()
        plt.show()

    def plot_element_axial_force(self, element, forces):
        """Plots the axial force for a single element."""
        plt.figure(figsize=(8, 6))
        plt.title(f'Axial Force for Element {element.id}')
        plt.grid(True)

        # Extract coordinates
        x0, z0 = 0, 0
        x1, z1 = element.L, 0

        # Get axial force components
        f_elem = forces[element.id]
        N_real = np.real(f_elem[0, :])
        N_imag = np.imag(f_elem[0, :])

        num_points = len(N_real)
        s = np.linspace(0, 1, num_points)

        x = x0 + (x1 - x0) * s
        z = z0 + (z1 - z0) * s

        dx, dz = x1 - x0, z1 - z0
        length = np.sqrt(dx**2 + dz**2)
        perp_x, perp_z = dz / length, dx / length

        axial_x = x + N_real * perp_x
        axial_z = z - N_real * perp_z
        axial_x_imag = x + N_imag * perp_x
        axial_z_imag = z - N_imag * perp_z

        plt.plot(axial_x, axial_z, 'b-', linewidth=1, label='Real Part')
        plt.plot(axial_x_imag, axial_z_imag, 'g--', linewidth=1, label='Imaginary Part')

        plt.legend()
        plt.show()

    def plot_element_shear_force(self, element, forces):
        """Plots the shear force for a single element."""
        plt.figure(figsize=(8, 6))
        plt.title(f'Shear Force for Element {element.id}')
        plt.grid(True)

        # Extract coordinates
        x0, z0 = 0, 0
        x1, z1 = element.L, 0

        # Get shear force components
        f_elem = forces[element.id]
        V_real = np.real(f_elem[1, :])
        V_imag = np.imag(f_elem[1, :])

        num_points = len(V_real)
        s = np.linspace(0, 1, num_points)

        x = x0 + (x1 - x0) * s
        z = z0 + (z1 - z0) * s

        dx, dz = x1 - x0, z1 - z0
        length = np.sqrt(dx**2 + dz**2)
        perp_x, perp_z = dz / length, dx / length

        shear_x = x - V_real * perp_x
        shear_z = z + V_real * perp_z
        shear_x_imag = x - V_imag * perp_x
        shear_z_imag = z + V_imag * perp_z

        plt.plot(shear_x, shear_z, 'b-', linewidth=1, label='Real Part')
        plt.plot(shear_x_imag, shear_z_imag, 'g--', linewidth=1, label='Imaginary Part')

        plt.legend()
        plt.show()

    def plot_element_shear_stress(self, element, stresses):
        """Plots the shear stress for a single element in local coordinates."""
        plt.figure(figsize=(8, 6))
        plt.title(f'Shear Stress for Element {element.id}')
        plt.grid(True)

        x0, z0 = 0, 0
        x1, z1 = element.L, 0  # Local coordinate system

        sigma_elem = stresses[element.id]
        tau_real = np.real(sigma_elem[1, :])
        tau_imag = np.imag(sigma_elem[1, :])

        num_points = len(tau_real)
        s = np.linspace(0, 1, num_points)

        x = x0 + (x1 - x0) * s
        z = np.zeros_like(x)  # Element lies along the x-axis in local coordinates

        dx, dz = x1 - x0, z1 - z0
        length = np.sqrt(dx**2 + dz**2)
        perp_x, perp_z = dz / length, dx / length  # Perpendicular unit vector

        shearstress_x = x - tau_real * perp_x
        shearstress_z = z + tau_real * perp_z
        shearstress_x_imag = x - tau_imag * perp_x
        shearstress_z_imag = z + tau_imag * perp_z

        plt.plot(shearstress_x, shearstress_z, 'b-', linewidth=1, label='Real Part')
        plt.plot(shearstress_x_imag, shearstress_z_imag, 'g--', linewidth=1, label='Imaginary Part')

        plt.legend()
        plt.show()

    def plot_element_axial_stress(self, element, stresses):
        """Plots the axial stress for a single element in local coordinates."""
        plt.figure(figsize=(8, 6))
        plt.title(f'Axial Stress for Element {element.id}')
        plt.grid(True)

        x0, z0 = 0, 0
        x1, z1 = element.L, 0  # Local coordinate system

        sigma_elem = stresses[element.id]
        sigma_xx_real = np.real(sigma_elem[0, :])
        sigma_xx_imag = np.imag(sigma_elem[0, :])

        num_points = len(sigma_xx_real)
        s = np.linspace(0, 1, num_points)

        x = x0 + (x1 - x0) * s
        z = np.zeros_like(x)  # Element lies along the x-axis in local coordinates

        dx, dz = x1 - x0, z1 - z0
        length = np.sqrt(dx**2 + dz**2)
        perp_x, perp_z = dz / length, dx / length  # Perpendicular unit vector

        axialstress_x = x + sigma_xx_real * perp_x
        axialstress_z = z - sigma_xx_real * perp_z
        axialstress_x_imag = x + sigma_xx_imag * perp_x
        axialstress_z_imag = z - sigma_xx_imag * perp_z

        plt.plot(axialstress_x, axialstress_z, 'b-', linewidth=1, label='Real Part')
        plt.plot(axialstress_x_imag, axialstress_z_imag, 'g--', linewidth=1, label='Imaginary Part')

        plt.legend()
        plt.show()

    def plot_element_bending_stress(self, element, stresses):
        """Plots the bending stress for a single element in local coordinates."""
        plt.figure(figsize=(8, 6))
        plt.title(f'Bending Stress for Element {element.id}')
        plt.grid(True)

        x0, z0 = 0, 0
        x1, z1 = element.L, 0  # Local coordinate system

        sigma_elem = stresses[element.id]
        sigma_yy_real = np.real(sigma_elem[-1, :])
        sigma_yy_imag = np.imag(sigma_elem[-1, :])

        num_points = len(sigma_yy_real)
        s = np.linspace(0, 1, num_points)

        x = x0 + (x1 - x0) * s
        z = np.zeros_like(x)  # Element lies along the x-axis in local coordinates

        dx, dz = x1 - x0, z1 - z0
        length = np.sqrt(dx**2 + dz**2)
        perp_x, perp_z = dz / length, dx / length  # Perpendicular unit vector

        bendingstress_x = x - sigma_yy_real * perp_x
        bendingstress_z = z + sigma_yy_real * perp_z
        bendingstress_x_imag = x - sigma_yy_imag * perp_x
        bendingstress_z_imag = z + sigma_yy_imag * perp_z

        plt.plot(bendingstress_x, bendingstress_z, 'b-', linewidth=1, label='Real Part')
        plt.plot(bendingstress_x_imag, bendingstress_z_imag, 'g--', linewidth=1, label='Imaginary Part')

        plt.legend()
        plt.show()
    
    # %% 
    def ShowStructure(self, title='Structure Configuration'):
        plt.title(title)
        plt.xlabel('X Coordinate')
        plt.ylabel('Z Coordinate')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.show()