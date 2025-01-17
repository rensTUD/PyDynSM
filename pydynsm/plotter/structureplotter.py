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
    
    def PlotNodes(self, nodes):
        plt.figure(figsize=(10, 6))
        for node in nodes:
            # marker = 'o' if len(node.constrained_dofs) == 0 else '^'
            # color = 'blue' if len(node.constrained_dofs) == 0 else 'red'
            # plt.scatter(node.x, node.z, color=color, marker=marker, label=f'{"Free" if len(node.constrained_dofs) == 0 else "Fixed"} Node' if f'{"Free" if len(node.constrained_dofs) == 0 else "Fixed"} Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.scatter(node.x, node.z, color='red', marker='o', label=f'Node: {node.id}')
            plt.text(node.x, node.z, f'{node.id}', fontsize=9, ha='right')

    def PlotElements(self, elements):
        for element in elements:
            x_values = [element.nodes[0].x, element.nodes[1].x]
            z_values = [element.nodes[0].z, element.nodes[1].z]
            plt.plot(x_values, z_values, 'k-', linewidth=2, label=f'Element: {element.id}' if 'Element' not in plt.gca().get_legend_handles_labels()[1] else "")
            mid_x = sum(x_values) / 2
            mid_z = sum(z_values) / 2
            plt.text(mid_x, mid_z, element.id, fontsize=9, color='black')
            
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
        plt.axis('equal')

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
            plt.text(x_disp_real, z_disp_real, f'{node_id}', fontsize=9, ha='right', color='blue')

        for node_id, (x_disp_imag, z_disp_imag) in displaced_nodes_imag.items():
            plt.scatter(x_disp_imag, z_disp_imag, color='green', marker='x', label='Imaginary Part of Displaced Node' if 'Imaginary Part of Displaced Node' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(x_disp_imag, z_disp_imag, f'{node_id}', fontsize=9, ha='right', color='green')

        # Add legend and show the plot
        plt.legend()
        plt.show()

    def ShowStructure(self, title='Structure Configuration'):
        plt.title(title)
        plt.xlabel('X Coordinate')
        plt.ylabel('Z Coordinate')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.show()