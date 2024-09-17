# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 19:52:53 2024

@author: rensv
"""

import matplotlib.pyplot as plt

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

    def ShowStructure(self, title='Structure Configuration'):
        plt.title(title)
        plt.xlabel('X Coordinate')
        plt.ylabel('Z Coordinate')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.show()