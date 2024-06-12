# %% import dependencies

import numpy as np

from ...elements import ElementFactory

# %% class definition
class Element:
    """
    A class to represent an element in a structural system.

    Attributes
    ----------
    ne : int
        Class variable to count the number of elements.
    id : int
        Unique identifier for the element.
    nodes : list of Node
        List of nodes connected to this element.
    dofs : dict
        Dictionary of DOFs for each node in the element.
    """
    
    # initialise number of elements (ne) as class-variable, to be updated by every new class object
    ne = 0
    
    def __init__(self, nodes):
        """
        Initializes the Element with its connected nodes.

        Parameters
        ----------
        nodes : list of Node
            List of nodes connected to this element.
        """
        
        # Element number
        self.id = Element.ne
        Element.ne += 1 # increment element number
        
        # reference the nodes
        self.nodes = nodes
        
        # connect the element to the node
        for node in nodes:
            node.connect_element(self)  # Connect this element to its nodes
                
        # assign same dofs of node on initialisation
        self.dofs = {node.id: node.dofs.copy() for node in nodes}
        
        # calculate geometrical properties
        self.geometrical_properties()
        
    def fix_dof(self, node, *dofs):
        """
        Fixes specified DOFs at a given node in the element.

        Parameters
        ----------
        node : Node
            The node where DOFs are to be fixed.
        *dofs : str
            DOFs to be fixed at the node.

        Raises
        ------
        ValueError
            If any of the specified DOFs are not in the node's current configuration.
        """
        
        self.prescribe_dof(node, **{dof: 0 for dof in dofs})
                
    def free_dof(self, node, *dofs):
        """
        Frees specified DOFs at a given node in the element.

        Parameters
        ----------
        node : Node
            The node where DOFs are to be freed.
        *dofs : str
            DOFs to be freed at the node.

        Raises
        ------
        ValueError
            If any of the specified DOFs are not in the node's current configuration.
        """
        
        self.prescribe_dof(node, **{dof: None for dof in dofs})
                
    def prescribe_dof(self, node, **dofs):
        """
        Prescribes specified DOFs at a given node in the element to given values.

        Parameters
        ----------
        node : Node
            The node where DOFs are to be prescribed.
        **dofs : dict
            DOFs and their values to be prescribed at the node.

        Raises
        ------
        ValueError
            If any of the specified DOFs are not in the node's current configuration.
        """
        try:
            for dof, value in dofs.items():
                if dof not in node.dofs:
                    raise ValueError(f"DOF '{dof}' is not available in the node's current configuration. Available DOFs: {node.dof_configurations[node.config]}")
                self.dofs[node.id][dof] = value
        except ValueError as e:
            print(e)              

    def update_node_dofs(self, node, changes):
        """
        Updates the DOFs of a given node in the element based on changes.

        Parameters
        ----------
        node : Node
            The node whose DOFs are to be updated.
        changes : dict
            Dictionary of DOFs and their new values.
        """
        if node in self.nodes:    
            for dof, value in changes.items():           
                self.dofs[node.id][dof] = value
                
    def geometrical_properties(self):
        '''
        Determines the following geometrical properties of the element:
            - length
            - full rotation matrix (12x12)
    
        the full displacement vector would be as follows:            
            u = [u_x_l, u_z_l, u_y_l, phi_x_l, phi_z_l, phi_y_l, u_x_r, u_z_r, u_y_r, phi_x_r, phi_z_r, phi_y_r]            
            
        Parameters
        ----------
        none
        '''
        
        ## 1 - calculate length of beam (L) and orientation 
        xl, zl, yl = self.nodes[0].get_coords() # x, z and y of left node

        xr, zr, yr = self.nodes[0].get_coords() # x, z and y of left node

        # get the length of the element
        l = np.sqrt((xr-xl)**2 + (zr-zl)**2 + (yr-yl)**2)           
        self.L = l
        
        # get the angles 
        Cx = (xr-xl)/l
        Cz = (zr-zl)/l
        Cy = (yr-yl)/l
        Cxy = np.sqrt(Cx**2 + Cy**2)
        
        ## Determining rotation

        # get the rotation about the local x-axis according to nonstructural node:
        
        # no nonstructnode is given, thus create one on the same x,y plane and above and outside of it
        nonstructnode = [xl, zl+2, yl+2]
            
        # first vector for determining global plane
        P1 = np.array(nonstructnode) - np.array([xl, zl, yl])  
        # second vector for determining global plane
        P2 = np.array([xr, zr, yr]) - np.array([xl, zl, yl])  
        # normal of the plane determined by nonstructural node
        normal1 = np.cross(P1,P2)          
        # vectors for plane in local unrotated plane
        P3 = np.array([1, 0, 0])                     
        P4 = np.array([1, 1, 0])                     
        
        # calculate the rotation
        if Cxy == 0:  # if element is placed vertically
            YY = np.array([[0, Cz, 0],
                           [-Cz, 0, 0],
                           [0, 0, 1]])
        else:
            # first 2 rotation matrices for basic element rotation neglecting rotation
            # about the local x-axis
            Rb = np.array([[Cx/Cxy, 0, Cz/Cxy],
                           [0, 1, 0],
                           [-Cy/Cxy, 0, Cx/Cxy]])
            Rz = np.array([[Cxy, Cz, 0],
                           [-Cz, Cxy, 0],
                           [0, 0, 1]])      
            YY = Rz @ Rb
            
        # rotate the basic vectors
        P3 = YY.T @ P3; P3 = P3.T  
        P4 = YY.T @ P4; P4 = P4.T
        # determine normal of the element plane without rotation
        normal2 = np.cross(P4,P3)  
    
        # now calculate the angle between the 2 planes:
        alpha = np.degrees(np.arctan2(np.linalg.norm(np.cross(normal1,normal2)), np.dot(normal1,normal2)))
        
        if alpha != 0:
            print(f'alpha = {alpha} !')
    
        # and do the rotation about all axes
        if Cxy == 0: # check whether member will be placed vertically
            Y = np.array([[0, Cz, 0], [-Cz*np.cos(np.radians(alpha)), 0, np.sin(np.radians(alpha))], [Cz*np.sin(np.radians(alpha)), 0, np.cos(np.radians(alpha))]])
        else:
            Y = np.array([[Cx, Cz, Cy], [(-Cx*Cz*np.cos(np.radians(alpha))-Cy*np.sin(np.radians(alpha)))/Cxy, Cxy*np.cos(np.radians(alpha)), (-Cz*Cy*np.cos(np.radians(alpha))+Cx*np.sin(np.radians(alpha)))/Cxy], [(Cx*Cz*np.sin(np.radians(alpha))-Cy*np.cos(np.radians(alpha)))/Cxy, -Cxy*np.sin(np.radians(alpha)), (Cz*Cy*np.sin(np.radians(alpha))+Cx*np.cos(np.radians(alpha)))/Cxy]])
    
        # full rotation matrix (12x12)
        self.R = np.block([[Y, np.zeros((3,9))], [np.zeros((3,3)), Y, np.zeros((3,6))], [np.zeros((3,6)), Y, np.zeros((3,3))], [np.zeros((3,9)), Y]])
        
    
                        
        
    
    