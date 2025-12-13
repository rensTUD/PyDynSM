# %% import dependencies

import numpy as np
from collections import defaultdict
from ..elements import ElementFactory
from ..sections import SectionFactory
from .dofs import DOFContainer, DOF




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
    section : Section, optional
        Cross-sectional geometry (set via SetSection())
    element_types : dict
        Dictionary of element types added to this element (set via SetElementType())
    
    Examples
    --------
    New two-step workflow:
    
    >>> # Step 1: Set cross-sectional geometry
    >>> element.SetSection('rectangle', {'width': 0.2, 'height': 0.3})
    >>> 
    >>> # Step 2: Set element type with material properties
    >>> element.SetElementType('EulerBernoulli Beam', E=210e9, rho=7850, ksi=0.01)
    >>> 
    >>> # Or with additional element-specific properties:
    >>> element.SetElementType('EulerBernoulli Beam Foundation', 
    ...                        E=210e9, rho=7850, ksi=0.01, kd=1e8, cd=0)
    
    Available section types:
    - 'rectangle': requires {'width', 'height'}
    - 'circle': requires {'diameter'}
    - 'hollow_circle': requires {'outer_diameter', 'inner_diameter'}
    - 'i_section': requires {'flange_width', 'flange_thickness', 'web_height', 'web_thickness'}
    """
    
    # initialise number of elements (ne) as class-variable, to be updated by every new class object
    ne = 0
    
    # dof mapping in full space and its reverse
    dof_mapping = {
        'x': 0,
        'z': 1,
        'y': 2,
        'phi_x': 3,
        'phi_z': 4,
        'phi_y': 5
    }
    
    # Ndof is essentially the length of dof_mapping, which is 6
    maxNdof = len(dof_mapping)
    
    reverse_dof_mapping = {v: k for k, v in dof_mapping.items()}
    
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
                
        # Initialize DOFContainers for global DOFs, 
        self.dof_containers = {}
        
        # initalise constraint type dict - options: monolithic - independent 
        self.constraint_types = {node.id: {} for node in nodes}
        
        # copying DOFs from nodes
        for node in nodes:
            # Load nodal DOFContainer
            node_dof_container = node.dof_container
            # intialise element DOFContainer
            element_dof_container = DOFContainer()
            
            # Copy DOFs from node's DOFContainer to element's DOFContainer
            for dof_name, node_dof in node_dof_container.dofs.items():
                # Create a new DOF instance with the same name and value
                element_dof = DOF(name=dof_name, value=node_dof.value, index=node_dof.index)
                element_dof_container.dofs[dof_name] = element_dof
                
                # set dof to monolithic initially
                self.constraint_types[node.id][dof_name] = 'monolithic'
                
            self.dof_containers[node.id] = element_dof_container
        
        # set local dofs containers as empty first
        self.local_dof_container = {node.id: DOFContainer() for node in nodes}        
        
        # calculate geometrical properties
        self.geometrical_properties()
        
        # Keep track of what kind of elements are included in the form of a dictionary to be able to use them later
        self.element_types = {}
        
        # Store section instance (set via SetSection)
        self.section = None
        
        # initialise empty list with nodal forces due to element loads
        self.element_loads = defaultdict(dict)

# %% Nodes and Dofs

    def get_element_node_dof_indices_global(self):
        return [dof.index for node in self.nodes for dof in self.dof_containers[node.id].dofs.values()]        
    
    def get_node_dof_indices_global(self):
        return [dof_index for node in self.nodes for dof_index in node.get_dof_indices()]   

    def check_dofs(self, dofs):
        """
        Checks which DOFs are present and returns their numerical values.

        Parameters
        ----------
        dofs : list of str
            List of DOFs to check (e.g., ['x', 'z']).

        Returns
        -------
        numerical_dofs : list of int
            List of numerical values corresponding to the DOFs.
        """
        numerical_dofs = []
        for dof_name in dofs:
            if dof_name in self.dof_mapping:
                numerical_dofs.append(self.dof_mapping[dof_name])
            else:
                raise ValueError(f"DOF '{dof_name}' is not recognized.")
        return numerical_dofs

    def get_element_type_dofs(self, element_type):
        '''
        Translates the dofs of an element to the nodal dofs;
        
        Example:
            Rod element has dof = 0 (axial only)
            This will give: nodal_dofs = [0,3]
        '''
        
        element_dofs = self.check_dofs(element_type.dofs)
        
        maxNdof = Element.maxNdof
        
        element_type_dofs = []
        
        # Process each DOF to map it to the nodes' corresponding DOFs
        for i in range(2):
            for dof in element_dofs:
                # Calculate the DOF index for the node
                node_dof_index = dof + i * maxNdof
                element_type_dofs.append(node_dof_index)
                
        return element_type_dofs
    
    
    def get_full_element_dof_indices_global(self) -> list:
        """
        Returns the list of global DOF indices for the element based on the DOFs stored in element.dof_containers
        Based on the numbering provided by class variable 'dof_mapping'
    
        Returns
        -------
        nodal_dofs : list of int
            List of global DOF indices for the element's DOFs.
        """
        element_dof_indices = []
        Ndof = Element.maxNdof  
        
        for i, node in enumerate(self.nodes):
            element_dof_container = self.dof_containers[node.id]
            for dof_name, dof in element_dof_container.dofs.items():
                # Calculate local index if global index is not assigned
                local_dof_index = self.dof_mapping[dof_name] + i * Ndof
                element_dof_indices.append(local_dof_index)
        
        return element_dof_indices
    
    def get_full_element_dof_indices_local(self) -> list:
        """
        Returns the list of global DOF indices for the element based on the DOFs stored in element.local_dof_container
        Based on the numbering provided by class variable 'dof_mapping'
    
        Returns
        -------
        nodal_dofs : list of int
            List of local DOF indices for the element's DOFs.
        """
        element_dof_indices = []
        Ndof = Element.maxNdof  
        
        for i, node in enumerate(self.nodes):
            element_dof_container = self.local_dof_container[node.id]
            for dof_name, dof in element_dof_container.dofs.items():
                # Calculate local index if global index is not assigned
                local_dof_index = self.dof_mapping[dof_name] + i * Ndof
                element_dof_indices.append(local_dof_index)
        
        return element_dof_indices
    
    def get_specific_dof_indices_global(self, dofs) -> list:
        """
        Returns the indices of the specified DOFs within the element’s global displacement vector (u_nodes_global).
    
        This method helps you extract specific DOFs (e.g., 'x', 'z') from the `u_nodes_global` vector in the order
        they appear for this element.
    
        Example:
        --------
        Suppose the element has the following DOFs in order:
            ['x', 'z', 'phi_y', 'x', 'z', 'phi_y']  # for node1 and node2
    
        Then:
            get_specific_dof_indices_global(['z']) → [1, 4]
            get_specific_dof_indices_global(['x']) → [0, 3]
            get_specific_dof_indices_global(['x', 'z']) → [0, 1, 3, 4]
    
        Use these indices to extract the corresponding displacements from `u_nodes_global`, like:
            indices = get_specific_dof_indices_global(['z'])
            uz = u_nodes_global[indices]
        """
        specific_dof_indices = []
        index_counter = 0
    
        for node in self.nodes:
            element_dof_container = self.dof_containers[node.id]  
            for dof_name in element_dof_container.dofs.keys():
                # Check if the current DOF is in the list of required DOFs
                if dof_name in dofs:
                    specific_dof_indices.append(index_counter)
                index_counter += 1
    
        return specific_dof_indices
    
    def get_specific_dof_indices_local(self, dofs) -> list:
        """
        Returns the indices of the specified DOFs within the element’s global displacement vector (u_nodes_global).
    
        This method helps you extract specific DOFs (e.g., 'x', 'z') from the `u_nodes_global` vector in the order
        they appear for this element.
    
        Example:
        --------
        Suppose the element has the following DOFs in order:
            ['x', 'z', 'phi_y', 'x', 'z', 'phi_y']  # for node1 and node2
    
        Then:
            get_specific_dof_indices_global(['z']) → [1, 4]
            get_specific_dof_indices_global(['x']) → [0, 3]
            get_specific_dof_indices_global(['x', 'z']) → [0, 1, 3, 4]
    
        Use these indices to extract the corresponding displacements from `u_nodes_global`, like:
            indices = get_specific_dof_indices_global(['z'])
            uz = u_nodes_global[indices]
        """
        specific_dof_indices = []
        index_counter = 0
    
        for node in self.nodes:
            element_dof_container = self.local_dof_container[node.id]  # Use local DOF container
            for dof_name in element_dof_container.dofs.keys():
                # Check if the current DOF is in the list of required DOFs
                if dof_name in dofs:
                    specific_dof_indices.append(index_counter)
                index_counter += 1
    
        return specific_dof_indices

# %% stiffness

    def Stiffness(self, omega):
        """
        Function to determine the global stiffness matrix of the element, 
        assembling it from its local element contributions.
        """
        # Get global and local DOF indices
        global_dof_indices = self.get_full_element_dof_indices_global()
        local_dof_indices = self.get_full_element_dof_indices_local()
    
        # Initialize full 12x12 stiffness matrix
        K_full = np.zeros((2 * Element.maxNdof, 2 * Element.maxNdof), dtype=complex)
    
        # Loop over all element types and accumulate stiffness contributions
        for element_type_name, element_type in self.element_types.items():
            # Get element-specific DOFs
            element_type_dofs = self.get_element_type_dofs(element_type)
            
            # Add element stiffness to full 12x12 matrix
            K_full[np.ix_(element_type_dofs, element_type_dofs)] += element_type.LocalStiffness(omega)
    
        # Apply the full rotation matrix to transform into the global coordinate system
        K_global_full = (self.R.T @ K_full) @ self.R  
    
        # Extract the relevant global DOFs **after** applying rotation
        K_global = K_global_full[np.ix_(global_dof_indices, global_dof_indices)]
    
        return K_global


# %% loads

    def AddDistributedLoad(self, **loads) -> None:
        '''
        Adds distributed loads to the element.
    
        Parameters
        ----------
        **loads : dict
            Keyword arguments where the key is the DOF (degree of freedom) and the value is the load magnitude.
        '''
        for dof, load in loads.items():
            if dof in self.element_loads:
                print(f'Distributed load for DOF {dof} is already present and will be overwritten. '
                      f'Old value: {self.element_loads[dof]}, New value: {load}')
            self.element_loads[dof] = load
            
    def RemoveDistributedLoad(self, dof: str) -> None:
        '''
        Removes the distributed load for a specific degree of freedom (DOF) from the element.
    
        Parameters
        ----------
        dof : str
            The degree of freedom (e.g., 'x', 'y', 'z', etc.) for which the distributed load should be removed.
        '''
        if dof in self.element_loads:
            removed_load = self.element_loads.pop(dof)
            print(f'Distributed load for DOF {dof} (value: {removed_load}) has been removed.')
        else:
            print(f'No distributed load found for DOF {dof}. Nothing was removed.')
            

    def EvaluateDistributedLoad(self, element_loads, omega):
        '''
        Evaluates the distributed load
        '''
        # get the local dof indices - i.e. the indices of dofs present in the element per node
        local_dof_indices = self.get_full_element_dof_indices_global()
        
        # number of dofs in the current Element
        Ndof = len(local_dof_indices)
        
        # initialise local nodal load vector
        q_loc = np.zeros(Ndof, complex)
        
        
        for element_type_name, element_type in self.element_types.items():
            # get dofs of specific element type
            dofs = element_type.dofs
            
            # get the load into a vector belonging to only the dofs of the current element_type. Evaluate with omega if load is a lambda function
            q_evaluated = [
                element_loads.get(dof,0)(omega) if callable(element_loads.get(dof,0)) else element_loads.get(dof,0)
                for dof in dofs
                ]
            
            # pass through the loads that are specific for that element
            q_loc += self.FullDistributedLoad(element_type, q_evaluated, omega)[np.ix_(local_dof_indices)]
        
        # get global force vector
        q_glob = self.R[np.ix_(local_dof_indices,local_dof_indices)].T @ q_loc
        
        return q_glob   

    def FullDistributedLoad(self, element_type, q, omega):
        '''
        Get the nodal loads from the distributed load per element type and put in the full load vector
        '''
        # get the element_type dofs
        element_type_dofs = self.get_element_type_dofs(element_type)
        
        # initialise maxNdof empty complex vector 
        q_full = np.zeros(2*Element.maxNdof,complex)
        
        # set directly
        q_full[np.ix_(element_type_dofs)] = element_type.LocalDistributedLoad(q, omega)

        
        return q_full             

#%% properties of the Element    
            
    def geometrical_properties(self,alpha = 0):
        '''
        Determines the following geometrical properties of the element:
            - length
            - full rotation matrix (12x12)
    
        the full displacement vector would be as follows:            
            u = [u_x_l, u_z_l, u_y_l, phi_x_l, phi_z_l, phi_y_l, u_x_r, u_z_r, u_y_r, phi_x_r, phi_z_r, phi_y_r]            
            
        Parameters
        ----------
        alpha : rotation in degrees of the cross-section. I ASSUME A ROTATION OF 0 MEANS THE CROSS SECTION IS PERFECTLY VERTICAL? NEED TO CHECK
        '''
        
        ## Calculate length of beam (L) and orientation 
        xl, zl, yl = self.nodes[0].get_coords() # x, z and y of left node

        xr, zr, yr = self.nodes[1].get_coords() # x, z and y of left node

        # get the length of the element
        l = np.sqrt((xr-xl)**2 + (zr-zl)**2 + (yr-yl)**2)           
        self.L = l
        
        # get the angles 
        Cx = (xr-xl)/l
        Cz = (zr-zl)/l
        Cy = (yr-yl)/l
        Cxy = np.sqrt(Cx**2 + Cy**2)
        
        ## Determining rotation if alpha is not given
        
        if alpha is None:
            # get the rotation about the local x-axis according to nonstructural node
            # no nonstructnode is given, thus create one on the same x,y plane and above and outside of it - NEED DEFINITION OF NONSTRUCTNODE
            nonstructnode = []
                
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
        
# %% Setting / removing a section        
    def SetSection(self, section_type, props):
        '''
        Set the cross-sectional geometry for the element.
        
        This method creates a section instance from the provided section type and dimensions.
        The section computes geometric properties (A, I_y, I_z, W_y, W_z) from dimensions.
        Material properties are NOT stored in the section - they are provided separately
        when calling SetElementType().
        
        Parameters
        ----------
        section_type : str
            Type of section (e.g., 'rectangle', 'circle', 'hollow_circle', 'i_section')
        props : dict
            Dictionary containing section dimensions (e.g., {'width': 0.2, 'height': 0.3}
            for rectangle, or {'diameter': 0.1} for circle)
        
        Examples
        --------
        >>> element.SetSection('rectangle', {'width': 0.2, 'height': 0.3})
        >>> element.SetSection('circle', {'diameter': 0.1})
        >>> element.SetSection('hollow_circle', {'outer_diameter': 0.2, 'inner_diameter': 0.15})
        >>> element.SetSection('i_section', {'flange_width': 0.2, 'flange_thickness': 0.02,
        ...                                   'web_height': 0.3, 'web_thickness': 0.01})
        '''
        try:
            # Create section instance via SectionFactory
            section = SectionFactory.CreateSection(section_type, **props)
            
            # Store section instance
            self.section = section
            
            print(f'Successfully set section of type: {section_type} to Element {self.id}')
        except Exception as e:
            print(f'Exception occurred while setting section - {e}')
            raise
    
    def SetElementType(self, element_type, **material_and_element_props):
        '''
        Set the element type with material and element-specific properties.
        
        This method creates an element instance using the stored section and provided
        material properties. The section must be set first using SetSection().
        
        Parameters
        ----------
        element_type : str
            Type of element (e.g., 'Rod', 'EulerBernoulli Beam', etc.)
        **material_and_element_props : dict
            Material properties (E, rho, ksi, G, nu) and element-specific properties
            (kd, cd, T, k, etc.) passed as keyword arguments.
        
        Examples
        --------
        >>> element.SetSection('rectangle', {'width': 0.2, 'height': 0.3})
        >>> element.SetElementType('EulerBernoulli Beam', E=210e9, rho=7850, ksi=0.01)
        >>> element.SetElementType('EulerBernoulli Beam Foundation', E=210e9, rho=7850,
        ...                        ksi=0.01, kd=1e8, cd=0)
        '''
        # Check that section is set
        if self.section is None:
            raise ValueError("Section must be set first using SetSection() before setting element type.")
        
        try:
            # Create element instance by calling element's __init__ with:
            # - section=self.section (section object passed directly)
            # - L=self.L (element length)
            # - **material_and_element_props (material and element-specific properties)
            element = ElementFactory.CreateElement(
                element_type,
                section=self.section,
                L=self.L,
                **material_and_element_props
            )
            
            # Add local DOFs to each node's local DOFContainer
            for node in self.nodes:
                local_container = self.local_dof_container[node.id]
                for dof_name in element.dofs:
                    if not local_container.has_dof(dof_name):
                        local_container.set_dof(dof_name)
            
            # check influence of the local dofs of the element on the global dofs of the current node
            self.apply_global_constraints_to_local_dofs()
            
            # store the element
            self.element_types[element_type] = element
            
            print(f'Successfully added element of type: {element_type} to Element {self.id}')
        except Exception as e:
            print(f'Exception occurred - {e}')
            raise
    
    def RemoveSection(self):
        '''
        Remove the stored section and optionally remove associated element types.
        
        This removes the section instance and all element types that depend on it.
        '''
        try:
            # Remove section
            self.section = None
            
            # Remove all element types (they depend on the section)
            self.element_types.clear()
            
            print(f'Successfully removed section from Element {self.id}')
        except Exception as e:
            print(f'An error has occurred whilst trying to remove a section - {e}')
    
    def RemoveElementType(self, element_type):
        '''
        Remove a specific element type from the element.
        
        Parameters
        ----------
        element_type : str
            Type of element to remove
        '''
        try:
            del self.element_types[element_type]
            print(f'Successfully removed element type {element_type} from Element {self.id}')
        except KeyError:
            print(f'Element type {element_type} not found in Element {self.id}')
        except Exception as e:
            print(f'An error has occurred whilst trying to remove element type - {e}')

                
# %% Prescribing dofs   
    def fix_dof(self, node, *dofs):
        """
        Fixes specified DOFs at a given node in the element.
    
        Parameters
        ----------
        node : Node
            The node where DOFs are to be fixed.
        *dofs : str
            DOFs to be fixed at the node.
        """
        
        self.prescribe_dof(node, **{dof_name: 0 for dof_name in dofs})
                
    def free_dof(self, node, *dofs):
        """
        Frees specified DOFs at a given node in the element.
    
        Parameters
        ----------
        node : Node
            The node where DOFs are to be freed.
        *dofs : str
            DOFs to be freed at the node.
        """
        self.prescribe_dof(node, **{dof_name: None for dof_name in dofs})
                
        
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
            If any of the specified DOFs cannot be added due to lack of influence by local DOFs.
        """
        try:
            dof_container = self.dof_containers[node.id]  # Access the node's DOFContainer
            for dof_name, value in dofs.items():
                # Check if the DOF is already in the node's DOFContainer
                if dof_container.has_dof(dof_name):
                    # Update the DOF with the new value
                    self._update_dof(node, dof_name, value)
                else:
                    # Get list of global dofs affected by local dofs
                    global_dofs = self.map_local_to_global_dofs(self.local_dof_container[node.id].dofs.keys())
                    # Flatten the list of global DOFs influenced by local DOFs
                    influenced_global_dofs = [gdof for gdofs in global_dofs.values() for gdof in gdofs]
                    
                    # Check if the global DOF is influenced by local DOFs
                    if dof_name in influenced_global_dofs:
                        self._handle_missing_dof(node, dof_name, value)
                    else:
                        print(f"DOF '{dof_name}' cannot be added, as it is not affected by local DOFs in the current element configuration.")
        except ValueError as e:
            print(e)
    
    
    def _update_dof(self, node, dof_name, value):
        """
        Updates the DOF for a given node, also updating local DOFs if necessary.
    
        Parameters
        ----------
        node : Node
            The node to update.
        dof_name : str
            The DOF to be updated.
        value : Any
            The value to prescribe to the DOF.
        """
        # Update the global DOF in the element's DOFContainer
        dof_container = self.dof_containers[node.id]
        dof_container.set_dof(dof_name, value=value)
    
        # Update the local DOFs based on the new global DOF
        local_dofs = self.map_global_to_local_dofs([dof_name])
    
        for local_dof_name in local_dofs.get(dof_name, []):
            if self.local_dof_container[node.id].has_dof(local_dof_name):
                self.local_dof_container[node.id].set_dof_value(local_dof_name, value=value)

    def couple_dof(self, node, *dofs):
        """
        Ensures that the specified DOFs are monolithically connected to the node.
        This means the DOF will be fully constrained across all elements at this node.
    
        Parameters
        ----------
        node : Node
            The node where the DOFs should be coupled.
        *dofs : str
            DOFs that should be monolithically connected.
        """
        for dof_name in dofs:
            if dof_name in self.dof_containers[node.id].dofs:
                self.constraint_types[node.id][dof_name] = "monolithic"
    
    
    def decouple_dof(self, node, *dofs):
        """
        Allows the specified DOFs to move freely in the element, preventing monolithic connection.
    
        Parameters
        ----------
        node : Node
            The node where the DOFs should be decoupled.
        *dofs : str
            DOFs that should be independent in the element.
        """
        for dof_name in dofs:
            if dof_name in self.dof_containers[node.id].dofs:
                self.constraint_types[node.id][dof_name] = "independent"

    
    def _handle_missing_dof(self, node, dof_name, value):
        """
        Handles the case where a global DOF is not present in the element but can be influenced by local DOFs.
    
        Parameters
        ----------
        node : Node
            The node to update.
        dof_name : str
            The missing DOF.
        value : Any
            The value to prescribe to the DOF.
        """
        # Check for user input
        response = input(f"DOF '{dof_name}' is not in the element. Would you like to add it? (yes/no): ").strip().lower()
        if response == 'no':
            return
    
        # If answer is yes, add missing DOF if influenced by local DOFs
        dof_container = self.dof_containers[node.id]
        dof_container.set_dof(dof_name, value=value)
        print(f"DOF '{dof_name}' added to the element and prescribed to {value}.")
    
        # Update local DOFs based on the new global DOF
        local_dofs = self.map_global_to_local_dofs([dof_name])
    
        for local_dof_name in local_dofs.get(dof_name, []):
            if self.local_dof_container[node.id].has_dof(local_dof_name):
                self.local_dof_container[node.id].set_dof_value(local_dof_name, value=value)
        print(f"Local DOFs updated to reflect new global DOF '{dof_name}'.")
         

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
            dof_container = self.dof_containers[node.id]
            for dof_name, value in changes.items():
                if dof_container.has_dof(dof_name):
                    dof_container.set_dof(dof_name, value=value)
                    
# %% Element Displacements

    def Displacements(self, u_nodes_global, omega, num_points=20, local_axes = False):
        """
        Gets the displacements over the local axis of the element evaluated at num_points.
    
        Result is structured as (example in 2D config):
            u_elem = [u_x(s)
                      u_z(s)
                      phi_y(s)]
    
        Where:
            - u_x(s) = displacement in x
            - u_z(s) = displacement in z
            - phi_y(s) = rotation about y axis
    
        Parameters
        ----------
        u_nodes_global : np.ndarray
            Global nodal displacements of the element.
        omega : float
            Frequency (used by element types).
        num_points : int, optional
            Number of evaluation points along the element (default is 20).
        local_axes : bool, optional
            Determines whether returned response is in global axes or local axes system
        """
    
        # Step 1 — Extract local/global DOF indices
        local_dof_indices = self.get_full_element_dof_indices_local()
        global_dof_indices = self.get_full_element_dof_indices_global()
        Ndof = Element.maxNdof
    
        # Step 2 — Project global nodal displacements into local frame
        u_local = self.R[np.ix_(local_dof_indices, global_dof_indices)] @ u_nodes_global
    
        # Step 3 — Initialize displacement array in local frame
        u_elem = np.array([np.zeros(num_points, dtype=complex) for _ in range(Ndof)])
    
        # Step 4 — Collect local DOFs that contribute via element_type
        present_local_dofs = set()
        for element_type in self.element_types.values():
            present_local_dofs.update(element_type.dofs)
    
        # Step 5 — Accumulate displacements from element_types
        for element_type in self.element_types.values():
            dofs = element_type.dofs
            specific_dof_indices = self.get_specific_dof_indices_local(dofs)
            u_elem_contribution = element_type.LocalElementDisplacements(
                u_local[np.ix_(specific_dof_indices)], omega, num_points
            )
            for i, dof in enumerate(dofs):
                global_dof_index = self.dof_mapping[dof]
                if u_elem_contribution[i] is not None:
                    u_elem[global_dof_index] += u_elem_contribution[i]
    
        # Step 6 — Determine expected DOFs from node configuration
        node_left, node_right = self.nodes
        node_dof_set = set(node_left.dof_container.dofs.keys()).union(set(node_right.dof_container.dofs.keys()))
        expected_linear_dofs = [dof for dof in ['x', 'z', 'y'] if dof in node_dof_set]
    
        # Step 7 — Interpolate missing transverse DOFs (e.g., 'z')
        if 'z' in expected_linear_dofs and 'z' not in present_local_dofs:
            z_indices = self.get_specific_dof_indices_global(['z'])
            if len(z_indices) == 2:
                # get global displacements and transform to local
                uz_local = self.R[np.ix_([1, 7],global_dof_indices)] @ u_nodes_global
                uz_local_l = uz_local[0]
                uz_local_r = uz_local[1]
                
                # linearly interpolate to ensure rigid body motion and rotation (no curvature!)
                ξ = np.linspace(0, 1, num_points)
                uz_interp = (1 - ξ) * uz_local_l + ξ * uz_local_r
    
                z_index = self.dof_mapping['z']
                u_elem[z_index] += uz_interp
    
        # Step 8 — Axial consistency check: if 'x' DOF is missing but expected
        if 'x' in expected_linear_dofs and 'x' not in present_local_dofs:
            x_indices = self.get_specific_dof_indices_global(['x'])
            if len(x_indices) == 2:
                # get global displacements and transform to local
                ux_local = self.R[np.ix_([0, 6],global_dof_indices)] @ u_nodes_global
                ux_local_l = ux_local[0]
                ux_local_r = ux_local[1]
    
                if not np.isclose(ux_local_l, ux_local_r, atol=1e-6):
                    print(f"⚠️  Element {self.id}: x-displacements of nodes differ but element has no axial DOF — check model consistency.")
        
        # Transform displacements back to global frame if asked for
        if not local_axes:
            u_global = self.R[:6, :6].T @ u_elem
            return u_global
        # Or return displacements in local coordinate frame
        if local_axes:
            return u_elem
    
# %% Element Forces
    def Forces(self, u_nodes_global, omega, num_points=20):
        """
        Gets the Forces over the local axis of the element evaluated at num_points.
    
        Result is structured as (example in 2D config):
            u_elem = [Nx(s)
                      Vz(s)
                      Myy(s)]
    
        Where:
            - Nx(s) = axial force along x axis 
            - Vz(s) = shear force along z axis 
            - Myy(s) = bending moment about y-y axis
    
        Input:
            u_nodes_global: [u_x_left
                             u_z_left
                             phi_y_left
                             u_x_right
                             u_z_right
                             phi_y_right]
        """
        # Get all element dof indices present based on dof mapping
        local_dof_indices = self.get_full_element_dof_indices_local()
        global_dof_indices = self.get_full_element_dof_indices_global()
        
        # Convert the global nodal displacements to the local element coordinate system
        u_local = self.R[np.ix_(local_dof_indices, global_dof_indices)] @ u_nodes_global
        
        # Initialize empty list to hold displacements for each DOF at specified points
        F_elem = np.array([np.zeros(num_points, dtype=complex) for _ in range(Element.maxNdof)])
        
        # Loop over all element types
        for element_type_name, element_type in self.element_types.items():
            # Get the DOFs of the element type in the context of the nodes
            dofs = element_type.dofs
            specific_dof_indices = self.get_specific_dof_indices_local(dofs)
            
            # Calculate local displacements for the specific element type
            F_elem_contribution = element_type.LocalElementForces(
                u_local[np.ix_(specific_dof_indices)], omega, num_points
            )
            
            # Add the contribution of each element type to the full local elment displacements
            for i, dof in enumerate(dofs):
                global_dof_index = self.dof_mapping[dof]
                if F_elem_contribution[i] is not None: 
                    F_elem[global_dof_index] += F_elem_contribution[i]
                    
        F_global =  F_elem
   
        return F_global
    
# %% Element Stresses
    def Stresses(self, u_nodes_global, omega, num_points=20):
            """
            Gets the Stresses over the local axis of the element evaluated at num_points.
        
            Result is structured as (example in 2D config):
                sigma_elem = [sigmaxx(s)
                              tau(s)]
        
            Where:
                - sigmaxx(s) = axial stress along x axis 
                - tau(s) = averaged shear stress
        
            Input:
                u_nodes_global: [u_x_left
                                 u_z_left
                                 phi_y_left
                                 u_x_right
                                 u_z_right
                                 phi_y_right]
            """
            # Get all element dof indices present based on dof mapping
            local_dof_indices = self.get_full_element_dof_indices_local()
            global_dof_indices = self.get_full_element_dof_indices_global()
            
            # Convert the global nodal displacements to the local element coordinate system
            u_local = self.R[np.ix_(local_dof_indices, global_dof_indices)] @ u_nodes_global
            
            # Initialize empty list to hold displacements for each DOF at specified points
            sigma_elem = np.array([np.zeros(num_points, dtype=complex) for _ in range(Element.maxNdof)])
            
            # Loop over all element types
            for element_type_name, element_type in self.element_types.items():
                # Get the DOFs of the element type in the context of the nodes
                dofs = element_type.dofs
                specific_dof_indices = self.get_specific_dof_indices_local(dofs)
                
                # Calculate local displacements for the specific element type
                sigma_elem_contribution = element_type.LocalElementStresses(
                    u_local[np.ix_(specific_dof_indices)], omega, num_points
                )
                
                # Add the contribution of each element type to the full local elment displacements
                for i, dof in enumerate(dofs):
                    global_dof_index = self.dof_mapping[dof]
                    if sigma_elem_contribution[i] is not None: 
                        sigma_elem[global_dof_index] += sigma_elem_contribution[i]
                        
            sigma_global = sigma_elem
            
            return sigma_global
        
# %% help functions    
    def create_local_dofs_vector(self, node):
        """
        Creates a vector with ones where the local DOFs are present and zeros elsewhere.

        Returns
        -------
        numpy.ndarray
            A vector of size equal to the number of DOFs, with ones indicating the presence of local DOFs.
        """
        dof_vector = np.zeros(len(self.dof_mapping))
        for dof_name in self.local_dof_container[node.id].dofs:
            if dof_name in self.dof_mapping:
                index = self.dof_mapping[dof_name]
                dof_vector[index] = 1
        return dof_vector

    def local_to_global_dof(self, local_dof_vector):
        """
        Transforms a local DOF vector into the global DOF vector using the rotation matrix.
    
        Parameters
        ----------
        local_dof_vector : numpy.ndarray
            The local DOF vector.
    
        Returns
        -------
        global_dof_vector : numpy.ndarray
            The global DOF vector.            
        """
        
        global_dof_vector = self.R[:6,:6].T @ local_dof_vector 
        return global_dof_vector
    
    def global_to_local_dof(self, global_dof_vector):
        """
        Transforms a global DOF vector into the local DOF vector using the inverse of the rotation matrix.
    
        Parameters
        ----------
        global_dof_vector : numpy.ndarray
            The global DOF vector.
    
        Returns
        -------
        local_dof_vector : numpy.ndarray
            The local DOF vector.
        """
        
        # Apply the rotation matrix transpose (which is the inverse for orthogonal matrices)
        local_dof_vector = self.R[:6, :6] @ global_dof_vector
        return local_dof_vector
    
    def map_local_to_global_dofs(self,local_dof_names):
        """
        Maps local DOFs to affected global DOFs.
    
        Parameters
        ----------
        local_dof_names : list or iterable of str
            Names of local DOFs.
    
        Returns
        -------
        dict
            Mapping from local DOF names to lists of affected global DOF names.
        """
        affected_global_dofs = {}
        for local_dof_name in local_dof_names:
            local_dof_vector = np.zeros(len(self.dof_mapping))
            local_dof_vector[self.dof_mapping[local_dof_name]] = 1
    
            # Transform to global DOF space
            global_dof_vector = self.local_to_global_dof(local_dof_vector)
    
            # Identify the affected global DOFs
            affected_indices = np.nonzero(global_dof_vector)[0]
            global_dof_names = [self.reverse_dof_mapping[i] for i in affected_indices]
            affected_global_dofs[local_dof_name] = global_dof_names
        return affected_global_dofs
    
    def map_global_to_local_dofs(self, global_dof_names):
        """
        Maps global DOFs to affected local DOFs.
    
        Parameters
        ----------
        global_dof_names : list
            List of global DOF names.
    
        Returns
        -------
        dict
            Mapping from global DOF names to lists of affected local DOF names.
        """
        affected_local_dofs = {}
        for global_dof_name in global_dof_names:
            global_dof_vector = np.zeros(len(self.dof_mapping))
            global_dof_vector[self.dof_mapping[global_dof_name]] = 1
    
            # Transform to local DOF space
            local_dof_vector = self.global_to_local_dof(global_dof_vector)
    
            # Identify the affected local DOFs
            affected_indices = np.nonzero(local_dof_vector)[0]
            local_dof_names = [self.reverse_dof_mapping[i] for i in affected_indices]
            affected_local_dofs[global_dof_name] = local_dof_names
        return affected_local_dofs
    
    def apply_global_constraints_to_local_dofs(self):
        """
        Maps global DOF constraints to local DOFs and applies these constraints to the local DOFs
        of the element.
        """
        for node in self.nodes:
            # Access global DOFs from the node's DOFContainer
            global_container = node.dof_container
            global_dof_names = [dof_name for dof_name in global_container.dofs]
            local_dof_map = self.map_global_to_local_dofs(global_dof_names)
            
            local_dof_container = self.local_dof_container[node.id]
            
            # Apply constraints from global DOFs to local DOFs
            for global_dof_name, local_dof_names in local_dof_map.items():
                global_dof = global_container.get_dof(global_dof_name)
                global_dof_value = global_dof.value
    
                if global_dof_value is not None:
                    for local_dof_name in local_dof_names:
                        if local_dof_container.has_dof(local_dof_name):
                            local_dof_container.set_dof(local_dof_name, value=global_dof_value)

    def apply_global_dof_change(self, node, global_dof_name, global_dof_value):
        """
        Applies the change in a specific global DOF to both the global and local DOFs
        of the element.
        
        Parameters
        ----------
        node : Node
            The node whose global DOF is being changed.
        global_dof_name : str
            The name of the global DOF that has changed.
        global_dof_value : float
            The new value of the global DOF (fixed, free, or prescribed).
        """
        # Update the global DOF in the node's DOFContainer
        dof_container = self.dof_containers[node.id]
        if dof_container.has_dof(global_dof_name):
            dof_container.set_dof_value(global_dof_name, global_dof_value)
            print(f"Global DOF '{global_dof_name}' for node {node.id} updated to {global_dof_value}.")
        
        # Map the changed global DOF to the local DOFs
        local_dof_map = self.map_global_to_local_dofs([global_dof_name])
        
        local_dof_container = self.local_dof_container[node.id]
        
        # Apply the change to the local DOFs
        for local_dof_name in local_dof_map.get(global_dof_name, []):
            if local_dof_container.has_dof(local_dof_name):
                local_dof_container.set_dof_value(local_dof_name, global_dof_value)
                print(f"Local DOF '{local_dof_name}' updated to {global_dof_value}.")
    
    def check_and_handle_global_dofs(self,dofs):
        """
        Checks whether the affected global DOFs are present in self.dofs. If not, prompts for user action.

        Returns
        -------
        bool
            True if all affected global DOFs are present or user chooses to continue, False otherwise.
        """
        affected_dofs = self.check_global_dof_influence(dofs)
        missing_dofs = []

        for local_dof, global_dofs in affected_dofs.items():
            for global_dof in global_dofs:
                for node_id, node_dofs in self.dofs.items():
                    if global_dof not in node_dofs:
                        missing_dofs.append((node_id, global_dof))

        # if all dofs present, do nothing and return 
        if not missing_dofs:
            return 

        print(f"Missing global DOFs: {missing_dofs}")
        for node_id, missing_dof in missing_dofs:    
           # Check if the missing DOF is present in the other node's DOFs
           other_node_id = next(key for key in self.dofs.keys() if key != node_id)
           if missing_dof in self.dofs[other_node_id]:
               print(f"The missing DOF {missing_dof} at node_id {node_id} is present in the other node_id {other_node_id}.")
               action = input(f"The global DOF {missing_dof} at node_id {node_id} is missing. It is present in node_id {other_node_id}. Do you want to do something?: ").strip().lower()
           else:
               action = input(f"The global DOF {missing_dof} at node_id {node_id} is missing. Do you want to do something?: ").strip().lower()
           
           if action == 'yes':
               return True
    
        return False