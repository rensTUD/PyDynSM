# %% import dependencies

import numpy as np

from ...elements import ElementFactory

# %% class definition
class Element:
    ne = 0
    
    # degrees of freedom [x, z, phi]
    Ndof = 3
    
    dof_mapping = {
        'x': 0,
        'z': 1,
        'phi_y': 2
    }
    
    def Clear():
        Element.ne = 0
        
    
    def __init__ (self, nodes):
                
        # assign nodes
        self.nodes = nodes
        
        # determine length of the element
        self.L = np.sqrt((nodes[1].x - nodes[0].x)**2.0 + (nodes[1].z - nodes[0].z)**2.0)
        
        # determine the rotation matrix and assign
        dx = nodes[1].x - nodes[0].x
        dz = nodes[1].z - nodes[0].z

        self.cos = dx / self.L
        self.sin = dz / self.L

        R = np.zeros ((6,6))

        R[0,0] = R[1,1] = R[3,3] = R[4,4] = self.cos
        R[0,1] = R[3,4] = -self.sin
        R[1,0] = R[4,3] =  self.sin
        R[2,2] = R[5,5] = 1.0
        
        # set rotation matrices
        self.R  = R
        self.Rt = np.transpose(R)
        
        # increase element number total
        Element.ne += 1
        
        self.name = f"Element {Element.ne}"
        
        # Keep track of what kind of elements are included in the form of a dictionary to be able to use them later
        self.element_types = {}
        
        # initialise empty list with nodal forces due to element loads
        self.element_nodal_loads = []
        


    def SetSection(self, element_type, props):
        '''
        This function serves to set the elements, we will have to think how we do this properly but this is an initial set-up
        
        '''
        # TODO - need to handle errors correctly when not the correct parameters are given
        
        # add L to props
        props['L'] = self.L
        
        try:
            
            #  try to load the element from the element factory                
            element = ElementFactory.CreateElement(element_type, **props)
            
            # store the element
            self.element_types[element_type] = element    
        except Exception as e:
            print(f'Exception occurred - {e}')
        else:
            print(f'Successfully added element of type: {element_type}')
    
    def RemoveSection(self,element_type):
        '''
        does what it says
        '''        
        try:
            del self.element_types[element_type]
        except Exception as e:
            print(f'An error has occurred whilst trying to remove a section - {e}')

    def GlobalDofs(self):
        return np.hstack ((self.nodes[0].dofs, self.nodes[1].dofs))

    def Stiffness(self, omega):
        '''
        function to determine the global stiffness matrix of the element, as created from its local elements
        '''
        
        # intialise empty stiffness matrix
        k_loc = np.zeros( (6, 6), dtype=complex) 
        
        # loop over all present elements and add their contribution to the full matrix
        for element_type, element in self.element_types.items():
            # k_loc += element.FullStiffness(omega)
            k_loc += self.FullStiffness(element,omega)
        
        # return the full global stiffness matrix
        k_glob = ( self.Rt @ k_loc ) @ self.R 
                
        return k_glob

    def AddDistributedLoad(self, q):
        '''
        function to add distributed load to the element. 
        
        q = [q_x;
             q_z;
             q_phi;]
        
        input should be lambda functions depending on omega or 0
        '''
        
        # TODO - THIS OVERWRITES THE OTHER TODOS: FOR NOW JUST ASSIGN THE LOADS HERE. HAVE ANOTHER FUNCTION EVALUATE THEM
        
    
        self.element_nodal_loads.append(q)
    
    def EvaluateDistributedLoad(self, q, omega):
        '''
        Evaluates the distributed load
        '''
        # initialise local nodal load
        q_loc = np.zeros(2*Element.Ndof, complex)
        
        # evaluate, if necessary, the nodal loads and put in array
        q_evaluated = np.array([load(omega) if callable(load) else load for load in q])
        
        for element_type, element in self.element_types.items():
            # get dofs of specific element type
            dofs = self.check_dofs(element.dofs)
            
            # pass through the loads that are specific for that element
            # q_loc += element.FullDistributedLoad(q_evaluated[dofs], omega)
            q_loc += self.FullDistributedLoad(element, q_evaluated[dofs], omega)
        
        q_glob = self.Rt @ q_loc
        
        return q_glob
    
    def Displacements(self, u_nodes_global, omega, num_points = 2):
        '''
        Gets the displacements over the local axis of the element evaluated at num_points
        
        result is structured as:
            u_elem = [u(s)
                 w(s)
                 phi(s)]
    
        Where:
            - u(s) = axial displacement
            - w(s) = transverse displacement
            - phi(s) = rotation
        '''
        
        # intialise empty list with size Ndof to hold the local displacements of the elemnent
        u_elem = [None] * Element.Ndof
        
        
        # loop over all element types
        for element_type, element in self.element_types.items():
            # get element displacement vector
            # u_elem_1 = element.FullElementDisplacements(u_nodes_global, omega, num_points)
            u_elem_1 = self.FullElementDisplacements(element, u_nodes_global, omega, num_points)
            
            # loop over the displacements and assign whenever result is not None
            for i, u_elem_i in enumerate(u_elem_1):
                if u_elem_i is not None:
                    u_elem[i] = u_elem_i
        
        # Handle None values in u_elem (e.g. when no axial motion is present) and return that
        return [np.zeros(num_points) if u_elem_i is None else u_elem_i for u_elem_i in u_elem]

# %% Element stuff

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
        for dof in dofs:
            if dof in self.dof_mapping:
                numerical_dofs.append(self.dof_mapping[dof])
            else:
                raise ValueError(f"DOF '{dof}' is not recognized.")
        return numerical_dofs

    def GetNodalDofs(self, element):
        '''
        Translates the dofs of an element to the nodal dofs;
        
        Example:
            Rod element has dof = 0 (axial only)
            This will give: nodal_dofs = [0,3]
        '''
        
        dofs = self.check_dofs(element.dofs)
        
        nodal_dofs = []
        # Process each DOF to map it to the current and next node's corresponding DOFs
        for dof in dofs:
            # Add the DOF for the current node
            nodal_dofs.append(dof)
        # After processing each DOF for the current node, add their counterparts in the next node
        for dof in dofs:
            nodal_dofs.append(dof + Element.Ndof)
                
        return nodal_dofs

    def FullStiffness(self, element, omega):
        '''
        Function that assembles the full stiffness matrix based on the local stiffness matrix. 
        
        For example, it will translate the 2x2 K_local of the rod to the full 6x6 matrix which it is in 2D
        
        '''
        
        # get nodal dofs
        nodal_dofs = self.GetNodalDofs(element)
        
        # initialise 6x6 empty complex matrix
        K_full = np.zeros((6,6), complex)
        
        # calculate local stiffness matrix of the element
        K_local = element.LocalStiffness(omega)
        
        # assign to correct locations in full matrix
        for i, row in enumerate(nodal_dofs):
            for j, col in enumerate(nodal_dofs):
                K_full[row, col] = K_local[i, j]
                
        return K_full 
    
    def FullDistributedLoad(self, element, q, omega):
        
        # get nodal dofs
        nodal_dofs = self.GetNodalDofs(element)
        
        # initialise 6x1 empty complex vector 
        q_full = np.zeros(6,complex)
        
        # calculate local load vector
        q_local = element.LocalDistributedLoad(q, omega)
        
        # assign to correct locations in full vector
        for i, row in enumerate(nodal_dofs):
            q_full[row] = q_local[i]
        
        return q_full
    
    def FullDisplacement(self, u_node, omega):
        
        # initialise 6x1 empty complex vector 
        u_full = np.zeros(6,complex)
        
        # calculate local load vector
        u_local = self.LocalDisplacements(u_node, omega)
        
        # assign to correct locations in full vector
        for i, row in enumerate(self.nodal_dofs):
            u_full[row] = u_local[i]
        
        return u_full
    
    def FullElementDisplacements(self, element, u_nodes_global, omega, num_points):
        
        # intilise full u_elem
        u_elem_full = [None] * Element.Ndof
        
        # calculate local u_elem
        u_elem_local = element.LocalElementDisplacements(u_nodes_global, omega, num_points)
        
        # assign
        for i, dof in enumerate(self.dofs):
            u_elem_full[dof] = u_elem_local[i]
        
        return u_elem_full
    
                        
        
    
    