# %% import dependencies

import numpy as np

# import the defines elements
from .elements import Rod_1D
from .elements import EB_Beam

# %% class definition
class Element:
    ne = 0
    
    def clear():
        Element.ne = 0
        
    def __init__ (self, nodes, assembler=None, Rod_1D=Rod_1D, EB_Beam=EB_Beam):
        
        # assign element classes through dependency injection 
        self.Rod_1D = Rod_1D
        self.EB_Beam = EB_Beam
                
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
        
        # Keep track of what kind of elements are included in the form of a dictionary to be able to use them later
        self.element_types = {}
        
        # add element to the assembler
        if assembler is not None:
            assembler.RegisterElement(self)

    def SetSection (self, element_type, props):
        '''
        This function serves to set the elements, we will have to think how we do this properly but this is an initial set-up
        
        '''
        # TODO - need to handle errors correctly when not the correct parameters are given
        
        # add L to props
        props['L'] = self.L
        
        try:
            # check for element and accordingly assign the properties 
            if element_type == 'Rod':
                # unpack props
                element = Rod_1D(**props)
            elif element_type == 'EB Beam':
                element = self.EB_Beam(**props)
        
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
            k_loc += element.FullStiffness(omega)
        
        # return the full global stiffness matrix
        k_glob = np.matmul ( np.matmul ( self.Rt, k_loc ), self.R )
                
        return k_glob

    def AddDistributedLoad(self, q, omega):
        '''
        function to add distributed load to the element. 
        
        q = [q_rod;
             q_EB_V;
             q_EB_M;]
        '''
        
        q_loc = np.zeros(6, complex)
        
        for i, q_el in enumerate(q):
            if q_el != 0:
                if i == 0:
                    q_loc += self.element_types['Rod'].FullDistributedLoad(q_el, omega)
                elif i == 1:                    
                    q_loc += self.element_types['EB Beam'].FullDistributedLoad(q_el, omega)
        
        # get global element load
        q_glob = np.matmul(self.Rt, q_loc)
        
        # TODO - think of changing the code such that we get a clear separation between omega. Thus that we can create a loop that will define all variables per omega. Would be useful
        # assign to nodes
        self.nodes[0].add_load( q_glob[0:3] )
        self.nodes[1].add_load( q_glob[3:6] )
        
    def displacements ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        w = C[3] * np.cos(beta_b * x) + C[2] * np.sin(beta_b * x) + C[0] * np.cosh(beta_b * x) + C[1] * np.sinh(beta_b * x) - q / EI / beta_b**4
       
        return w
    
    def rotations ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        dw_dx = -C[3] * beta_b * np.sin(beta_b * x) + C[2] * beta_b * np.cos(beta_b * x) + C[0] * beta_b * np.sinh(beta_b * x) + C[1] * beta_b * np.cosh(beta_b * x)
       
        return -dw_dx
    
    def bending_moments ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        d2w_dx2 = -C[3] * beta_b ** 2 * np.cos(beta_b * x) - C[2] * beta_b ** 2 * np.sin(beta_b * x) + C[0] * beta_b ** 2 * np.cosh(beta_b * x) + C[1] * beta_b ** 2 * np.sinh(beta_b * x)
        
        M = -EI * d2w_dx2

        return M
    
    def shear_forces ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        d3w_dx3 = C[3] * beta_b ** 3 * np.sin(beta_b * x) - C[2] * beta_b ** 3 * np.cos(beta_b * x) + C[0] * beta_b ** 3 * np.sinh(beta_b * x) + C[1] * beta_b ** 3 * np.cosh(beta_b * x)
        
        return -EI * d3w_dx3



class Element_old:
    ne = 0

    def clear():
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
        
        self.R  = R
        self.Rt = np.transpose(R)
        
        # increase element number total
        Element.ne += 1

    def set_section (self, props):
        
        if 'EA' in props:
            self.EA = props['EA']
        else:
            self.EA = 1.e20
            
        if 'ksi' in props: # MODIFIED
            self.ksi = props['ksi'] # MODIFIED
        else: # MODIFIED
            self.ksi = 0.01  # MODIFIED
            
        if 'rhoA' in props:  # MODIFIED
            self.rhoA = props['rhoA']  # MODIFIED
        else:  # MODIFIED
            self.rhoA = 1.e20  # MODIFIED
            
        if 'EI' in props:
            self.EI = props['EI']
        else:
            self.EI = 1.e20
            
        if 'omega' in props:  # MODIFIED
            self.omega = props['omega']  # MODIFIED
        else:   # MODIFIED
            self.omega = 1.e20  # MODIFIED
        
        if 'k_1_r' in props:  # MODIFIED
            self.k_1_r = props['k_1_r']  # MODIFIED
        else:   # MODIFIED
            self.k_1_r = 1.e20  # MODIFIED
            
        if 'k_2_r' in props:  # MODIFIED
            self.k_2_r = props['k_2_r']  # MODIFIED
        else:   # MODIFIED
            self.k_2_r = 1.e20  # MODIFIED
            
        if 'c_1_r' in props:  # MODIFIED
            self.c_1_r = props['c_1_r']  # MODIFIED
        else:   # MODIFIED
            self.c_1_r = 1.e20  # MODIFIED
            
        if 'c_2_r' in props:  # MODIFIED
            self.c_2_r = props['c_2_r']  # MODIFIED
        else:   # MODIFIED
            self.c_2_r = 1.e20  # MODIFIED
        
        if 'k_1_b' in props:  # MODIFIED
            self.k_1_b = props['k_1_b']  # MODIFIED
        else:   # MODIFIED
            self.k_1_b = 1.e20  # MODIFIED
            
        if 'k_2_b' in props:  # MODIFIED
            self.k_2_b = props['k_2_b']  # MODIFIED
        else:   # MODIFIED
            self.k_2_b = 1.e20  # MODIFIED
            
        if 'c_1_b' in props:  # MODIFIED
            self.c_1_b = props['c_1_b']  # MODIFIED
        else:   # MODIFIED
            self.c_1_b = 1.e20  # MODIFIED
            
        if 'c_2_b' in props:  # MODIFIED
            self.c_2_b = props['c_2_b']  # MODIFIED
        else:   # MODIFIED
            self.c_2_b = 1.e20  # MODIFIED
            

    def global_dofs  (self):
        return np.hstack ((self.nodes[0].dofs, self.nodes[1].dofs))

    def stiffness ( self ):
        
        k = np.zeros ((6, 6), dtype=complex) # MODIFIED
        
        ksi = self.ksi # MODIFIED
        EA = self.EA * (1 + 2j * ksi) # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        c_r = (EA/rhoA) ** 0.5  # MODIFIED
        beta_r = omega / c_r  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        k_1_r = self.k_1_r # MODIFIED
        k_2_r = self.k_2_r # MODIFIED
        c_1_r = self.c_1_r # MODIFIED
        c_2_r = self.c_2_r # MODIFIED
        k_1_b = self.k_1_b # MODIFIED
        k_2_b = self.k_2_b # MODIFIED
        c_1_b = self.c_1_b # MODIFIED
        c_2_b = self.c_2_b # MODIFIED

        # Extension contribution

        k[0,0] = EA * beta_r * np.cos(beta_r * L) / np.sin(beta_r * L) + k_1_r + 1j * omega * c_1_r   # MODIFIED
        k[3,3] = EA * beta_r * np.cos(beta_r * L) / np.sin(beta_r * L) + k_2_r + 1j * omega * c_2_r  # MODIFIED
        k[3,0] = - EA * beta_r / np.sin(beta_r * L) # MODIFIED
        k[0,3] = - EA * beta_r / np.sin(beta_r * L) # MODIFIED

        # Bending contribution
        
        K_beam = np.array([[-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 3 * (np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b ** 2 * (-np.cosh(beta_b * L) + np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * beta_b ** 3 * (np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)],[EI * beta_b ** 2 * (np.cosh(beta_b * L) - np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.sinh(beta_b * L) + np.sin(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),-EI * np.sinh(beta_b * L) * np.sin(beta_b * L) * beta_b ** 2 / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1),EI * beta_b * (-np.cosh(beta_b * L) * np.sin(beta_b * L) + np.sinh(beta_b * L) * np.cos(beta_b * L)) / (np.cosh(beta_b * L) * np.cos(beta_b * L) - 1)]])

        k[1,1] = K_beam[0,0] + (k_1_b + 1j*omega*c_1_b)
        k[1,2] = K_beam[0,1]
        k[1,4] = K_beam[0,2]
        k[1,5] = K_beam[0,3]
        k[2,1] = K_beam[1,0]
        k[2,2] = K_beam[1,1]
        k[2,4] = K_beam[1,2]
        k[2,5] = K_beam[1,3]
        k[4,1] = K_beam[2,0]
        k[4,2] = K_beam[2,1]
        k[4,4] = K_beam[2,2] + (k_2_b + 1j*omega*c_2_b)
        k[4,5] = K_beam[2,3]
        k[5,1] = K_beam[3,0]
        k[5,2] = K_beam[3,1]
        k[5,4] = K_beam[3,2]
        k[5,5] = K_beam[3,3]
        
        #k[1,1] = k[4,4] =  12.0 * EI / L / L / L
        #k[1,4] = k[4,1] = -12.0 * EI / L / L / L
        #k[1,2] = k[2,1] = k[1,5] = k[5,1] = -6.0 * EI / L / L
        #k[2,4] = k[4,2] = k[4,5] = k[5,4] = 6.0 * EI / L / L
        #k[2,2] = k[5,5] = 4.0 * EI / L
        #k[2,5] = k[5,2] = 2.0 * EI / L

        return np.matmul ( np.matmul ( self.Rt, k ), self.R )

    def add_distributed_load ( self, q ):

        L = self.L

        self.q = np.array( q )
        # MODIFIED TO INCLUDE DYNAMIC LOAD UNIFORM DISTRIBUTED:
        el = np.array([(np.cos(beta_r*L) - 1.0)*q[0]/(np.sin(beta_r*L)*beta_r), q[1]*((np.cos(beta_b*L) - 1.0)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)) , -q[1]*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)) , (np.cos(beta_r*L) - 1.0)*q[0]/(np.sin(beta_r*L)*beta_r) , q[1]*((np.cos(beta_b*L) - 1)*np.sinh(beta_b*L) + np.cosh(beta_b*L)*np.sin(beta_b*L) - np.sin(beta_b*L))/(beta_b*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)) , q[1]*(np.sinh(beta_b*L)*np.sin(beta_b*L) - np.cosh(beta_b*L) + np.cos(beta_b*L))/(beta_b**2.0*(np.cosh(beta_b*L)*np.cos(beta_b*L) - 1.0)) ])
        
       # #el = [ 0.5*q[0]*l, 0.5*q[1]*l, -1.0/12.0*q[1]*l*l, 0.5*q[0]*l, 0.5*q[1]*l, 1.0/12.0*q[1]*l*l ]
    
        eg = np.matmul ( self.Rt, el)

        self.nodes[0].add_load ( eg[0:3] )
        self.nodes[1].add_load ( eg[3:6] )
        
    def displacements ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        w = C[3] * np.cos(beta_b * x) + C[2] * np.sin(beta_b * x) + C[0] * np.cosh(beta_b * x) + C[1] * np.sinh(beta_b * x) - q / EI / beta_b**4
       
        return w
    
    def rotations ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        dw_dx = -C[3] * beta_b * np.sin(beta_b * x) + C[2] * beta_b * np.cos(beta_b * x) + C[0] * beta_b * np.sinh(beta_b * x) + C[1] * beta_b * np.cosh(beta_b * x)
       
        return -dw_dx
    
    def bending_moments ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        d2w_dx2 = -C[3] * beta_b ** 2 * np.cos(beta_b * x) - C[2] * beta_b ** 2 * np.sin(beta_b * x) + C[0] * beta_b ** 2 * np.cosh(beta_b * x) + C[1] * beta_b ** 2 * np.sinh(beta_b * x)
        
        M = -EI * d2w_dx2

        return M
    
    def shear_forces ( self, u_global, num_points=2 ):
        
        ksi = self.ksi # MODIFIED
        rhoA = self.rhoA  # MODIFIED 
        EI = self.EI * (1 + 2j * ksi) # MODIFIED
        L = self.L
        omega = self.omega  # MODIFIED
        beta_b = (omega**2 * rhoA / EI) ** 0.25  # MODIFIED
        q = self.q[1]

        x = np.linspace ( 0.0, L, num_points )
        M  = np.zeros(num_points)

        ul = np.matmul ( self.R, u_global )
        
        ul = np.concatenate((ul[1:3],ul[4:]))
        
        A = np.array([[(np.sin(beta_b * L) * np.sinh(beta_b * L) + np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.cos(beta_b * L) - np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) - np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(-np.cos(beta_b * L) * np.sinh(beta_b * L) - np.sin(beta_b * L) * np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(1 + np.sin(beta_b * L) * np.sinh(beta_b * L) - np.cos(beta_b * L) * np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(np.sin(beta_b * L) + np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.sin(beta_b * L) * np.cosh(beta_b * L) + np.cos(beta_b * L) * np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) + 1) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.sin(beta_b * L) - np.sinh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.cos(beta_b * L) - np.cosh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2],[(np.cos(beta_b * L) * np.cosh(beta_b * L) - np.sin(beta_b * L) * np.sinh(beta_b * L) - 1) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(np.sin(beta_b * L) * np.cosh(beta_b * L) - np.cos(beta_b * L) * np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2,(-np.cos(beta_b * L) + np.cosh(beta_b * L)) / (2 * np.cos(beta_b * L) * np.cosh(beta_b * L) - 2),(-np.sin(beta_b * L) + np.sinh(beta_b * L)) / beta_b / (np.cos(beta_b * L) * np.cosh(beta_b * L) - 1) / 2]])
        
        C = A @ (ul + np.array([1/(EI*beta_b**4),0,1/(EI*beta_b**4),0]) * q)
        
        d3w_dx3 = C[3] * beta_b ** 3 * np.sin(beta_b * x) - C[2] * beta_b ** 3 * np.cos(beta_b * x) + C[0] * beta_b ** 3 * np.sinh(beta_b * x) + C[1] * beta_b ** 3 * np.cosh(beta_b * x)
        
        return -EI * d3w_dx3