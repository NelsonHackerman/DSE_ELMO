import numpy as np
import matplotlib.pyplot as plt
import scipy
from Truss_force_and_stiffness import apply_force, apply_boundary_conditions, calculate_stiffness_matrix, calculate_element_forces, visualize_truss, check_stress, calculate_mass, calculate_mass_matrix, solve_eigenvalue_problem, calculate_natural_frequencies, visualize_mode_shape 
global E,v,rho,Sig_tu,Sig_ty,g,SF,l_ax_com,l_ax_ten,l_lat,f_ax,f_lat
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,morb,mtot,l_ax_com,l_ax_ten,l_lat,m_bend_u,m_bend_l,l_eq_ten_u,l_eq_ten_l,l_eq_com_u,l_eq_com_l,f_ax,f_lat, ru,rl,lu,ll,dl,du



class TrussStructure:

    def __init__(self):
        pass

    def DefineGeometry(self, baseShape, radius, height, rings):
        self.columns = []
        ringPitch = height / (rings+1)

        if baseShape == 'square':
            self.baseNodes = [[radius, radius, 0], [-radius, radius, 0], [-radius, -radius, 0], [radius, -radius, 0]]
        else:
            basePitch = (2 * np.pi) / baseShape
            self.baseNodes = []
            for i in range(baseShape):
                baseNode = [radius * np.cos(basePitch*i), radius * np.sin(basePitch*i), 0]
                self.baseNodes.append(baseNode)
            
        for baseNode in self.baseNodes:
            column = []
            for i in range(rings+2):
                column.append([baseNode[0], baseNode[1], ringPitch*i])
            self.columns.append(column)
        self.nodes = [node for column in self.columns for node in column]
        self.elements = []
        self.t=[]
        self.r=[]
        self.A=[]
        # #vertical connections between nodes in the same column
        for column in self.columns:
            for i in range(len(column)-1):
                element = [column[i], column[i+1]]
                t=0.002
                r=0.06
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.r.append(r)
                self.A.append(A)
                self.elements.append(element)

        #horizonal connections between nodes on the same vertical level
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])-1):
                element = [self.columns[i][j], self.columns[(i+1) % (len(self.columns))][j]]
                self.elements.append(element)
                t=0.002
                r=0.06
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.r.append(r)
                self.A.append(A)
                
        
        # #diagonal connections
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])-1):
                element = [self.columns[i][j], self.columns[(i+1) % len(self.columns)][j+1]]
                self.elements.append(element)
                t=0.002
                r=0.06
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.r.append(r)
                self.A.append(A)
        
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])-1):
                element = [self.columns[i][j], self.columns[(i-1) % len(self.columns)][j+1]]
                self.elements.append(element)
                t=0.002
                r=0.06
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.r.append(r)
                self.A.append(A)
        
        self.length=[]
        for element in self.elements:
            node1, node2 = element
            length = np.sqrt((node2[0] - node1[0])**2 + (node2[1] - node1[1])**2 + (node2[2] - node1[2])**2)
            self.length.append(length)
        self.node_index = {tuple(node): i for i, node in enumerate(self.nodes)}


truss = TrussStructure()
truss.DefineGeometry(4, ru, lu, 1)
#print(truss.length)
#fig = plt.figure()
#ax = plt.axes(projection ='3d')


for element in truss.elements:
    x = [element[0][0], element[1][0]]
    y = [element[0][1], element[1][1]]
    z = [element[0][2], element[1][2]]
    #ax.plot3D(x, y, z, 'black')

#plt.show()

# --- External Force ---
force_vector = [0, 0, l_eq_com_u*morb]  # Example: 10 kN force in the negative z-direction
bottom_node_indices = [truss.node_index[tuple(np.float64(x) for x in column[0])] for column in truss.columns] # Apply at top nodes
F = apply_force(truss,force_vector, bottom_node_indices)

# --- Calculate Stiffness Matrix and Apply Boundary Conditions ---
K = calculate_stiffness_matrix(truss,E)
K, F = apply_boundary_conditions(truss,K, F)

# --- Solve for Displacements ---
u = np.linalg.solve(K, F)

# --- Calculate Element Forces ---
element_forces, element_stresses = calculate_element_forces(truss,E, u)
for i, force in enumerate(element_forces):
    print(f"Element {i+1} force: {force:.2f} N, stress factor : {element_stresses[i]/Sig_ty*SF:.2f} ")
stress_check=check_stress(element_stresses, Sig_ty/SF,truss)

# --- Visualize the Truss ---
visualize_truss(truss,stress_check,element_forces)
mass=calculate_mass(truss)
print('Mass: ',mass)
print(truss.nodes)

M=calculate_mass_matrix(truss)
eigenvalues,eigenvectors=solve_eigenvalue_problem(truss, E)
frequencies=calculate_natural_frequencies(truss, E)
print('frequencies ',frequencies)
visualize_mode_shape(truss, eigenvectors, scale_factor=50)
print(truss.nodes)