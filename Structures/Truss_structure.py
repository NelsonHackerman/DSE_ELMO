import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.linalg import eigh
from numpy.linalg import solve
from Truss_force_and_stiffness import calculate_element_forces, visualize_truss, check_stress, calculate_mass, calculate_natural_frequencies
global E,v,rho,Sig_tu,Sig_ty,g,SF,l_ax_com,l_ax_ten,l_lat,f_ax,f_lat,Sig_c
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,morb,mtot,l_ax_com,l_ax_ten,l_lat,m_bend_u,m_bend_l,l_eq_ten_u,l_eq_ten_l,l_eq_com_u,l_eq_com_l,f_ax,f_lat, ru,rl,lu,ll,dl,du,Sig_c



class TrussStructure:

    def __init__(self):
        pass
    def DefineGeometry(self, baseShape, radius, height, rings,vt,vr,ht1,hr1,ht2,hr2,dt,dr):
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
        self.Avert=[]
        # #vertical connections between nodes in the same column
        for column in self.columns:
            for i in range(len(column)-1):
                element = [column[i], column[i+1]]
                t=vt
                r=vr
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.r.append(r)
                self.A.append(A)
                self.elements.append(element)
                if i==0:
                    self.Avert.append(A)

        #horizonal connections between nodes on the same vertical level
        for i in range(len(self.columns)):
            for j in range(1,len(self.columns[0])):
                element = [self.columns[i][j], self.columns[(i+1) % (len(self.columns))][j]]
                self.elements.append(element)
                if j==0:
                    t=ht1
                    r=hr1
                else:
                    t=ht2
                    r=hr2
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.r.append(r)
                self.A.append(A)
                
        
        # #diagonal connections
        intersection_nodes = []
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])-1):
                # Calculate intersection point of diagonals
                point1 = self.columns[i][j]
                point2 = self.columns[(i+1) % len(self.columns)][j+1]
                point3 = self.columns[i][j]
                point4 = self.columns[(i-1) % len(self.columns)][j+1]
                
                # Calculate midpoint of both diagonals (they intersect at the same point)
                mid_x = (point1[0] + point2[0]) / 2
                mid_y = (point1[1] + point2[1]) / 2
                mid_z = (point1[2] + point2[2]) / 2
                
                intersection_nodes.append([mid_x, mid_y, mid_z])

        # Add intersection nodes to main node list
        self.nodes.extend(intersection_nodes)
        print('nodes ',len(self.nodes))
        intersection_idx = 0
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])-1):
                start = self.columns[i][j]
                end = self.columns[(i+1) % len(self.columns)][j+1]
                mid_point = [(start[0] + end[0])/2, (start[1] + end[1])/2, (start[2] + end[2])/2]
                

                # Replace one diagonal with two elements meeting at midpoint
                self.elements.append([start, mid_point])
                self.elements.append([mid_point, end])
                t=dt
                r=dr
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.t.append(t)
                self.r.append(r)
                self.r.append(r)
                self.A.append(A)
                self.A.append(A)
        
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])-1):
                start = self.columns[i][j]
                end = self.columns[(i-1) % len(self.columns)][j+1]
                mid_point = [(start[0] + end[0])/2, (start[1] + end[1])/2, (start[2] + end[2])/2]
                

                # Replace one diagonal with two elements meeting at midpoint
                self.elements.append([start, mid_point])
                self.elements.append([mid_point, end])
                t=dt
                r=dr
                A=np.pi*(r**2-(r-t)**2)
                self.t.append(t)
                self.t.append(t)
                self.r.append(r)
                self.r.append(r)
                self.A.append(A)
                self.A.append(A)
        
        self.length=[]
        for element in self.elements:
            node1, node2 = element
            length = np.sqrt((node2[0] - node1[0])**2 + (node2[1] - node1[1])**2 + (node2[2] - node1[2])**2)
            self.length.append(length)
        self.node_index = {tuple(node): i for i, node in enumerate(self.nodes)}


truss = TrussStructure()
truss.DefineGeometry(4, ru, lu, 0,0.002,0.04,0.001,0.01,0.001,0.03,0.001,0.04)


for element in truss.elements:
    x = [element[0][0], element[1][0]]
    y = [element[0][1], element[1][1]]
    z = [element[0][2], element[1][2]]


# --- External Force ---

force_vector = [0, 0, -l_eq_com_u*morb] #-2812*9.81*(7.1+5.1*2.5)] 
force_vector2=[0, 0, -l_eq_ten_u*morb]
print('force vector: ',force_vector)
# --- Calculate Element Forces ---
element_forces, element_stresses = calculate_element_forces(truss,force_vector)
for i, force in enumerate(element_forces):
    print(f"Element {i+1} force: {force:.2f} N, stress factor : {element_stresses[i]/Sig_ty*SF:.2f} ")
stress_check=check_stress(element_stresses,truss)

# --- Visualize the Truss ---
visualize_truss(truss,stress_check,element_forces)
mass=calculate_mass(truss)
print('Mass: ',mass)

frequencies=calculate_natural_frequencies(truss, E, morb,lu)
print('frequencies ',frequencies)
totalmass=mass*1/0.6347
print('total mass: ',totalmass)
#element_forces, element_stresses = calculate_element_forces(truss,force_vector2)
#for i, force in enumerate(element_forces):
    #print(f"Element {i+1} force: {force:.2f} N, stress factor : {element_stresses[i]/Sig_ty*SF:.2f} ")
#stress_check=check_stress(element_stresses,truss)







truss2 = TrussStructure()
truss2.DefineGeometry(8, rl, ll, 1,0.003,0.03,0.001,0.03,0.001,0.03,0.002,0.03)


for element in truss2.elements:
    x = [element[0][0], element[1][0]]
    y = [element[0][1], element[1][1]]
    z = [element[0][2], element[1][2]]


# --- External Force ---
force_vector = [0, 0, -l_eq_com_u*mtot]  
print('force vector: ',force_vector)
# --- Calculate Element Forces ---
element_forces, element_stresses = calculate_element_forces(truss2,force_vector)
for i, force in enumerate(element_forces):
    print(f"Element {i+1} force: {force:.2f} N, stress factor : {element_stresses[i]/Sig_ty*SF:.2f} ")
stress_check=check_stress(element_stresses,truss2)

# --- Visualize the Truss ---
visualize_truss(truss2,stress_check,element_forces)
mass=calculate_mass(truss2)
print('Mass: ',mass)

frequencies=calculate_natural_frequencies(truss2, E, mtot,ll)
print('frequencies ',frequencies)
totalmass=mass*1/0.6347
print('total mass: ',totalmass)

