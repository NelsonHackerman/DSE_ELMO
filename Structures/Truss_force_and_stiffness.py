import numpy as np
import matplotlib.pyplot as plt
import scipy
import calfem.core as cfc
import calfem.utils as cfu
import calfem.vis_mpl as cfv
from scipy.linalg import eigh
from numpy.linalg import solve
#from Structures.Truss_structure import TrussStructure
global E,v,rho,Sig_tu,Sig_ty,g,SF,l_ax_com,l_ax_ten,l_lat,f_ax,f_lat,Sig_c
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,morb,mtot,l_ax_com,l_ax_ten,l_lat,m_bend_u,m_bend_l,l_eq_ten_u,l_eq_ten_l,l_eq_com_u,l_eq_com_l,f_ax,f_lat, ru,rl,lu,ll,dl,du,Sig_c


def calculate_element_forces(truss, external_force):
    """
    Calculate forces in truss elements using static equilibrium method
    
    Parameters:
    truss: TrussStructure object containing geometry and properties
    external_force: [Fx, Fy, Fz] external force vector applied to each top node
    
    Returns:
    element_forces: Array of axial forces in each element
    element_stresses: Array of axial stresses in each element
    """
    n_nodes = len(truss.nodes)
    n_elements = len(truss.elements)
    for i in range(len(external_force)):
        external_force[i]=external_force[i]/len(truss.baseNodes)
    
    
    # Create coefficient matrix A and force vector b for Ax = b
    # Each row represents equilibrium equation (3 per node: Fx, Fy, Fz)
    # Each column represents one element force
    A = np.zeros((3 * n_nodes, n_elements))
    b = np.zeros(3 * n_nodes)
    
    # Fill coefficient matrix A
    for i, node in enumerate(truss.nodes):
        node_tuple = tuple(node)
        # Find all elements connected to this node
        for j, element in enumerate(truss.elements):
            if list(node_tuple) == element[0] or list(node_tuple) == element[1]:
                # Calculate direction cosines
                node1 = element[0]
                node2 = element[1]
                dx = node2[0] - node1[0]
                dy = node2[1] - node1[1]
                dz = node2[2] - node1[2]
                L = np.sqrt(dx**2 + dy**2 + dz**2)
                
                # Direction depends on whether the node is start or end of element
                sign = -1 if list(node_tuple) == element[0] else 1
                
                # Add direction cosines to coefficient matrix
                A[3*i][j] = sign * dx/L      # x-component
                A[3*i+1][j] = sign * dy/L    # y-component
                A[3*i+2][j] = sign * dz/L    # z-component
    
    # Fill force vector b with external forces
    top_z = max(node[2] for node in truss.nodes)
    for i, node in enumerate(truss.nodes):
        if abs(node[2] - top_z) < 1e-10:  # Node is at top level
            b[3*i:3*i+3] = external_force
    
    # Remove equations for support nodes (base nodes)
    free_rows = []
    for i, node in enumerate(truss.nodes):
        if abs(node[2]) > 1e-10:  # Not a base node
            free_rows.extend([3*i, 3*i+1, 3*i+2])
    
    A_red = A[free_rows, :]
    b_red = b[free_rows]
    
    # Solve system of equations using least squares to handle potential redundancy
    element_forces = np.linalg.lstsq(A_red, b_red, rcond=None)[0]
    
    # Calculate stresses
    element_stresses = element_forces / np.array(truss.A)
    
    return element_forces, element_stresses



def calculate_mass(truss):
    truss.mass=[]
    for i in range (len(truss.elements)):
        elem_mass = truss.length[i]*truss.A[i]*rho
        truss.mass.append(elem_mass)
    mass=np.sum(truss.mass)
    
    return mass


def check_stress(element_stresses,truss):
    """
    Checks if the stress in each element is within the allowable limit.

    Args:
        element_stresses: A list of element stresses.
        allowable_stress: The allowable stress for the material.

    Returns:
        A list of booleans, where True indicates that the stress in the
        corresponding element exceeds the allowable limit.
    """
    stress_exceeded = []
    for i,stress in enumerate(element_stresses):
        if stress > Sig_tu/1.25:
            stress_exceeded.append(True)
            print('Tension exceeded')
        elif -stress>Sig_c:
            stress_exceeded.append(True)
            print('Compression exceeded')
        
        elif -stress>.6*E*truss.t[i]/truss.r[i]*(1-0.901*(1-np.exp(-1/16*np.sqrt(truss.r[i]/truss.t[i])))):
            stress_exceeded.append(True)
            print('buckles')
            
        else:
            stress_exceeded.append(False)
    #print(stress_exceeded)
    return stress_exceeded

def visualize_truss(truss, stress_check, element_forces=None):
    """
    Visualizes the truss structure with stress indicators and element labels.
    
    Args:
        truss: The truss structure object
        stress_check: List of boolean values indicating stress limit violations
        element_forces: List of element forces
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot elements
    for i, element in enumerate(truss.elements):
        x = [element[0][0], element[1][0]]
        y = [element[0][1], element[1][1]]
        z = [element[0][2], element[1][2]]
        
        # Determine color based on stress and force
        if element_forces is not None and abs(element_forces[i]) < 1e-6:
            color = 'b'  # Blue for zero-force members
        elif stress_check[i]:
            color = 'r'  # Red for above stress limit
        else:
            color = 'g'  # Green for below stress limit
        
        # Plot element
        ax.plot(x, y, z, color=color, linewidth=2)
        
        # Add element labels
        if element_forces is not None:
            mid_x = (element[0][0] + element[1][0]) / 2
            mid_y = (element[0][1] + element[1][1]) / 2
            mid_z = (element[0][2] + element[1][2]) / 2
            
            label_text = f"E{i+1}"
            ax.text(mid_x, mid_y, mid_z, label_text, color='black', fontsize=8)
    
    # Plot nodes
    node_coords = np.array(truss.nodes)
    ax.scatter(node_coords[:, 0], node_coords[:, 1], node_coords[:, 2], 
              c='k', s=50)
    
    # Set equal aspect ratio
    max_range = np.array([
        node_coords[:, 0].max() - node_coords[:, 0].min(),
        node_coords[:, 1].max() - node_coords[:, 1].min(),
        node_coords[:, 2].max() - node_coords[:, 2].min()
    ]).max() / 2.0
    
    mid_x = (node_coords[:, 0].max() + node_coords[:, 0].min()) * 0.5
    mid_y = (node_coords[:, 1].max() + node_coords[:, 1].min()) * 0.5
    mid_z = (node_coords[:, 2].max() + node_coords[:, 2].min()) * 0.5
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # Labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Truss Structure')
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='b', label='Zero Force'),
        Line2D([0], [0], color='r', label='Above Stress Limit'),
        Line2D([0], [0], color='g', label='Below Stress Limit')
    ]
    ax.legend(handles=legend_elements)
    
    # Adjust view angle
    ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    plt.show()


def calculate_natural_frequencies(truss, E, spacecraft_mass,l):

    """
    Retrieves cross-sectional areas of horizontal elements at a given level.

    Args:
        truss: A TrussStructure object.
        level: The level (ring index, starting from 0) to query.

    Returns:
        A list of cross-sectional areas (floats) for horizontal elements at that level.
    """
    areas = []
    frequencies=[]
    rs=[]
    num_columns = len(truss.columns)
    num_nodes_per_level = len(truss.columns[0])

    # Initialize areas and radii for vertical and diagonal beams
    vertical_areas = []
    vertical_radii = []
    vertical_elements=[]
    diagonal_areas = []
    diagonal_radii = []
    diagonal_elements=[]

    # Iterate through the elements in the first level only
    for i in range(num_columns):  # Iterate over columns
        for j in range(1):  # First level only
            # Get the current element
            element = truss.elements[j]
            node1, node2 = element

            # Extract node coordinates
            x1, y1, z1 = node1[0],node1[1],node1[2]
            x2, y2, z2 = node2[0],node2[1],node2[2]

            # Calculate the angle of the element with respect to the horizontal
            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1

            angle = np.arctan2(np.sqrt(dy**2 + dz**2), dx)  # 3D angle (dx is horizontal reference)

            # Check for vertical (90 degrees or π/2 radians)
            if np.isclose(angle, np.pi / 2):
                vertical_areas.append(truss.A[j])
                vertical_radii.append(truss.r[j])
                vertical_elements.append(truss.elements[j])

            # Check for diagonal (not horizontal or vertical)
            elif not np.isclose(angle, 0) and not np.isclose(angle, np.pi / 2):
                diagonal_areas.append(truss.A[j]/np.cos(angle))
                diagonal_radii.append(truss.r[j])
                diagonal_elements.append(truss.elements[j])

    # Output cross-sectional areas
    vertical_area_total = np.sum(vertical_areas)
    diagonal_area_total = np.sum(diagonal_areas)

    print(f"Vertical Beam Area (First Level): {vertical_area_total}")
    print(f"Diagonal Beam Area (First Level): {diagonal_area_total}")

    area=vertical_area_total+diagonal_area_total
    
    f_axx=0.25*np.sqrt(area*E/(spacecraft_mass*l))
    frequencies.append(f_axx)
    structure_center = np.mean(truss.nodes, axis=0)
    
    #Calculate the center of the element
    
    
    #for i in range(len(truss.Avert)):
        #element_center = np.mean(truss.elements, axis=0)
        #d=np.linalg.norm(np.array(element_center) - np.array(structure_center))
        #iss=truss.Avert[i]*d**2
        #iis.append(iss)
    #print('no nodes ',len(truss.Avert))
    iis=[]
    for i in range (len(vertical_areas)):
        element_center = np.mean(vertical_elements[i], axis=0)
        d=np.linalg.norm(np.array(element_center) - np.array(structure_center))
        iis.append(vertical_areas[i]*d**2)
    for i in range (len(diagonal_areas)):
        element_center = np.mean(diagonal_elements[i], axis=0)
        d=np.linalg.norm(np.array(element_center) - np.array(structure_center))
        iis.append(diagonal_areas[i]*d**2)
        
    I=np.sum(iis)
    f_latt=0.56*np.sqrt(I*E/(spacecraft_mass*l**3))
    frequencies.append(f_latt)
    return frequencies