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


def apply_force(truss, force_vector, node_indices):
    """
    Modified force application function that distributes forces more evenly.
    
    Args:
        truss: The truss structure object
        force_vector: The force vector [Fx, Fy, Fz] to be applied
        node_indices: List of node indices where force should be applied
        
    Returns:
        F: Global force vector with distributed forces
    """
    num_nodes = len(truss.nodes)
    F = np.zeros(3 * num_nodes)
    
    # Distribute the force equally among all specified nodes
    num_application_points = len(node_indices)
    distributed_force = np.array(force_vector) / num_application_points
    
    for idx in node_indices:
        # Apply distributed force components to global force vector
        F[3*idx:3*idx+3] += distributed_force
    
    return F


def calculate_stiffness_matrix(truss, E):
    """
    Calculates the global stiffness matrix of the truss.

    Args:
        E: Young's modulus of the material.
        A: Cross-sectional area of the elements.

    Returns:
        The global stiffness matrix (K).
    """
    num_nodes = len(truss.nodes)
    K = np.zeros((3 * num_nodes, 3 * num_nodes))

    for i, element in enumerate(truss.elements):
        node1, node2 = element
        idx1, idx2 = truss.node_index[tuple(node1)], truss.node_index[tuple(node2)]

        x1, y1, z1 = node1
        x2, y2, z2 = node2
        L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        c_x, c_y, c_z = (x2 - x1) / L, (y2 - y1) / L, (z2 - z1) / L

        k_local = (E * truss.A[i] / L) * np.array([
            [c_x**2, c_x*c_y, c_x*c_z, -c_x**2, -c_x*c_y, -c_x*c_z],
            [c_x*c_y, c_y**2, c_y*c_z, -c_x*c_y, -c_y**2, -c_y*c_z],
            [c_x*c_z, c_y*c_z, c_z**2, -c_x*c_z, -c_y*c_z, -c_z**2],
            [-c_x**2, -c_x*c_y, -c_x*c_z, c_x**2, c_x*c_y, c_x*c_z],
            [-c_x*c_y, -c_y**2, -c_y*c_z, c_x*c_y, c_y**2, c_y*c_z],
            [-c_x*c_z, -c_y*c_z, -c_z**2, c_x*c_z, c_y*c_z, c_z**2]
        ])

        dof_indices = [3 * idx1, 3 * idx1 + 1, 3 * idx1 + 2, 3 * idx2, 3 * idx2 + 1, 3 * idx2 + 2]
        for i in range(6):
            for j in range(6):
                K[dof_indices[i], dof_indices[j]] += k_local[i, j]

    return K

def apply_boundary_conditions(truss, K, F):
    """
    Applies boundary conditions (fixed supports) to the stiffness matrix and force vector.

    Args:
        K: The global stiffness matrix.
        F: The global force vector.

    Returns:
        K_reduced: The reduced stiffness matrix after applying boundary conditions.
        F_reduced: The reduced force vector after applying boundary conditions.
    """
    # Fix the bottom nodes (all DOFs constrained)
    for i in range(len(truss.columns)):
        top_node = tuple(truss.columns[i][-1])
        idx = truss.node_index[top_node]
        for j in range(3):  # x, y, z directions
            dof = 3 * idx + j
            K[dof, :] = 0
            K[:, dof] = 0
            K[dof, dof] = 1
            F[dof] = 0

    return K, F


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
def calculate_mass_matrix(truss, spacecraft_mass, lumped_mass=False):
    """
    Calculates the global mass matrix of the truss with proper mass distribution.
    
    Args:
        truss: The truss structure object
        spacecraft_mass: Mass of the spacecraft to be added to top nodes
        lumped_mass: If True, use lumped mass matrix, else consistent mass
        
    Returns:
        M: Global mass matrix
    """
    num_nodes = len(truss.nodes)
    M = np.zeros((3 * num_nodes, 3 * num_nodes))
    rho = 2700  # density of aluminum (kg/m^3)

    # Calculate structural mass matrix
    for i, element in enumerate(truss.elements):
        node1, node2 = element
        idx1, idx2 = truss.node_index[tuple(node1)], truss.node_index[tuple(node2)]
        
        # Element length
        dx = node2[0] - node1[0]
        dy = node2[1] - node1[1]
        dz = node2[2] - node1[2]
        length = np.sqrt(dx**2 + dy**2 + dz**2)
        
        if lumped_mass:
            # Lumped mass matrix
            element_mass = rho * truss.A[i] * length / 2
            for idx in [idx1, idx2]:
                for j in range(3):
                    M[3*idx + j, 3*idx + j] += element_mass
        else:
            # Consistent mass matrix
            element_mass = rho * truss.A[i] * length
            m_local = (element_mass / 6) * np.array([
                [2, 0, 0, 1, 0, 0],
                [0, 2, 0, 0, 1, 0],
                [0, 0, 2, 0, 0, 1],
                [1, 0, 0, 2, 0, 0],
                [0, 1, 0, 0, 2, 0],
                [0, 0, 1, 0, 0, 2]
            ])
            
            # Transform to global coordinates if needed
            # (simplified here - add transformation if needed)
            
            # Assemble into global matrix
            dofs = [3*idx1, 3*idx1+1, 3*idx1+2, 3*idx2, 3*idx2+1, 3*idx2+2]
            for p in range(6):
                for q in range(6):
                    M[dofs[p], dofs[q]] += m_local[p, q]

    # Add spacecraft mass to top nodes only
    top_nodes = set(tuple(column[-1]) for column in truss.columns)
    mass_per_top_node = spacecraft_mass / len(top_nodes)
    
    for node in top_nodes:
        idx = truss.node_index[node]
        for j in range(3):
            M[3*idx + j, 3*idx + j] += mass_per_top_node

    return M

def solve_eigenvalue_problem(truss, E, spacecraft_mass, lumped_mass=False):
    """
    Solves the eigenvalue problem with improved stability.
    
    Args:
        truss: The truss structure object
        E: Young's modulus
        spacecraft_mass: Mass of the spacecraft
        lumped_mass: Whether to use lumped mass matrix
        
    Returns:
        eigenvalues: Natural frequencies squared
        eigenvectors: Mode shapes
        free_dofs: Array of free degrees of freedom
    """
    # Calculate matrices
    K = calculate_stiffness_matrix(truss, E)
    M = calculate_mass_matrix(truss, spacecraft_mass, lumped_mass)
    
    # Apply boundary conditions
    num_nodes = len(truss.nodes)
    fixed_dofs = []
    
    # Fix bottom nodes
    for column in truss.columns:
        bottom_node = tuple(column[0])
        idx = truss.node_index[bottom_node]
        fixed_dofs.extend([3*idx, 3*idx + 1, 3*idx + 2])
    
    fixed_dofs = np.array(fixed_dofs)
    all_dofs = np.arange(3 * num_nodes)
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)
    
    # Reduce matrices
    K_reduced = K[np.ix_(free_dofs, free_dofs)]
    M_reduced = M[np.ix_(free_dofs, free_dofs)]
    
    # Check matrix properties
    if not np.allclose(K_reduced, K_reduced.T):
        print("Warning: K_reduced is not symmetric")
    if not np.allclose(M_reduced, M_reduced.T):
        print("Warning: M_reduced is not symmetric")
    
    # Solve eigenvalue problem
    try:
        eigenvalues, eigenvectors = scipy.linalg.eigh(K_reduced, M_reduced)
        
        # Check for negative eigenvalues
        if np.any(eigenvalues < -1e-10):  # allowing for small numerical errors
            print("Warning: Negative eigenvalues found")
            
    except np.linalg.LinAlgError:
        print("Error in eigenvalue solution - check matrix conditioning")
        raise
        
    return eigenvalues, eigenvectors, free_dofs

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



def visualize_mode_shape(truss, eigenvectors, free_dofs,scale_factor=50):
    """
    Visualizes the mode shapes of the truss structure.
    
    Args:
        truss: The truss structure object
        eigenvectors: Eigenvectors from modal analysis
        scale_factor: Factor to scale the mode shape displacements
    """
    # Create subplots for the first few modes
    num_modes = min(6, eigenvectors.shape[1])  # Show up to 6 modes
    fig = plt.figure(figsize=(15, 10))
    
    for mode in range(num_modes):
        ax = fig.add_subplot(2, 3, mode + 1, projection='3d')
        
        # Get the mode shape vector
        mode_shape = np.zeros(len(truss.nodes) * 3)
        mode_shape[free_dofs] = eigenvectors[:, mode]
        
        # Calculate displaced positions
        displaced_nodes = []
        for i, node in enumerate(truss.nodes):
            orig_pos = np.array(node)
            # Get displacement for this node
            disp = mode_shape[i*3:(i+1)*3]
            # Apply scaled displacement
            new_pos = orig_pos + disp * scale_factor
            displaced_nodes.append(new_pos)
        
        displaced_nodes = np.array(displaced_nodes)
        
        # Plot original structure (grey)
        for element in truss.elements:
            x = [element[0][0], element[1][0]]
            y = [element[0][1], element[1][1]]
            z = [element[0][2], element[1][2]]
            ax.plot(x, y, z, 'gray', alpha=0.3)
        
        # Plot displaced structure (blue)
        for i, element in enumerate(truss.elements):
            idx1 = truss.node_index[tuple(element[0])]
            idx2 = truss.node_index[tuple(element[1])]
            
            x = [displaced_nodes[idx1, 0], displaced_nodes[idx2, 0]]
            y = [displaced_nodes[idx1, 1], displaced_nodes[idx2, 1]]
            z = [displaced_nodes[idx1, 2], displaced_nodes[idx2, 2]]
            ax.plot(x, y, z, 'b')
        
        # Plot nodes
        ax.scatter(displaced_nodes[:, 0], displaced_nodes[:, 1], 
                  displaced_nodes[:, 2], c='r', s=50)
        
        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'Mode {mode + 1}')
        
        # Set equal aspect ratio
        max_range = np.array([
            displaced_nodes[:, 0].max() - displaced_nodes[:, 0].min(),
            displaced_nodes[:, 1].max() - displaced_nodes[:, 1].min(),
            displaced_nodes[:, 2].max() - displaced_nodes[:, 2].min()
        ]).max() / 2.0
        
        mid_x = (displaced_nodes[:, 0].max() + displaced_nodes[:, 0].min()) * 0.5
        mid_y = (displaced_nodes[:, 1].max() + displaced_nodes[:, 1].min()) * 0.5
        mid_z = (displaced_nodes[:, 2].max() + displaced_nodes[:, 2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        # Adjust view angle for better visualization
        ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    plt.show()

def calculate_stiffness_matrix_calfem(truss, E):
    """Calculates the global stiffness matrix using CALFEM."""
    num_nodes = len(truss.nodes)
    K = np.zeros((3 * num_nodes, 3 * num_nodes))

    for i, element in enumerate(truss.elements):
        A = truss.A[i]
        node1, node2 = element
        idx1, idx2 = truss.node_index[tuple(node1)], truss.node_index[tuple(node2)]

        ex = [node1[0], node2[0]]
        ey = [node1[1], node2[1]]
        ez = [node1[2], node2[2]]
        ep = [E, A]

        # Calculate element stiffness matrix
        k_local = cfc.bar3e(ex, ey, ez, ep)

        # --- Debugging Prints ---
        

        # Define degrees of freedom for the element
        dof_indices = np.array([3*idx1, 3*idx1+1, 3*idx1+2, 3*idx2, 3*idx2+1, 3*idx2+2])

        # --- Debugging Prints ---
       

        # Add element stiffness to global stiffness matrix
        K = cfc.assem(dof_indices, K, k_local) 
    return K

def calculate_mass_matrix_calfem(truss, spacecraft_mass, lumped_mass=False):
    """Calculates the global mass matrix using CALFEM (lumped)."""
    num_nodes = len(truss.nodes)
    M = np.zeros((3 * num_nodes, 3 * num_nodes))

    for i, element in enumerate(truss.elements):
        node1, node2 = element
        idx1, idx2 = truss.node_index[tuple(node1)], truss.node_index[tuple(node2)]
        L = truss.length[i]
        A = truss.A[i]

        # Use lumped mass approximation for simplicity
        m_local = rho * A * L / 2 * np.eye(6)

        # Assemble into the global mass matrix
        dof_indices = np.array([3 * idx1, 3 * idx1 + 1, 3 * idx1 + 2, 3 * idx2, 3 * idx2 + 1, 3 * idx2 + 2])
        M = cfc.assem(dof_indices, M, m_local)

    # Add spacecraft mass to top nodes
    num_top_nodes = len(truss.columns[-1])
    top_node_mass = spacecraft_mass / num_top_nodes
    for column in truss.columns:
        top_node_index = truss.node_index[tuple(column[-1])]
        M[3 * top_node_index : 3 * top_node_index + 3, 3 * top_node_index : 3 * top_node_index + 3] += top_node_mass * np.eye(3)
    #print(f"Shape of M in calculate_mass_matrix: {M.shape}")
    return M

def apply_boundary_conditions2(truss, K, F):
    """
    Modified boundary conditions with improved stability.
    
    Args:
        truss: The truss structure object
        K: The global stiffness matrix
        F: The global force vector
        
    Returns:
        K_reduced: The reduced stiffness matrix
        F_reduced: The reduced force vector
        free_dofs: Array of free degrees of freedom
    """
    num_nodes = len(truss.nodes)
    fixed_dofs = []
    
    # Find and fix bottom nodes
    # for column in truss.columns:
    #     bottom_node = tuple(column[0])
    #     idx = truss.node_index[bottom_node]
    #     fixed_dofs.extend([3*idx, 3*idx + 1, 3*idx + 2])
    
    # fixed_dofs = np.array(fixed_dofs)
    # all_dofs = np.arange(3 * num_nodes)
    # free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    

    # Iterate through all columns to fix base nodes' X and Y DOFs
    for column in truss.columns:
        bottom_node = tuple(column[0])
        idx = truss.node_index.get(bottom_node)
        if idx is not None:
            fixed_dofs.extend([3*idx, 3*idx + 1])  # Fix X and Y DOFs only
        else:
            print(f"Warning: Node {bottom_node} not found in node_index.")

    # Remove any duplicate DOFs
    fixed_dofs = np.unique(fixed_dofs)

    # Define all DOFs
    num_nodes = len(truss.nodes)
    all_dofs = np.arange(3 * num_nodes)

    # Define free DOFs
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    print(f"Number of fixed DOFs: {len(fixed_dofs)}")
    print(f"Number of free DOFs: {len(free_dofs)}")
    print(f"Fixed DOFs indices: {fixed_dofs}")
                                
    # Apply matrix conditioning
    K_reduced = K[np.ix_(free_dofs, free_dofs)]
    F_reduced = F[free_dofs]
    
    # Improve numerical stability
    scale_factor = 1.0 / np.max(np.abs(K_reduced))
    K_reduced = K_reduced * scale_factor
    F_reduced = F_reduced * scale_factor
    
    return K_reduced, F_reduced, free_dofs

def expand_displacement_vector(u_reduced, free_dofs, total_dofs):
    """
    Expands the reduced displacement vector to full size.
    
    Args:
        u_reduced: Reduced displacement vector
        free_dofs: Array of free degrees of freedom
        total_dofs: Total number of degrees of freedom
        
    Returns:
        u_full: Full-size displacement vector
    """
    u_full = np.zeros(total_dofs)
    u_full[free_dofs] = u_reduced
    return u_full

def apply_boundary_conditions_eigenproblem(truss, K, M):
    """
    Applies boundary conditions for eigenvalue problem.
    
    Args:
        truss: The truss structure object
        K: The global stiffness matrix
        M: The global mass matrix
        
    Returns:
        K_reduced: The reduced stiffness matrix
        M_reduced: The reduced mass matrix
        free_dofs: Array of free degrees of freedom
    """
    # num_nodes = len(truss.nodes)
    # fixed_dofs = []
    # for i in range(len(truss.columns)):
    #     bottom_node = tuple(truss.columns[i][0])
    #     idx = truss.node_index[bottom_node]
    #     for j in range(3):  # x, y, z directions
    #         fixed_dofs.append(3 * idx + j)

    # all_dofs = np.arange(3 * num_nodes)
    # free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    # K_reduced = K[np.ix_(free_dofs, free_dofs)]
    # M_reduced = M[np.ix_(free_dofs, free_dofs)]
    
    # return K_reduced, M_reduced, free_dofs
    

    num_nodes = len(truss.nodes)
    fixed_dofs = []

    # Identify bottom nodes (assuming bottom nodes have the minimum z-coordinate)
    min_z = min(node[2] for node in truss.nodes)
    bottom_nodes = [i for i, node in enumerate(truss.nodes) if node[2] == min_z]

    # Fix degrees of freedom (DOFs) for bottom nodes in x, y, z directions
    for idx in bottom_nodes:
        for j in range(3):  # x, y, z directions
            fixed_dofs.append(3 * idx + j)

    # Determine free DOFs
    all_dofs = np.arange(3 * num_nodes)
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    # Reduce matrices to free DOFs
    K_reduced = K[np.ix_(free_dofs, free_dofs)]
    M_reduced = M[np.ix_(free_dofs, free_dofs)]

    return K_reduced, M_reduced, free_dofs



def solve_eigenvalue_problem2(truss, E, spacecraft_mass, lumped_mass=False):
    """
    Solves the eigenvalue problem to find natural frequencies and mode shapes.
    
    Args:
        truss: The truss structure object
        E: Young's modulus
        spacecraft_mass: The mass of the spacecraft
        lumped_mass: Whether to use lumped mass matrix
        
    Returns:
        eigenvalues, eigenvectors, and free_dofs
    """
    K = calculate_stiffness_matrix_calfem(truss, E)
    M = calculate_mass_matrix_calfem(truss, spacecraft_mass)

    # Apply boundary conditions using the new function
    K_reduced, M_reduced, free_dofs = apply_boundary_conditions_eigenproblem(truss, K, M)

    # Solve the eigenvalue problem
    eigenvalues, eigenvectors = eigh(K_reduced, M_reduced)
    return eigenvalues, eigenvectors, free_dofs


