import numpy as np
import matplotlib.pyplot as plt
import scipy
#from Structures.Truss_structure import TrussStructure
global E,v,rho,Sig_tu,Sig_ty,g,SF,l_ax_com,l_ax_ten,l_lat,f_ax,f_lat
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,morb,mtot,l_ax_com,l_ax_ten,l_lat,m_bend_u,m_bend_l,l_eq_ten_u,l_eq_ten_l,l_eq_com_u,l_eq_com_l,f_ax,f_lat, ru,rl,lu,ll,dl,du


def apply_force(Truss, force_vector, node_indices):
    """
    Applies an external force to specified nodes.

    Args:
        force_vector: A list or array representing the force vector (Fx, Fy, Fz).
        node_indices: A list of node indices where the force is applied.
    """
    num_nodes = len(Truss.nodes)
    F = np.zeros(3 * num_nodes)

    for node_index in node_indices:
        idx = Truss.node_index[tuple(Truss.nodes[node_index])]
        F[3 * idx: 3 * idx + 3] = force_vector

    return F

def calculate_stiffness_matrix(Truss, E):
    """
    Calculates the global stiffness matrix of the truss.

    Args:
        E: Young's modulus of the material.
        A: Cross-sectional area of the elements.

    Returns:
        The global stiffness matrix (K).
    """
    num_nodes = len(Truss.nodes)
    K = np.zeros((3 * num_nodes, 3 * num_nodes))

    for i, element in enumerate(Truss.elements):
        node1, node2 = element
        idx1, idx2 = Truss.node_index[tuple(node1)], Truss.node_index[tuple(node2)]

        x1, y1, z1 = node1
        x2, y2, z2 = node2
        L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        c_x, c_y, c_z = (x2 - x1) / L, (y2 - y1) / L, (z2 - z1) / L

        k_local = (E * Truss.A[i] / L) * np.array([
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

def apply_boundary_conditions(Truss, K, F):
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
    for i in range(len(Truss.columns)):
        top_node = tuple(Truss.columns[i][-1])
        idx = Truss.node_index[top_node]
        for j in range(3):  # x, y, z directions
            dof = 3 * idx + j
            K[dof, :] = 0
            K[:, dof] = 0
            K[dof, dof] = 1
            F[dof] = 0

    return K, F

def calculate_element_forces(truss, E, u):
        """
        Calculates the internal forces and stresses in each element.

        Args:
            E: Young's modulus.
            A: Cross-sectional area.
            u: The displacement vector.

        Returns:
            A list of element forces (tension is positive, compression is negative).
            A list of element stresses.
        """
        element_forces = []
        element_stresses = []
        for i, element in enumerate(truss.elements):
            node1, node2 = element
            idx1, idx2 = truss.node_index[tuple(node1)], truss.node_index[tuple(node2)]

            x1, y1, z1 = node1
            x2, y2, z2 = node2
            L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            c_x, c_y, c_z = (x2 - x1) / L, (y2 - y1) / L, (z2 - z1) / L

            # Extract displacements for the element
            u_local = np.array([u[3 * idx1], u[3 * idx1 + 1], u[3 * idx1 + 2],
                                u[3 * idx2], u[3 * idx2 + 1], u[3 * idx2 + 2]])

            # Calculate element force
            force = (E * truss.A[i] / L) * np.array([-c_x, -c_y, -c_z, c_x, c_y, c_z]).dot(u_local)
            element_forces.append(force)

            # Calculate element stress (simplified)
            stress = force / (truss.A[i])
            element_stresses.append(stress)
            force_vector = force * np.array([c_x, c_y, c_z])  # Force vector along the element
            if not np.allclose(force_vector, force * np.array([c_x, c_y, c_z]), atol=1e-6):
                print(f"Warning: Force is not aligned with element for element connecting nodes {node1} and {node2}")
                print(f"  Force vector: {force_vector}")
                print(f"  Expected direction: {[c_x, c_y, c_z]}")

        return element_forces, element_stresses

def calculate_mass(truss):
    truss.mass=[]
    for i in range (len(truss.elements)):
        elem_mass = truss.length[i]*truss.A[i]*rho
        truss.mass.append(elem_mass)
    mass=np.sum(truss.mass)
    
    return mass


def check_stress(element_stresses, allowable_stress,truss):
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
        if stress > allowable_stress:
            stress_exceeded.append(True)
        
        elif -stress>.6*E*truss.t[i]/truss.r[i]*(1-0.901*(1-np.exp(-1/16*np.sqrt(truss.r[i]/truss.t[i])))):
            stress_exceeded.append(True)
            print('buckles')
            
        else:
            stress_exceeded.append(False)
    print(stress_exceeded)
    return stress_exceeded

def visualize_truss(Truss,stress_check, element_forces=None):
    """
    Visualizes the truss structure using matplotlib.

    Args:
        element_forces: (Optional) A list of element forces. If provided, zero-force
                        members will be colored red.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i, element in enumerate(Truss.elements):
        x = [element[0][0], element[1][0]]
        y = [element[0][1], element[1][1]]
        z = [element[0][2], element[1][2]]

        if element_forces is not None and abs(element_forces[i]) < 1e-6:
            color = 'b'  # Blue for zero-force members
        elif stress_check[i]==True:
            color = 'r'  # Red for above stress limit
        else:
            color = 'g'  # Green for below stress limit

        ax.plot(x, y, z, color=color)
        if element_forces is not None:
                # Calculate the midpoint of the element
                mid_x = (element[0][0] + element[1][0]) / 2
                mid_y = (element[0][1] + element[1][1]) / 2
                mid_z = (element[0][2] + element[1][2]) / 2

                # Add a text label near the midpoint
                label_text = f"E{i+1}"
                ax.text(mid_x, mid_y, mid_z, label_text, color='black', fontsize=8)
    # Plot nodes
    node_coords = np.array(Truss.nodes)
    ax.scatter(node_coords[:, 0], node_coords[:, 1], node_coords[:, 2], c='k')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Truss Structure')
    plt.show()
    
def calculate_mass_matrix(Truss, lumped_mass=False, spacecraft_mass=morb):
    """
    Calculates the global mass matrix of the truss.

    Args:
        lumped_mass: If True, use a lumped mass matrix. If False, use a consistent mass matrix.
        spacecraft_mass: The mass of the spacecraft to be added to the top nodes.

    Returns:
        The global mass matrix (M).
    """
    num_nodes = len(Truss.nodes)
    M = np.zeros((3 * num_nodes, 3 * num_nodes))

    for i, element in enumerate(Truss.elements):
        node1, node2 = element
        idx1, idx2 = Truss.node_index[tuple(node1)], Truss.node_index[tuple(node2)]

        if lumped_mass:
            # Lumped mass matrix (simpler, less accurate)
            m_local = rho * Truss.A[i] * Truss.length[i] / 2 * np.eye(6)
        else:
            # Consistent mass matrix (more accurate)
            m_local = (rho * Truss.A[i] * Truss.length[i] / 420) * np.array([
            [140, 0, 0, 70, 0, 0],
            [0, 140, 0, 0, 70, 0],
            [0, 0, 140, 0, 0, 70],
            [70, 0, 0, 140, 0, 0],
            [0, 70, 0, 0, 140, 0],
            [0, 0, 70, 0, 0, 140]
        ])

    dof_indices = [3 * idx1, 3 * idx1 + 1, 3 * idx1 + 2, 3 * idx2, 3 * idx2 + 1, 3 * idx2 + 2]
    for p in range(6):
        for q in range(6):
            M[dof_indices[p], dof_indices[q]] += m_local[p, q]

    # Add spacecraft mass to top nodes
    num_top_nodes = len(Truss.columns[0])
    top_node_mass = spacecraft_mass / num_top_nodes
    for column in Truss.columns:
        top_node_index = Truss.node_index[tuple(column[-1])]
        M[3 * top_node_index, 3 * top_node_index] += top_node_mass  # x-direction
        M[3 * top_node_index + 1, 3 * top_node_index + 1] += top_node_mass  # y-direction
        M[3 * top_node_index + 2, 3 * top_node_index + 2] += top_node_mass  # z-direction

    return M

def solve_eigenvalue_problem(truss, E, lumped_mass=False, spacecraft_mass=morb):
    """
    Solves the eigenvalue problem to find natural frequencies and mode shapes.

    Args:
        E: Young's modulus.
        lumped_mass: If True, use a lumped mass matrix. If False, use a consistent mass matrix.
        spacecraft_mass: The mass of the spacecraft.

    Returns:
        A tuple containing:
            eigenvalues: The eigenvalues (squared angular natural frequencies).
            eigenvectors: The corresponding eigenvectors (mode shapes).
    """
    K = calculate_stiffness_matrix(truss,E)
    M = calculate_mass_matrix(truss, lumped_mass, spacecraft_mass)

    # Apply boundary conditions (fix bottom nodes)
    num_nodes = len(truss.nodes)
    fixed_dofs = []
    for i in range(len(truss.columns)):
        bottom_node = tuple(truss.columns[i][0])
        idx = truss.node_index[bottom_node]
        for j in range(3):  # x, y, z directions
            fixed_dofs.append(3 * idx + j)

    all_dofs = np.arange(3 * num_nodes)
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    K_reduced = K[np.ix_(free_dofs, free_dofs)]
    M_reduced = M[np.ix_(free_dofs, free_dofs)]

    # Solve the eigenvalue problem
    eigenvalues, eigenvectors = scipy.linalg.eigh(K_reduced, M_reduced)

    return eigenvalues, eigenvectors

def calculate_natural_frequencies(truss, E, lumped_mass=False, spacecraft_mass=morb):
    """
    Calculates the natural frequencies of the truss.

    Args:
        E: Young's modulus.
        lumped_mass: If True, use a lumped mass matrix. If False, use a consistent mass matrix.
        spacecraft_mass: The mass of the spacecraft.

    Returns:
        A list of natural frequencies in Hz.
    """
    eigenvalues, _ = truss.solve_eigenvalue_problem(E, lumped_mass, spacecraft_mass)
    angular_frequencies = np.sqrt(np.maximum(eigenvalues, 0))  # Ensure non-negative
    frequencies = angular_frequencies / (2 * np.pi)
    return frequencies

def visualize_mode_shape(truss, mode_shape, scale_factor=50):
    """
    Visualizes a specific mode shape of the truss.

    Args:
        mode_shape: The eigenvector corresponding to the mode shape.
        scale_factor: A scaling factor to exaggerate the deformation.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the undeformed truss in grey
    for element in truss.elements:
        x = [element[0][0], element[1][0]]
        y = [element[0][1], element[1][1]]
        z = [element[0][2], element[1][2]]
        ax.plot(x, y, z, color='grey', alpha=0.5)

    # Deform the nodes based on the mode shape
    deformed_nodes = []
    for i, node in enumerate(truss.nodes):
        displacement = mode_shape[3*i:3*i+3]
        deformed_node = (np.array(node) + scale_factor * displacement).tolist()
        deformed_nodes.append(deformed_node)

    # Plot the deformed truss in blue
    for element in truss.elements:
        node1, node2 = element
        idx1 = truss.node_index[tuple(node1)]
        idx2 = truss.node_index[tuple(node2)]
        x = [deformed_nodes[idx1][0], deformed_nodes[idx2][0]]
        y = [deformed_nodes[idx1][1], deformed_nodes[idx2][1]]
        z = [deformed_nodes[idx1][2], deformed_nodes[idx2][2]]
        ax.plot(x, y, z, color='b')

        # Plot nodes
        node_coords = np.array(deformed_nodes)
        ax.scatter(node_coords[:, 0], node_coords[:, 1], node_coords[:, 2], c='k')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Truss Mode Shape')
        plt.show()


