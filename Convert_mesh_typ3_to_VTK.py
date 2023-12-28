# -*- coding: utf-8 -*-

#convertion des maillages 3D pour le FVCA6
#ce code permet de convertir un maillage de type typ3 ou msh en vtk
#le fichier vtk de sortie est lisible sur paraview 
#---------------------------------------------------------------------------------------------

import vtk
import numpy as np

#-----------Function to check if the nodes are oriented counterclockwise (anti-horaire)-------
def is_counterclockwise(node_indices, coords):
    a = np.array(coords[node_indices[0] - 1])  # Subtract 1 from the node index.
    b = np.array(coords[node_indices[1] - 1])  # Subtract 1 from the node index.
    c = np.array(coords[node_indices[2] - 1])  # Subtract 1 from the node index.
    cross_product = np.cross(b - a, c - a)
    return cross_product[2] < 0  # Z-component of the cross product (produit vectoriel)

filename = "dkershaw_1.msh"
output_filename = filename.replace(".msh", ".vtk")
print("Converting ", filename)

read_vertices = False
coords = []
read_cells = False
read_new_cell_connectivity = False
cell_connectivity = []
cell_codes = []

#----------------------read the mesh file typ3 of msh-----------------------------------------------
with open(filename, 'r', encoding='ISO-8859-1') as f:
    for line in f:
        # remove end of line character
        line = line[:-1]
        infos = line.split()

        if infos and infos[0] == "Vertices":
            nb_vertices = int(infos[1])
            read_vertices = True
        elif infos and (infos[0] == "Volumes->faces" or infos[0] == "Volumes->Faces"):
            # stop reading node coords
            read_vertices = False
        elif read_vertices:
            # read node coords
            coords_i = [float(v) for v in infos]
            coords.append(coords_i)  # Add the coordinates to the coords list

        elif infos and (infos[0] == "Volumes->Vertices" or infos[0] == "Volumes->Verticess"):
            # start reading cells connectivity
            nb_cells = int(infos[1])
            read_cells = True
            read_new_cell_connectivity = True

        elif infos and (infos[0] == "Faces->Edgess" or infos[0] == "Faces->Edges"):
            # stop reading cells connectivity
            read_cells = False
            read_faces = True
            nb_faces = int(infos[1])

        elif read_cells:
            values = [int(v) for v in infos]
            if read_new_cell_connectivity:
                nb_nodes_in_cell = values[0]
                cell_connectivity = values[1:]
                if len(cell_connectivity) < nb_nodes_in_cell:
                    read_new_cell_connectivity = False
                else:
                    read_new_cell_connectivity = True
            else:
                cell_connectivity += values
                if len(cell_connectivity) == nb_nodes_in_cell:
                    read_new_cell_connectivity = True
                else:
                    read_new_cell_connectivity = False

        if read_new_cell_connectivity and cell_connectivity:
            # Add the cell connectivity to the list of cells
            cell_codes.append(cell_connectivity)
            cell_connectivity = []

# Inside the loop that processes cell codes
for cell_nodes in cell_codes:
    if len(cell_nodes) > 3:  # Assuming faces with more than 3 points
        if not is_counterclockwise(cell_nodes, coords):
            cell_nodes.reverse()

# --------------Write the VTK file-----------------------------------------------------------------
with open(output_filename, 'w') as f:
    f.write('# vtk DataFile Version 2.0\n')
    f.write('Really cool data\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS {0} float\n'.format(nb_vertices))
    for v in range(nb_vertices):
        f.write('{0} {1} {2}\n'.format(coords[v][0], coords[v][1], coords[v][2]))

    #---- Write the cells----------------------------------------------------------------------------
    total_cell_indices = sum(len(cell) + 1 for cell in cell_codes)  # Add 1 for the cell size
    f.write('CELLS {0} {1}\n'.format(nb_cells, total_cell_indices))
    for cell in cell_codes:
        sorted_cell = cell[:]  # Make a copy of the cell to sort
        if len(sorted_cell) > 3 and not is_counterclockwise(sorted_cell, coords):
            sorted_cell.reverse()  # Reverse the cell if it's clockwise
        f.write('{0} '.format(len(sorted_cell)))
        f.write(' '.join(str(node - 1) for node in sorted_cell))  # Correct indexing to start from 0
        f.write('\n')

    # ---Write the cell types-------------------------------------------------------------------------
    f.write('CELL_TYPES {0}\n'.format(nb_cells))
    for nb_nodes_in_cell in [len(cell) for cell in cell_codes]:
        if nb_nodes_in_cell == 4:
            cell_type = vtk.VTK_TETRA
        elif nb_nodes_in_cell == 5:
            cell_type = vtk.VTK_PYRAMID
        elif nb_nodes_in_cell == 6:
            cell_type = vtk.VTK_WEDGE
        elif nb_nodes_in_cell == 8:
            cell_type = vtk.VTK_HEXAHEDRON       #there are 2 vtk cell_type with 8 nodes
        #elif nb_nodes_in_cell == 8:
            #cell_type = vtk.VTK_VOXEL
        elif nb_nodes_in_cell == 9:
            cell_type = vtk.VTK_QUADRATIC_TETRA
        elif nb_nodes_in_cell == 10:
            cell_type = vtk.VTK_PENTAGONAL_PRISM 
        elif nb_nodes_in_cell == 12:
            cell_type = vtk.VTK_HEXAGONAL_PRISM
        elif nb_nodes_in_cell == 20:
           cell_type = vtk.VTK_QUADRATIC_HEXAHEDRON
        # Add more cases for other cell types based on your logic
        else:
            cell_type = vtk.VTK_POLYHEDRON

        f.write('{0}\n'.format(cell_type))
    #-----------------------------------------------------------------------------------------------------------
    print("File generated:", output_filename)
    print("...done")
