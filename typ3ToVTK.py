"""
Karnon K.
25/06/2023
Copyright (c) 2023,Karnon K.
"""
import sys
import vtk
import numpy as np
import argparse

def verify(cellVertices, verticesCoords):
    nbVertices = len(cellVertices)

    meanVertex = np.array([0. for d in range(3)])

    for vertex in cellVertices:
        coord = np.array(verticesCoords[vertex - 1])
        meanVertex += coord

    meanVertex *= 1 / nbVertices
    print(meanVertex)

# options
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="file name", default="")
parser.add_argument("-p", "--polyhedron", help="enforce polyhedron cell format",
                    action="store_true", default=False)
args = parser.parse_args()

filename = args.file
outFileName = filename.replace(".msh", ".vtk")
print("Converting {0} into {1}".format(filename, outFileName))

verticesPart = False
cellsPart = False
facesPart = False
cellConnectedItemsReadingStarted = False
faceVerticesReadingStarted = False
verticesCoords = []
cellsConnectedItems = []
cellConnectedItems = []
nbCellConnectedItems = -1
facesVertices = []
faceVertices = []
nbFaceVertices = -1

with open(filename, 'r', encoding='ISO-8859-1') as f:
    for line in f:
        infos = line.split()

        if not infos:
            continue

        if infos[0] == "Vertices":
            verticesPart = True
            continue

        if infos[0].startswith("Volumes"):
            verticesPart = False
            if args.polyhedron:
                if infos[0].startswith("Volumes->faces"):
                    cellsPart = True
                if infos[0].startswith("Volumes->Vertices"):
                    cellsPart = False
            else:
                if infos[0].startswith("Volumes->Vertices"):
                    cellsPart = True
            continue

        if infos[0].startswith("Faces"):
            cellsPart = False
            if infos[0].startswith("Faces->Vertices"):
                facesPart = True
            if infos[0].startswith("Faces->Control"):
                facesPart = False
            continue

        if verticesPart:
            vertexCoord = [float(v) for v in infos]
            verticesCoords.append(vertexCoord)

        if cellsPart:
            values = [int(v) for v in infos]
            if cellConnectedItemsReadingStarted:
                cellConnectedItems += values
            else:
                cellConnectedItemsReadingStarted = True
                nbCellConnectedItems = int(infos[0])
                cellConnectedItems = values[1:]

            if len(cellConnectedItems) == nbCellConnectedItems:
                cellsConnectedItems.append(cellConnectedItems)
                cellConnectedItemsReadingStarted = False
                cellConnectedItems = []
                nbCellConnectedItems = -1

        if facesPart:
            values = [int(v) for v in infos]
            if faceVerticesReadingStarted:
                faceVertices += values
            else:
                faceVerticesReadingStarted = True
                nbFaceVertices = int(infos[0])
                faceVertices = values[1:]

            if len(faceVertices) == nbFaceVertices:
                facesVertices.append(faceVertices)
                faceVerticesReadingStarted = False
                faceVertices = []
                nbFaceVertices = -1

with open(outFileName, 'w') as f:
    nbVertices = len(verticesCoords)
    f.write('# vtk DataFile Version 2.0\n')
    f.write('Really cool data\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS {0} float\n'.format(nbVertices))
    for v in range(len(verticesCoords)):
        f.write('{0} {1} {2}\n'.format(verticesCoords[v][0], verticesCoords[v][1], verticesCoords[v][2]))

    #---- Write the cells----------------------------------------------------------------------------
    nbCells = len(cellsConnectedItems)
    if args.polyhedron:
        print("polyhedron chosen")
        totalCellIndices = 0
        cellSumIndices = []
        for cell in cellsConnectedItems:
            totalCellIndices += 2
            n = 1
            for face in cell:
                tmp = len(facesVertices[face-1]) + 1
                totalCellIndices += tmp
                n += tmp
            cellSumIndices.append(n)
        f.write('CELLS {0} {1}\n'.format(nbCells, totalCellIndices))
        for c in range(len(cellsConnectedItems)):
            cell = cellsConnectedItems[c]
            f.write('{0} {1}\n'.format(cellSumIndices[c], len(cell)))
            for face in cell:
                f.write('{0} '.format(len(facesVertices[face-1])))
                f.write(' '.join(str(node - 1) for node in facesVertices[face-1]))
                f.write('\n')
    else:
        totalCellIndices = sum(len(cell) + 1 for cell in cellsConnectedItems)
        f.write('CELLS {0} {1}\n'.format(nbCells, totalCellIndices))
        for cell in cellsConnectedItems:
            nbCellConnectedItems = len(cell)
            if(nbCellConnectedItems == 8):
                f.write('{0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(nbCellConnectedItems,
                                                                     cell[0] - 1,
                                                                     cell[4] - 1,
                                                                     cell[6] - 1,
                                                                     cell[2] - 1,
                                                                     cell[1] - 1,
                                                                     cell[5] - 1,
                                                                     cell[7] - 1,
                                                                     cell[3] - 1))
            else:
                f.write('{0} '.format(nbCellConnectedItems))
                f.write(' '.join(str(node - 1) for node in cell))  # Correct indexing to start from 0
            f.write('\n')


    # ---Write the cell types-------------------------------------------------------------------------
    f.write('CELL_TYPES {0}\n'.format(nbCells))
    if args.polyhedron:
        for nb_nodes_in_cell in [len(cell) for cell in cellsConnectedItems]:
            cell_type = vtk.VTK_POLYHEDRON
            f.write('{0}\n'.format(cell_type))
    else:
        for nb_nodes_in_cell in [len(cell) for cell in cellsConnectedItems]:
            if nb_nodes_in_cell == 4:
                cell_type = vtk.VTK_TETRA
            elif nb_nodes_in_cell == 5:
                cell_type = vtk.VTK_PYRAMID
            elif nb_nodes_in_cell == 6:
                cell_type = vtk.VTK_WEDGE
            elif nb_nodes_in_cell == 8:
                cell_type = vtk.VTK_HEXAHEDRON
            elif nb_nodes_in_cell == 9:
                cell_type = vtk.VTK_QUADRATIC_TETRA
            elif nb_nodes_in_cell == 10:
                cell_type = vtk.VTK_PENTAGONAL_PRISM
            elif nb_nodes_in_cell == 12:
                cell_type = vtk.VTK_HEXAGONAL_PRISM
            elif nb_nodes_in_cell == 20:
                cell_type = vtk.VTK_QUADRATIC_HEXAHEDRON
            f.write('{0}\n'.format(cell_type))
    #-----------------------------------------------------------------------------------------------------------
    print("File generated:", outFileName)

    #verify(cellsVertices[0], verticesCoords)

    print("...done")
