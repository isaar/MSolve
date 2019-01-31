using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;

//TODO: abstract this in order to be used with points in various coordinate systems
//TODO: perhaps the origin should be (0.0, 0.0) and the meshes could then be transformed. Abaqus does something similar with its
//      meshed parts during assembly
namespace ISAAR.MSolve.Preprocessor.Meshes.Custom
{
    /// <summary>
    /// Creates 3D meshes based on uniform rectilinear grids: the distance between two consecutive vertices for the same axis is 
    /// constant. This distance may be different for each axis though. For now the cells are hexahedral with 8 vertices. 
    /// (bricks in particular).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class UniformMeshGenerator3D : IMeshProvider2D<Node_v2, CellConnectivity_v2>
    {
        private readonly double minX, minY, minZ;
        private readonly double dx, dy, dz;
        private readonly int cellsPerX, cellsPerY, cellsPerZ;
        private readonly int verticesPerX, verticesPerY, verticesPerZ;

        public UniformMeshGenerator3D(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, 
            int cellsPerX, int cellsPerY, int cellsPerZ)
        {
            this.minX = minX;
            this.minY = minY;
            this.minZ = minZ;
            this.dx = (maxX - minX) / cellsPerX;
            this.dy = (maxY - minY) / cellsPerY;
            this.dz = (maxZ - minZ) / cellsPerZ;
            this.cellsPerX = cellsPerX;
            this.cellsPerY = cellsPerY;
            this.cellsPerZ = cellsPerZ;
            this.verticesPerX = cellsPerX + 1;
            this.verticesPerY = cellsPerY + 1;
            this.verticesPerZ = cellsPerZ + 1;
        }

        public bool StartIDsAt0 { get; set; } = true;

        /// <summary>
        /// Generates a uniform mesh with the dimensions and density defined in the constructor.
        /// </summary>
        /// <returns></returns>
        public (IReadOnlyList<Node_v2> vertices, IReadOnlyList<CellConnectivity_v2> cells) CreateMesh()
        {
            Node_v2[] vertices = CreateVertices();
            CellConnectivity_v2[] cells = CreateCells(vertices);
            return (vertices, cells);
        }

        private Node_v2[] CreateVertices()
        {
            var vertices = new Node_v2[verticesPerX * verticesPerY * verticesPerZ];
            int id = 0;
            int start = StartIDsAt0 ? 0 : 1;
            for (int k = 0; k < verticesPerZ; ++k)
            {
                for (int j = 0; j < verticesPerY; ++j)
                {
                    for (int i = 0; i < verticesPerX; ++i)
                    {
                        vertices[id] = new Node_v2 { ID = start + id, X = minX + i * dx, Y = minY + j * dy, Z = minZ + k * dz };
                        ++id;
                    }
                }
            }
            return vertices;
        }

        private CellConnectivity_v2[] CreateCells(Node_v2[] allVertices)
        {
            var cells = new CellConnectivity_v2[cellsPerX * cellsPerY * cellsPerZ];
            for (int k = 0; k < cellsPerZ; ++k)
            {
                for (int j = 0; j < cellsPerY; ++j)
                {
                    for (int i = 0; i < cellsPerX; ++i)
                    {
                        int cell = k * cellsPerX * cellsPerY + j * cellsPerX + i;
                        int firstVertex = k * verticesPerX * verticesPerY + j * verticesPerX + i;
                        var verticesOfCell = new int[]                                 // Node order for Hexa8
                        {
                            firstVertex,                                                    // (-1, -1, -1)
                            firstVertex + 1,                                                // ( 1, -1, -1)
                            firstVertex + verticesPerY + 1,                                 // ( 1,  1, -1)
                            firstVertex + verticesPerY,                                     // (-1,  1, -1)
                            firstVertex + verticesPerX * verticesPerY,                      // (-1, -1,  1)
                            firstVertex + verticesPerX * verticesPerY + 1,                  // ( 1, -1,  1)
                            firstVertex + verticesPerX * verticesPerY + verticesPerY + 1,   // ( 1,  1,  1)
                            firstVertex + verticesPerX * verticesPerY + verticesPerY        // (-1,  1,  1)
                        };
                        cells[cell] = new CellConnectivity_v2(CellType.Hexa8, 
                            verticesOfCell.Select(idx => allVertices[idx]).ToArray()); // row major
                    }
                }
            }
            return cells;
        }
    }
}
