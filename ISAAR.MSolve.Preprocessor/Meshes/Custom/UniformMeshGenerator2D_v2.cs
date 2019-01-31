using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;

//TODO: abstract this in order to be used with points in various coordinate systems
//TODO: perhaps the origin should be (0.0, 0.0) and the meshes could then be transformed. Abaqus does something similar with its
//      meshed parts during assembly
namespace ISAAR.MSolve.Preprocessor.Meshes.Custom
{
    /// <summary>
    /// Creates 2D meshes based on uniform rectilinear grids: the distance between two consecutive vertices for the same axis is 
    /// constant. This distance may be different for each axis though. For now the cells are quadrilateral with 4 vertices 
    /// (rectangles in particular).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class UniformMeshGenerator2D_v2 : IMeshProvider2D<Node_v2, CellConnectivity_v2>
    {
        private readonly double minX, minY;
        private readonly double dx, dy;
        private readonly int cellsPerX, cellsPerY;
        private readonly int verticesPerX, verticesPerY;

        public UniformMeshGenerator2D_v2(double minX, double minY, double maxX, double maxY, int cellsPerX, int cellsPerY)
        {
            this.minX = minX;
            this.minY = minY;
            this.dx = (maxX - minX) / cellsPerX;
            this.dy = (maxY - minY) / cellsPerY;
            this.cellsPerX = cellsPerX;
            this.cellsPerY = cellsPerY;
            this.verticesPerX = this.cellsPerX + 1;
            this.verticesPerY = this.cellsPerY + 1;
        }

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
            var vertices = new Node_v2[verticesPerY * verticesPerX];
            int id = 0;
            for (int j = 0; j < verticesPerY; ++j)
            {
                for (int i = 0; i < verticesPerX; ++i)
                {
                    vertices[id] = new Node_v2 { ID = id, X = minX + i * dx , Y = minY + j * dy };
                    ++id;
                }
            }
            return vertices;
        }

        private CellConnectivity_v2[] CreateCells(Node_v2[] allVertices)
        {
            var cells = new CellConnectivity_v2[cellsPerY * cellsPerX];
            for (int j = 0; j < cellsPerY; ++j)
            {
                for (int i = 0; i < cellsPerX; ++i)
                {
                    int cell = j * cellsPerX + i;
                    int firstVertex = j * verticesPerX + i;
                    Node_v2[] verticesOfCell = 
                    {
                        allVertices[firstVertex], allVertices[firstVertex+1],
                        allVertices[firstVertex + verticesPerX + 1], allVertices[firstVertex + verticesPerX]
                    };
                    cells[cell] = new CellConnectivity_v2(CellType.Quad4, verticesOfCell); // row major
                }
            }
            return cells;
        }
    }
}
