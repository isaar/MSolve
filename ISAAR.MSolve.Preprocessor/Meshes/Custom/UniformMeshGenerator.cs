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
    /// Creates meshes based on uniform rectilinear grids: the distance between two consecutive vertices for the same axis is 
    /// constant. This distance may be different for each axis though. For now the cells are quadrilateral with 4 vertices 
    /// (rectangles in particular).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class UniformMeshGenerator : IMeshProvider2D<Node2D, CellConnectivity2D>
    {
        private readonly double minX;
        private readonly double minY;
        private readonly double dx;
        private readonly double dy;
        private readonly int vertexRows, vertexColumns, cellRows, cellColumns;

        public UniformMeshGenerator(double minX, double minY, double maxX, double maxY, int cellsPerX, int cellsPerY)
        {
            this.minX = minX;
            this.minY = minY;
            this.dx = (maxX - minX) / cellsPerX;
            this.dy = (maxY - minY) / cellsPerY;
            this.cellRows = cellsPerY;
            this.cellColumns = cellsPerX;
            this.vertexRows = cellRows + 1;
            this.vertexColumns = cellColumns + 1;
        }

        /// <summary>
        /// Generates a uniform mesh with the dimensions and density defined in the constructor.
        /// </summary>
        /// <returns></returns>
        public (IReadOnlyList<Node2D> vertices, IReadOnlyList<CellConnectivity2D> cells) CreateMesh()
        {
            Node2D[] vertices = CreateVertices();
            CellConnectivity2D[] cells = CreateCells(vertices);
            return (vertices, cells);
        }

        private Node2D[] CreateVertices()
        {
            var vertices = new Node2D[vertexRows * vertexColumns];
            int id = 0;
            for (int row = 0; row < vertexRows; ++row)
            {
                for (int col = 0; col < vertexColumns; ++col)
                {
                    vertices[id] = new Node2D(id, minX + col * dx, minY + row * dy);
                    ++id;
                }
            }
            return vertices;
        }

        private CellConnectivity2D[] CreateCells(Node2D[] allVertices)
        {
            var cells = new CellConnectivity2D[cellRows * cellColumns];
            for (int row = 0; row < cellRows; ++row)
            {
                for (int col = 0; col < cellColumns; ++col)
                {
                    int firstVertex = row * vertexColumns + col;
                    Node2D[] verticesOfCell = { allVertices[firstVertex], allVertices[firstVertex+1],
                        allVertices[firstVertex + vertexColumns + 1], allVertices[firstVertex + vertexColumns] };
                    cells[row * cellColumns + col] = new CellConnectivity2D(CellType.Quad4, verticesOfCell); // row major
                }
            }
            return cells;
        }
    }
}
