using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;

//TODO: Not sure if this works for 3D cells too
namespace ISAAR.MSolve.Preprocessor.Meshes.GMSH
{
    /// <summary>
    /// Converts cell types and the order of their vertices from GMSH to MSolve.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class GmshCellFactory_v2
    {
        private static readonly IReadOnlyDictionary<int, CellType> gmshCellCodes;

        // Vertex order for cells. Index = gmsh order, value = MSolve order.
        private static readonly IReadOnlyDictionary<CellType, int[]> gmshCellConnectivity;

        static GmshCellFactory_v2()
        {
            var codes = new Dictionary<int, CellType>();
            codes.Add(2, CellType.Tri3);
            codes.Add(3, CellType.Quad4);
            codes.Add(9, CellType.Tri6);
            codes.Add(10, CellType.Quad9);
            codes.Add(16, CellType.Quad8);
            gmshCellCodes = codes;

            var connectivity = new Dictionary<CellType, int[]>();
            connectivity.Add(CellType.Tri3, new int[] { 0, 1, 2 });
            connectivity.Add(CellType.Quad4, new int[] { 0, 1, 2, 3 });
            connectivity.Add(CellType.Tri6, new int[] { 0, 1, 2, 3, 4, 5 });
            connectivity.Add(CellType.Quad9, new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 });
            connectivity.Add(CellType.Quad8, new int[] { 0, 1, 2, 3, 4, 5, 6, 7 });
            gmshCellConnectivity = connectivity;
        }

        private readonly IReadOnlyList<Node_v2> allVertices;

        public GmshCellFactory_v2(IReadOnlyList<Node_v2> allVertices)
        {
            this.allVertices = allVertices;
        }

        /// <summary>
        /// Returns true and a <see cref="CellConnectivity_v2"/> if the <paramref name="cellCode"/> corresponds to a valid 
        /// MSolve <see cref="CellType"/>. 
        /// Otherwise returns false and null.
        /// </summary>
        /// <param name="cellCode"></param>
        /// <param name="vertexIDs"> These must be 0-based</param>
        /// <param name="cell"></param>
        /// <returns></returns>
        public bool TryCreateCell(int cellCode, int[] vertexIDs, out CellConnectivity_v2 cell)
        {
            bool validCell = gmshCellCodes.TryGetValue(cellCode, out CellType type);
            if (validCell)
            {
                var cellVertices = new Node_v2[vertexIDs.Length];
                for (int i = 0; i < vertexIDs.Length; ++i)
                {
                    int msolveIndex = gmshCellConnectivity[type][i];
                    cellVertices[msolveIndex] = allVertices[vertexIDs[i]];
                }
                FixVerticesOrder(cellVertices);
                cell = new CellConnectivity_v2(type, cellVertices);
                return true;
            }
            else
            {
                cell = null;
                return false;
            }
        }

        /// <summary>
        /// If the order is clockwise, it is reversed. Not sure if it sufficient or required for second order elements.
        /// </summary>
        /// <param name="cellVertices"></param>
        private void FixVerticesOrder(Node_v2[] cellVertices)
        {
            //TODO: This only works for 2D cells
            // The area of the cell with clockwise vertices is negative!
            double cellArea = 0.0; // Actually double the area will be computed, but we only care about the sign here
            for (int i = 0; i < cellVertices.Length; ++i)
            {
                Node_v2 vertex1 = cellVertices[i];
                Node_v2 vertex2 = cellVertices[(i + 1) % cellVertices.Length];
                cellArea += vertex1.X * vertex2.Y - vertex2.X * vertex1.Y;
            }
            if (cellArea < 0) Array.Reverse(cellVertices);
            return;
        }
    }
}
