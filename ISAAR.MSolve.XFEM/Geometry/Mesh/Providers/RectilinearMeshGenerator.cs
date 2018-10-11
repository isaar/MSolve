using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh.Providers
{
    class RectilinearMeshGenerator
    {
        private readonly IReadOnlyList<double> coordinatesX;
        private readonly IReadOnlyList<double> coordinatesY;

        public int NodeRows { get; }
        public int NodeColumns { get; }

        public RectilinearMeshGenerator(double[,] meshSizeAlongX, double[,] meshSizeAlongY):
            this(PartitionLine(meshSizeAlongX), PartitionLine(meshSizeAlongY))
        {
        }

        public RectilinearMeshGenerator(IReadOnlyList<double> coordinatesX, IReadOnlyList<double> coordinatesY)
        {
            this.coordinatesX = coordinatesX;
            this.coordinatesY = coordinatesY;
            NodeRows = coordinatesY.Count;
            NodeColumns = coordinatesX.Count;
        }

        public (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateMesh()
        {
            XNode2D[] nodes = CreateNodes();
            List<XNode2D[]> elementConnectivity = CreateElements(nodes);
            return (nodes, elementConnectivity);
        }

        private XNode2D[] CreateNodes()
        {
            var nodes = new XNode2D[NodeColumns * NodeRows];
            int id = 0;
            for (int row = 0; row < NodeRows; ++row)
            {
                for (int col = 0; col < NodeColumns; ++col)
                {
                    nodes[id] = new XNode2D(id, coordinatesX[col], coordinatesY[row]);
                    ++id;
                }
            }
            return nodes;
        }

        private List<XNode2D[]> CreateElements(XNode2D[] allNodes)
        {
            int elementRows = NodeRows - 1;
            int elementColumns = NodeColumns - 1;
            var elementNodes = new List<XNode2D[]>();
            for (int row = 0; row < elementRows; ++row)
            {
                for (int col = 0; col < elementColumns; ++col)
                {
                    int firstNode = row * NodeColumns + col;
                    XNode2D[] nodesOfElement = { allNodes[firstNode], allNodes[firstNode+1],
                        allNodes[firstNode + NodeColumns + 1], allNodes[firstNode + NodeColumns] };
                    elementNodes.Add(nodesOfElement);
                }
            }
            return elementNodes;
        }

        private static List<double> PartitionLine(double[,] meshSizes)
        {
            int n = meshSizes.GetLength(0);
            var points = new List<double>();
            for (int i = 0; i < n - 1; ++i)
            {
                points.Add(meshSizes[i, 0]);
                points.AddRange(
                    PartitionSegment(meshSizes[i, 0], meshSizes[i, 1], meshSizes[i + 1, 0], meshSizes[i + 1, 1]));
            }
            points.Add(meshSizes[n - 1, 0]);
            return points;
        }

        /// <summary>
        /// x1 < x2
        /// </summary>
        private static List<double> PartitionSegment(double x1, double h1, double x2, double h2)
        {
            var points = new List<double>();
            double slope = (h2 - h1) / (x2 - x1);
            if (h1 <= h2)
            {
                double x = x1;
                double h = h1;
                while (true)
                {
                    x = x + h;
                    h = slope * (x - x1) + h1;
                    if (x < x2) points.Add(x);
                    else return points;
                }
            }
            else
            {
                double x = x2;
                double h = h2;
                //points.Add(x); Do not add the edges
                while (true)
                {
                    x = x - h;
                    h = slope * (x - x1) + h1;
                    if (x > x1) points.Add(x);
                    else
                    {
                        points.Reverse();
                        return points;
                    }
                }
            }
        }
    }
}
