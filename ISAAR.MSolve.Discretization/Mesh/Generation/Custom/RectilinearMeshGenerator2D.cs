using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Discretization.Mesh.Generation.Custom
{
    public class RectilinearMeshGenerator2D<TNode> : IMeshGenerator<TNode> where TNode : INode
    {
        private readonly IReadOnlyList<double> coordinatesX;
        private readonly IReadOnlyList<double> coordinatesY;

        public int NodeRows { get; }
        public int NodeColumns { get; }

        public RectilinearMeshGenerator2D(double[,] meshSizeAlongX, double[,] meshSizeAlongY) :
            this(PartitionLine(meshSizeAlongX), PartitionLine(meshSizeAlongY))
        {
        }

        public RectilinearMeshGenerator2D(IReadOnlyList<double> coordinatesX, IReadOnlyList<double> coordinatesY)
        {
            this.coordinatesX = coordinatesX;
            this.coordinatesY = coordinatesY;
            NodeRows = coordinatesY.Count;
            NodeColumns = coordinatesX.Count;
        }

        public (IReadOnlyList<TNode> nodes, IReadOnlyList<CellConnectivity<TNode>> elements) 
            CreateMesh(CreateNode<TNode> createNode)
        {
            TNode[] nodes = CreateNodes(createNode);
            IReadOnlyList<CellConnectivity<TNode>> elementConnectivity = CreateElements(nodes);
            return (nodes, elementConnectivity);
        }

        private TNode[] CreateNodes(CreateNode<TNode> createNode)
        {
            var nodes = new TNode[NodeColumns * NodeRows];
            int id = 0;
            for (int row = 0; row < NodeRows; ++row)
            {
                for (int col = 0; col < NodeColumns; ++col)
                {
                    nodes[id] = createNode(id, coordinatesX[col], coordinatesY[row], 0.0);
                    ++id;
                }
            }
            return nodes;
        }

        private IReadOnlyList<CellConnectivity<TNode>> CreateElements(TNode[] allNodes)
        {
            int elementRows = NodeRows - 1;
            int elementColumns = NodeColumns - 1;
            var elements = new List<CellConnectivity<TNode>>();
            for (int row = 0; row < elementRows; ++row)
            {
                for (int col = 0; col < elementColumns; ++col)
                {
                    int firstNode = row * NodeColumns + col;
                    TNode[] nodesOfElement = { allNodes[firstNode], allNodes[firstNode+1],
                        allNodes[firstNode + NodeColumns + 1], allNodes[firstNode + NodeColumns] };
                    elements.Add(new CellConnectivity<TNode>(CellType.Quad4, nodesOfElement));  // row major
                }
            }
            return elements;
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
