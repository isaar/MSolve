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
        private readonly double[] coordinatesX;
        private readonly double[] coordinatesY;

        public int NodeRows { get; }
        public int NodeColumns { get; }

        public RectilinearMeshGenerator(double[] coordinatesX, double[] coordinatesY)
        {
            this.coordinatesX = coordinatesX;
            this.coordinatesY = coordinatesY;
            NodeRows = coordinatesY.Length;
            NodeColumns = coordinatesX.Length;
        }

        public RectilinearMeshGenerator(double minX, double maxX, double minY, double maxY, 
            double[] normalizedCoordinatesX, double[] normalizedCoordinatesY)
        {
            NodeRows = normalizedCoordinatesY.Length;
            NodeColumns = normalizedCoordinatesX.Length;

            coordinatesX = new double[NodeRows];
            coordinatesY = new double[NodeColumns];

            for (int row = 0; row < NodeRows; ++row)
            {
                coordinatesX[row] = minX + (maxX - minX) * normalizedCoordinatesX[row];
            }
            for (int col = 0; col < NodeColumns; ++col)
            {
                coordinatesY[col] = minY + (maxY - minY) * normalizedCoordinatesY[col];
            }
            
        }

        public Tuple<XNode2D[], List<XNode2D[]>> CreateMesh()
        {
            XNode2D[] nodes = CreateNodes();
            List<XNode2D[]> elementNodes = CreateElements(nodes);
            return new Tuple<XNode2D[], List<XNode2D[]>>(nodes, elementNodes);
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
    }
}
