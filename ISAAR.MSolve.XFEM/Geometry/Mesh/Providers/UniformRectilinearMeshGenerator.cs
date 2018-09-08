using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh.Providers
{
    class UniformRectilinearMeshGenerator
    {
        private readonly double dx;
        private readonly double dy;
        private readonly int nodeRows, nodeColumns, elementRows, elementColumns;

        public UniformRectilinearMeshGenerator(double domainDim1, double domainDim2, int elementsPerDim1, int elementsPerDim2)
        {
            this.dx = domainDim1 / elementsPerDim1;
            this.dy = domainDim2 / elementsPerDim2;
            this.elementRows = elementsPerDim2;
            this.elementColumns = elementsPerDim1;
            this.nodeRows = elementRows + 1;
            this.nodeColumns = elementColumns + 1;
        }

        public (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateMesh()
        {
            XNode2D[] nodes = CreateNodes();
            List<XNode2D[]> elementConnectivity = CreateElements(nodes);
            return (nodes, elementConnectivity);
        }

        private XNode2D[] CreateNodes()
        {
            var nodes = new XNode2D[nodeRows * nodeColumns];
            int id = 0;
            for (int row = 0; row < nodeRows; ++row)
            {
                for (int col = 0; col < nodeColumns; ++col)
                {
                    nodes[id] = new XNode2D(id, col * dx, row * dy);
                    ++id;
                }
            }
            return nodes;
        }

        private List<XNode2D[]> CreateElements(XNode2D[] allNodes)
        {
            var elementNodes = new List<XNode2D[]>();
            for (int row = 0; row < elementRows; ++row)
            {
                for (int col = 0; col < elementColumns; ++col)
                {
                    int firstNode = row * nodeColumns + col;
                    XNode2D[] nodesOfElement = { allNodes[firstNode], allNodes[firstNode+1],
                        allNodes[firstNode + nodeColumns + 1], allNodes[firstNode + nodeColumns] };
                    elementNodes.Add(nodesOfElement);
                }
            }
            return elementNodes;
        }
    }
}
