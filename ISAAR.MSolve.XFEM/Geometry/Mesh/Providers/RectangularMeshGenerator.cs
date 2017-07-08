using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh.Providers
{
    class RectangularMeshGenerator
    {
        private readonly double dx;
        private readonly double dy;
        private readonly int nodeRows, nodeColumns, elementRows, elementColumns;

        public RectangularMeshGenerator(double domainDim1, double domainDim2, int elementsPerDim1, int elementsPerDim2)
        {
            this.dx = domainDim1 / elementsPerDim1;
            this.dy = domainDim2 / elementsPerDim2;
            this.elementRows = elementsPerDim1;
            this.elementColumns = elementsPerDim2;
            this.nodeRows = elementRows + 1;
            this.nodeColumns = elementColumns + 1;
        }

        public Tuple<XNode2D[], MeshFacePlaceholder2D[]> CreateMesh()
        {
            XNode2D[] nodes = CreateNodes();
            MeshFacePlaceholder2D[] elements = CreateElements(nodes);
            return new Tuple<XNode2D[], MeshFacePlaceholder2D[]>(nodes, elements);
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

        private MeshFacePlaceholder2D[] CreateElements(XNode2D[] nodes)
        {
            var elements = new MeshFacePlaceholder2D[elementRows * elementColumns];
            int id = 0;
            for (int row = 0; row < elementRows; ++row)
            {
                for (int col = 0; col < elementColumns; ++col)
                {
                    int firstNode = row * nodeColumns + col;
                    XNode2D[] elementNodes = { nodes[firstNode], nodes[firstNode+1],
                        nodes[firstNode + nodeColumns + 1], nodes[firstNode + nodeColumns] };
                    elements[id] = new MeshFacePlaceholder2D(elementNodes);
                    ++id;
                }
            }
            return elements;
        }
    }
}
