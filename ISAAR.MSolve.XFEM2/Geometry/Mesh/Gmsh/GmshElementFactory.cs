using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh
{
    class GmshElementFactory
    {
        private static readonly IReadOnlyDictionary<int, IsoparametricElementType2D> gmshElementCodes;

        // Node order for elements. Index = gmsh order, value = MSolve order.
        private static readonly IReadOnlyDictionary<IsoparametricElementType2D, int[]> gmshElementConnectivity;

        static GmshElementFactory()
        {
            var codes = new Dictionary<int, IsoparametricElementType2D>();
            codes.Add(2, IsoparametricElementType2D.Tri3);
            codes.Add(3, IsoparametricElementType2D.Quad4);
            codes.Add(10, IsoparametricElementType2D.Quad9);
            gmshElementCodes = codes;

            var connectivity = new Dictionary<IsoparametricElementType2D, int[]>();
            connectivity.Add(IsoparametricElementType2D.Tri3, new int[] { 0, 1, 2 });
            connectivity.Add(IsoparametricElementType2D.Quad4, new int[] { 0, 1, 2, 3 });
            connectivity.Add(IsoparametricElementType2D.Quad9, new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 });
            gmshElementConnectivity = connectivity;
        }

        private readonly IReadOnlyList<XNode2D> allNodes;

        public GmshElementFactory(IReadOnlyList<XNode2D> allNodes)
        {
            this.allNodes = allNodes;
        }

        /// <summary>
        /// Returns true and a GmshElement if the elementCode corresponds to a valid MSolve element type. 
        /// Otherwise returns false and null.
        /// </summary>
        /// <param name="elementCode"></param>
        /// <param name="nodeIDs"> These must be 0-based</param>
        /// <param name="element"></param>
        /// <returns></returns>
        public bool TryCreateElement(int elementCode, int[] nodeIDs, out GmshElement element)
        {
            IsoparametricElementType2D type;
            bool validElement = gmshElementCodes.TryGetValue(elementCode, out type);
            if (validElement)
            {
                XNode2D[] elementNodes = new XNode2D[nodeIDs.Length];
                for (int i = 0; i < nodeIDs.Length; ++i)
                {
                    int msolveIndex = gmshElementConnectivity[type][i];
                    elementNodes[msolveIndex] = allNodes[nodeIDs[i]];
                }
                ReverseNodeOrder(elementNodes);
                element = new GmshElement(type, elementNodes);
                return true;
            }
            else
            {
                element = null;
                return false;
            }
        }

        private void ReverseNodeOrder(XNode2D[] elementNodes)
        {
            // The area of the element with clockwise nodes is negative!
            double elementArea = 0.0; // Actually double the area, but we only care about the sign here
            for (int i = 0; i < elementNodes.Length; ++i)
            {
                XNode2D node1 = elementNodes[i];
                XNode2D node2 = elementNodes[(i + 1) % elementNodes.Length];
                elementArea += node1.X * node2.Y - node2.X * node1.Y;
            }

            if (elementArea < 0) Array.Reverse(elementNodes);
            return;
        }
    }
}
