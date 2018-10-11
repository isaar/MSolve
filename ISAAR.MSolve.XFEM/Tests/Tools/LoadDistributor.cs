using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class LoadDistributor
    {

        public LoadDistributor()
        {
        }

        public double[,] DistributeLoad(IReadOnlyList<XNode2D> edgeNodes, double tractionX, double tractionY)
        {
            double[] nodalAreasOfInfluence = new double[edgeNodes.Count];
            int lastNode = edgeNodes.Count - 1;
            nodalAreasOfInfluence[0] = 0.5 * DistanceOf(edgeNodes[0], edgeNodes[1]);
            for (int i = 1; i < lastNode; ++i)
            {
                double areaLeft = 0.5 * DistanceOf(edgeNodes[i - 1], edgeNodes[i]);
                double areaRight = 0.5 * DistanceOf(edgeNodes[i], edgeNodes[i+1]);
                nodalAreasOfInfluence[i] = areaLeft + areaRight;
            }
            nodalAreasOfInfluence[lastNode] = 0.5 * DistanceOf(edgeNodes[lastNode - 1], edgeNodes[lastNode]);

            double[,] nodalLoads = new double[edgeNodes.Count, 2];
            for (int i = 0; i < edgeNodes.Count; ++i)
            {
                nodalLoads[i, 0] = nodalAreasOfInfluence[i] * tractionX;
                nodalLoads[i, 1] = nodalAreasOfInfluence[i] * tractionY;
            }

            return nodalLoads;
        }

        private static double DistanceOf(XNode2D node1, XNode2D node2)
        {
            double dx = node2.X - node1.X;
            double dy = node2.Y - node1.Y;
            return Math.Sqrt(dx * dx + dy * dy);
        }
    }
}
