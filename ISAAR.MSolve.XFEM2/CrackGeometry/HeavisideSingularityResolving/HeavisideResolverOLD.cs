using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Interpolation;

//TODO: I suspect this is method is responsible for singular global matrices popping up randomly for many configurations.
//      Either this doesn't work as good as advertised or I messed it up. Try checking if there is at least  
//      least 1 GP on either side of the crack instead. Bonus points if the GPs are cached somehow (e.g. all std elements
//      have the same GPs in natural system, many enriched elements may also have the same active integration rule) for
//      when the integration actually happens.
namespace ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving
{
    class HeavisideResolverOLD: IHeavisideSingularityResolver
    {
        private readonly double relativeAreaTolerance;
        private readonly CartesianTriangulator triangulator;

        public HeavisideResolverOLD(double relativeAreaTolerance = 1e-4)
        {
            this.relativeAreaTolerance = relativeAreaTolerance;
            this.triangulator = new CartesianTriangulator();
        }

        /// <summary>
        /// Given a set of Heaviside enriched nodes, find which of them must not be enriched, in order to avoid the global
        /// stiffness matrix being singular.
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="heavisideNodes">They will not be altered.</param>
        /// <returns></returns>
        public ISet<XNode2D> FindHeavisideNodesToRemove(ISingleCrack crack, IMesh2D<XNode2D, XContinuumElement2D> mesh, 
            ISet<XNode2D> heavisideNodes)
        {
            var processedElements = new Dictionary<XContinuumElement2D, Tuple<double, double>>();
            var nodesToRemove = new HashSet<XNode2D>();
            foreach (var node in heavisideNodes)
            {
                double nodePositiveArea = 0.0;
                double nodeNegativeArea = 0.0;

                foreach (var element in mesh.FindElementsWithNode(node))
                {
                    bool alreadyProcessed = processedElements.TryGetValue(element, out Tuple<double, double> elementPosNegAreas);
                    if (!alreadyProcessed)
                    {
                        (double elementPosArea, double elementNegArea) = FindSignedAreasOfElement(crack, element);
                        elementPosNegAreas = new Tuple<double, double>(elementPosArea, elementNegArea);
                        processedElements[element] = elementPosNegAreas;
                    }
                    nodePositiveArea += elementPosNegAreas.Item1;
                    nodeNegativeArea += elementPosNegAreas.Item2;
                }

                if (crack.SignedDistanceOf(node) >= 0.0)
                {
                    double negativeAreaRatio = nodeNegativeArea / (nodePositiveArea + nodeNegativeArea);
                    if (negativeAreaRatio < relativeAreaTolerance) nodesToRemove.Add(node);
                }
                else
                {
                    double positiveAreaRatio = nodePositiveArea / (nodePositiveArea + nodeNegativeArea);
                    if (positiveAreaRatio < relativeAreaTolerance) nodesToRemove.Add(node);
                }
            }

            return nodesToRemove;
        }

        public ISet<XNode2D> FindHeavisideNodesToRemove(ISingleCrack crack, IReadOnlyList<XNode2D> heavisideNodes, IReadOnlyList<ISet<XContinuumElement2D>> nodalSupports)
        {
            throw new NotImplementedException();
        }

        private (double positiveArea, double negativeArea) FindSignedAreasOfElement(ISingleCrack crack, 
            XContinuumElement2D element)
        {
            #region Debug
            //if (element.Nodes[0].ID==154)
            //{
            //    Console.WriteLine();
            //}
            #endregion

            SortedSet<ICartesianPoint2D> triangleVertices = crack.FindTriangleVertices(element);
            IReadOnlyList<TriangleCartesian2D> triangles = triangulator.CreateMesh(triangleVertices);

            double positiveArea = 0.0;
            double negativeArea = 0.0;
            foreach (var triangle in triangles)
            {
                ICartesianPoint2D v0 = triangle.Vertices[0];
                ICartesianPoint2D v1 = triangle.Vertices[1];
                ICartesianPoint2D v2 = triangle.Vertices[2];
                double area = 0.5 * Math.Abs(v0.X * (v1.Y - v2.Y) + v1.X * (v2.Y - v0.Y) + v2.X * (v0.Y - v1.Y));

                // The sign of the area can be derived from any node with body level set != 0
                int sign = 0;
                foreach (var vertex in triangle.Vertices)
                {
                    XNode2D vertexAsNode = null;
                    foreach (var node in element.Nodes) // TODO: find a faster way to do this
                    {
                        if ((vertex.X == node.X) && (vertex.Y == node.Y))
                        {
                            vertexAsNode = node;
                            break;
                        }
                    }
                    if (vertexAsNode != null)
                    {
                        sign = Math.Sign(crack.SignedDistanceOf(vertexAsNode));
                        if (sign != 0) break;
                    }
                }

                // If no node with non-zero body level set is found, then find the body level set of its centroid
                if (sign == 0)
                {
                    // Report this instance in DEBUG messages. It should not happen with linear level sets and only 1 crack.
                    //if (reports)
                    //{
                    //    Console.WriteLine("--- DEBUG: Triangulation resulted in a triangle where no vertex is an element node. ---");
                    //}


                    var centroid = new CartesianPoint2D((v0.X + v1.X + v2.X) / 3.0, (v0.Y + v1.Y + v2.Y) / 3.0);
                    INaturalPoint2D centroidNatural = element.Interpolation.
                        CreateInverseMappingFor(element.Nodes).TransformCartesianToNatural(centroid);
                    EvaluatedInterpolation2D centroidInterpolation =
                        element.Interpolation.EvaluateAt(element.Nodes, centroidNatural);
                    sign = Math.Sign(crack.SignedDistanceOf(centroidNatural, element, centroidInterpolation));
                }

                if (sign > 0) positiveArea += area;
                else if (sign < 0) negativeArea += area;
                else throw new Exception(
                    "Even after finding the signed distance of its centroid, the sign of the area is unidentified");
            }

            return (positiveArea, negativeArea);
        }
    }
}
