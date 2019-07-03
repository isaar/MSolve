using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving
{
    class RelativeAreaResolver : IHeavisideSingularityResolver
    {
        private readonly double relativeAreaTolerance;
        private readonly Triangulator2D<CartesianPoint> triangulator;
        private readonly double zeroDistanceTolerance;

        public RelativeAreaResolver(double relativeAreaTolerance = 1E-4, double zeroDistanceTolerance = 1E-10)
        {
            this.relativeAreaTolerance = relativeAreaTolerance;
            this.zeroDistanceTolerance = zeroDistanceTolerance;
            this.triangulator = new Triangulator2D<CartesianPoint>((x, y) => new CartesianPoint(x, y));
        }

        /// <summary>
        /// Given a set of Heaviside enriched nodes, find which of them must not be enriched, in order to avoid the global
        /// stiffness matrix being singular.
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="heavisideNodes">They will not be altered.</param>
        /// <returns></returns>
        public ISet<XNode> FindHeavisideNodesToRemove(ISingleCrack crack, IMesh2D<XNode, XContinuumElement2D> mesh, 
            ISet<XNode> heavisideNodes)
        {
            var processedElements = new Dictionary<XContinuumElement2D, Tuple<double, double>>();
            var nodesToRemove = new HashSet<XNode>();
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

        /// <summary>
        /// Given a set of Heaviside enriched nodes, find which of them must not be enriched, in order to avoid the global
        /// stiffness matrix being singular.
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="heavisideNodes">They will not be altered.</param>
        /// <returns></returns>
        public ISet<XNode> FindHeavisideNodesToRemove(ISingleCrack crack, IReadOnlyList<XNode> heavisideNodes, 
            IReadOnlyList<ISet<XContinuumElement2D>> nodalSupports)
        {
            var processedElements = new Dictionary<XContinuumElement2D, Tuple<double, double>>();
            var nodesToRemove = new HashSet<XNode>();
            for (int i = 0; i < heavisideNodes.Count; ++i)
            {
                XNode node = heavisideNodes[i];
                double nodePositiveArea = 0.0;
                double nodeNegativeArea = 0.0;

                foreach (var element in nodalSupports[i])
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

        // TODO: I should really cache these somehow, so that they can be accessible from the crack object. They are used at various points.
        private (double positiveArea, double negativeArea) FindSignedAreasOfElement(ISingleCrack crack, 
            XContinuumElement2D element)
        {
            SortedSet<CartesianPoint> triangleVertices = crack.FindTriangleVertices(element);
            IReadOnlyList<Triangle2D<CartesianPoint>> triangles = triangulator.CreateMesh(triangleVertices);

            double positiveArea = 0.0;
            double negativeArea = 0.0;
            foreach (var triangle in triangles)
            {
                CartesianPoint v0 = triangle.Vertices[0];
                CartesianPoint v1 = triangle.Vertices[1];
                CartesianPoint v2 = triangle.Vertices[2];
                double area = 0.5 * Math.Abs(v0.X * (v1.Y - v2.Y) + v1.X * (v2.Y - v0.Y) + v2.X * (v0.Y - v1.Y));

                // The sign of the area can be derived from any node with body level set != 0
                int sign = 0;
                foreach (var vertex in triangle.Vertices)
                {
                    XNode vertexAsNode = null;
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
                        double distance = crack.SignedDistanceOf(vertexAsNode);
                        if (Math.Abs(distance) <= zeroDistanceTolerance) sign = 0;
                        else sign = Math.Sign(distance);
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


                    var centroid = new CartesianPoint((v0.X + v1.X + v2.X) / 3.0, (v0.Y + v1.Y + v2.Y) / 3.0);
                    NaturalPoint centroidNatural = element.Interpolation.
                        CreateInverseMappingFor(element.Nodes).TransformPointCartesianToNatural(centroid);
                    EvalInterpolation2D centroidInterpolation =
                        element.Interpolation.EvaluateAllAt(element.Nodes, centroidNatural);
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
