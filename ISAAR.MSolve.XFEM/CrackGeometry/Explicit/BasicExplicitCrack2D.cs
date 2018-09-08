using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;


namespace ISAAR.MSolve.XFEM.CrackGeometry.Explicit
{
    class BasicExplicitCrack2D: ISingleCrack
    {
        private static readonly bool reports = false;
        private static readonly IComparer<ICartesianPoint2D> pointComparer = new Point2DComparerXMajor();

        private readonly double tipEnrichmentAreaRadius;
        private readonly CartesianTriangulator triangulator;

        public ICartesianPoint2D CrackMouth { get; private set; }

        // TODO: Not too fond of the setters, but at least the enrichments are immutable. Perhaps I can pass their
        // parameters to a CrackDescription builder and construct them there, without involving the user 
        // (given how easy it is to forget the setters, it is a must).
        public IMesh2D<XNode2D, XContinuumElement2D> Mesh { get; set; }
        public CrackBodyEnrichment2D CrackBodyEnrichment { get; set; }
        public CrackTipEnrichments2D CrackTipEnrichments { get; set; }

        private List<ICartesianPoint2D> Vertices { get; }
        private List<DirectedSegment2D> Segments { get; }
        // Angles[i-1] is the angle of segment i w.r.t segment i-1, aka the crack growth angle.
        private List<double> Angles { get; }

        public ISet<XContinuumElement2D> ElementsModified => throw new NotImplementedException();

        public IReadOnlyList<IEnrichmentItem2D> Enrichments => throw new NotImplementedException();

        BiMesh2D ICrackDescription.Mesh => throw new NotImplementedException();

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesAll => throw new NotImplementedException();

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesNew => throw new NotImplementedException();

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesModified => throw new NotImplementedException();

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesNearModified => throw new NotImplementedException();

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode2D>> CrackTipNodesNew => throw new NotImplementedException();

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode2D>> CrackTipNodesOld => throw new NotImplementedException();

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode2D>> CrackBodyNodesRejected => throw new NotImplementedException();

        public IHeavisideSingularityResolver SingularityResolver => throw new NotImplementedException();

        public IReadOnlyList<ISingleCrack> SingleCracks => throw new NotImplementedException();

        public IReadOnlyDictionary<ICartesianPoint2D, IReadOnlyList<XContinuumElement2D>> CrackTipElements => throw new NotImplementedException();

        public IReadOnlyList<ICartesianPoint2D> CrackTips => throw new NotImplementedException();

        public IReadOnlyDictionary<ICartesianPoint2D, IPropagator> CrackTipPropagators => throw new NotImplementedException();

        private TipCoordinateSystem tipSystem;
        private List<XContinuumElement2D> tipElements;

        public BasicExplicitCrack2D(double tipEnrichmentAreaRadius = 0.0)
        {
            this.tipEnrichmentAreaRadius = tipEnrichmentAreaRadius;
            this.triangulator = new CartesianTriangulator();
            this.tipElements = new List<XContinuumElement2D>();

            Vertices = new List<ICartesianPoint2D>();
            Segments = new List<DirectedSegment2D>();
            Angles = new List<double>();
        }

        public ICartesianPoint2D GetCrackTip(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Single) return Vertices[Vertices.Count - 1];
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        public TipCoordinateSystem GetTipSystem(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Single) return tipSystem;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        public IReadOnlyList<XContinuumElement2D> GetTipElements(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Single) return tipElements;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        public void InitializeGeometry(ICartesianPoint2D crackMouth, ICartesianPoint2D crackTip)
        {
            CrackMouth = crackMouth;

            double dx = crackTip.X - crackMouth.X;
            double dy = crackTip.Y - crackMouth.Y;
            double tangentSlope = Math.Atan2(dy, dx);
            tipSystem = new TipCoordinateSystem(crackTip, tangentSlope);
            CrackTipEnrichments.TipSystem = tipSystem;

            Vertices.Add(crackMouth);
            Vertices.Add(crackTip);
            Segments.Add(new DirectedSegment2D(crackMouth, crackTip));
        }

        public void UpdateGeometry(double localGrowthAngle, double growthLength)
        {
            double globalGrowthAngle = AngleUtilities.Wrap(localGrowthAngle + tipSystem.RotationAngle);
            double dx = growthLength * Math.Cos(globalGrowthAngle);
            double dy = growthLength * Math.Sin(globalGrowthAngle);

            var oldTip = Vertices[Vertices.Count - 1];
            var newTip = new CartesianPoint2D(oldTip.X + dx, oldTip.Y + dy);
            Vertices.Add(newTip);
            Segments.Add(new DirectedSegment2D(oldTip, newTip));
            Angles.Add(localGrowthAngle); // These are independent of the global coordinate system
            tipSystem = new TipCoordinateSystem(newTip, globalGrowthAngle);
            CrackTipEnrichments.TipSystem = tipSystem;
        }

        public double SignedDistanceOf(XNode2D node)
        {
            return SignedDistanceOfPoint(node);
        }

        public double SignedDistanceOf(INaturalPoint2D point, XContinuumElement2D element,
             EvaluatedInterpolation2D interpolation)
        {
            return SignedDistanceOfPoint(interpolation.TransformPointNaturalToGlobalCartesian(point));
        }

        public SortedSet<ICartesianPoint2D> FindTriangleVertices(XContinuumElement2D element)
        {
            var polygon = ConvexPolygon2D.CreateUnsafe(element.Nodes);
            var triangleVertices = new SortedSet<ICartesianPoint2D>(element.Nodes, pointComparer);
            int nodesCount = element.Nodes.Count;

            foreach (var vertex in Vertices)
            {
                PolygonPointPosition position = polygon.FindRelativePositionOfPoint(vertex);
                if (position == PolygonPointPosition.Inside || position == PolygonPointPosition.OnEdge || 
                    position == PolygonPointPosition.OnVertex) triangleVertices.Add(vertex);
            }

            foreach (var crackSegment in Segments)
            {
                var segment = new LineSegment2D(crackSegment.Start, crackSegment.End);
                IReadOnlyList<ICartesianPoint2D> intersections = segment.IntersectionWith(polygon);
                foreach (var point in intersections)
                {
                    triangleVertices.Add(point);
                }
            }

            return triangleVertices;
        }

        public void UpdateEnrichments()
        {
            var bodyNodes = new HashSet<XNode2D>();
            var tipNodes = new HashSet<XNode2D>();
            tipElements.Clear();

            FindBodyAndTipNodesAndElements(bodyNodes, tipNodes);
            ApplyFixedEnrichmentArea(tipNodes, tipElements[0]);
            ResolveHeavisideEnrichmentDependencies(bodyNodes);

            ApplyEnrichmentFunctions(bodyNodes, tipNodes);
        }

        private void ApplyEnrichmentFunctions(HashSet<XNode2D> bodyNodes, HashSet<XNode2D> tipNodes)
        {
            // O(n) operation. TODO: This could be sped up by tracking the tip enriched nodes of each step.
            foreach (var node in Mesh.Vertices) node.EnrichmentItems.Remove(CrackTipEnrichments);
            foreach (var node in tipNodes)
            {
                double[] enrichmentValues = CrackTipEnrichments.EvaluateFunctionsAt(node);
                node.EnrichmentItems[CrackTipEnrichments] = enrichmentValues;
            }

            // Heaviside enrichment is never removed (unless the crack curves towards itself, but that creates a lot of
            // problems and cannot be modeled with LSM accurately). Thus there is no need to process each mesh node. 
            // TODO: It could be sped up by only updating the Heaviside enrichments of nodes that have updated body  
            // level sets, which requires tracking them.
            foreach (var node in bodyNodes)
            {
                double[] enrichmentValues = CrackBodyEnrichment.EvaluateFunctionsAt(node);
                node.EnrichmentItems[CrackBodyEnrichment] = enrichmentValues;
            }
        }

        /// <summary>
        /// If a fixed enrichment area is applied, all nodes inside a circle around the tip are enriched with tip 
        /// functions. They can still be enriched with Heaviside functions, if they do not belong to the tip 
        /// element(s).
        /// </summary>
        /// <param name="tipNodes"></param>
        /// <param name="tipElement"></param>
        private void ApplyFixedEnrichmentArea(HashSet<XNode2D> tipNodes, XContinuumElement2D tipElement)
        {
            if (tipEnrichmentAreaRadius > 0)
            {
                var enrichmentArea = new Circle2D(Vertices[Vertices.Count - 1], tipEnrichmentAreaRadius);
                foreach (var element in Mesh.FindElementsInsideCircle(enrichmentArea, tipElement))
                {
                    bool completelyInside = true;
                    foreach (var node in element.Nodes)
                    {
                        CirclePointPosition position = enrichmentArea.FindRelativePositionOfPoint(node);
                        if ((position == CirclePointPosition.Inside) || (position == CirclePointPosition.On))
                        {
                            tipNodes.Add(node);
                        }
                        else completelyInside = false;
                    }
                    if (completelyInside) element.EnrichmentItems.Add(CrackTipEnrichments);
                }

                #region alternatively
                /* // If there wasn't a need to enrich the elements, this is more performant
                foreach (var node in mesh.FindNodesInsideCircle(enrichmentArea, true, tipElement))
                {
                    tipNodes.Add(node); // Nodes of tip element(s) will not be included twice
                } */
                #endregion
            }
        }

        private ElementEnrichmentType CharacterizeElementEnrichment(XContinuumElement2D element)
        {
            var polygon = ConvexPolygon2D.CreateUnsafe(element.Nodes);
            int tipIndex = Vertices.Count - 1;

            // Check tip element
            PolygonPointPosition tipPosition = polygon.FindRelativePositionOfPoint(Vertices[tipIndex]);
            if (tipPosition == PolygonPointPosition.Inside || tipPosition == PolygonPointPosition.OnEdge ||
                    tipPosition == PolygonPointPosition.OnVertex)
            {
                PolygonPointPosition previousVertexPos = polygon.FindRelativePositionOfPoint(Vertices[tipIndex - 1]);
                if (previousVertexPos == PolygonPointPosition.Inside)
                {
                    // Problem with blending elements, if the tip element is also enriched with Heaviside. What happens
                    // after the crack tip? Based on the LSM, the signed distance of the blending element after the 
                    // crack tip should have a positive and negative region, however that element is not split by the 
                    // crack and  thus should not have discontinuity in the displacement field.
                    var builder = new StringBuilder();
                    builder.Append("Crack tip ");
                    builder.Append(Vertices[Vertices.Count - 1].ToString());
                    builder.Append(" and kink point ");
                    builder.Append(Vertices[Vertices.Count - 2].ToString());
                    builder.Append(" inside the same element with nodes: ");
                    foreach (var node in element.Nodes)
                    {
                        builder.Append(node.ToString());
                        builder.Append(' ');
                    }
                    throw new ArgumentException(builder.ToString());
                    //return ElementEnrichmentType.Both;
                }
                else return ElementEnrichmentType.Tip;
            }

            // Look at the other vertices 
            // (if a segment is inside an element, it will not be caught by checking the segment itself)
            bool previousVertexOnEdge = false;
            for (int v = 0; v < tipIndex; ++v)
            {
                PolygonPointPosition position = polygon.FindRelativePositionOfPoint(Vertices[v]);
                if (position == PolygonPointPosition.Inside) return ElementEnrichmentType.Heaviside;
                else if (position == PolygonPointPosition.OnEdge || position == PolygonPointPosition.OnVertex)
                {
                    if (previousVertexOnEdge) return ElementEnrichmentType.Heaviside;
                    else previousVertexOnEdge = true;
                }
                else previousVertexOnEdge = false;
            }

            // Look at each segment
            foreach (var crackSegment in Segments)
            {
                var segment = new LineSegment2D(crackSegment.Start, crackSegment.End);
                CartesianPoint2D intersectionPoint;
                foreach (var edge in polygon.Edges)
                {
                    LineSegment2D.SegmentSegmentPosition position = segment.IntersectionWith(edge, out intersectionPoint);
                    if (position == LineSegment2D.SegmentSegmentPosition.Intersecting)
                    {
                        // TODO: Perhaps the element should not be flagged as a Heaviside element, if the segment passes
                        // through 1 node only. To detect this, check if the intersection point coincides with an element
                        // node. If it does store it and go to the next edge. If a second intersection point (that does
                        // not coincide with the stored one) is found then it is a Heaviside element.
                        return ElementEnrichmentType.Heaviside;
                    }
                    else if (position == LineSegment2D.SegmentSegmentPosition.Overlapping)
                    {
                        return ElementEnrichmentType.Heaviside;
                    }
                }
            }

            // Then it must be a standard element
            return ElementEnrichmentType.Standard;
        }

        // Warning: TipElements must first be cleared in this iteration
        private void FindBodyAndTipNodesAndElements(HashSet<XNode2D> bodyNodes, HashSet<XNode2D> tipNodes)
        {
            var bothElements = new HashSet<XContinuumElement2D>();
            foreach (var element in Mesh.Cells)
            {
                element.EnrichmentItems.Clear();
                ElementEnrichmentType type = CharacterizeElementEnrichment(element);
                if (type == ElementEnrichmentType.Tip)
                {
                    tipElements.Add(element);
                    foreach (var node in element.Nodes) tipNodes.Add(node);
                    element.EnrichmentItems.Add(CrackTipEnrichments);
                }
                else if (type == ElementEnrichmentType.Heaviside)
                {
                    foreach (var node in element.Nodes) bodyNodes.Add(node);
                    element.EnrichmentItems.Add(CrackBodyEnrichment);
                }
                else if (type == ElementEnrichmentType.Both)
                {
                    tipElements.Add(element);
                    bothElements.Add(element);
                    foreach (var node in element.Nodes)
                    {
                        tipNodes.Add(node);
                        bodyNodes.Add(node);
                    }
                    element.EnrichmentItems.Add(CrackTipEnrichments);
                    element.EnrichmentItems.Add(CrackBodyEnrichment);
                }
            }

            // After all Heaviside nodes are aggregated remove the nodes of tip elements
            foreach (var element in tipElements)
            {
                foreach (var node in element.Nodes) bodyNodes.Remove(node);
            }
            foreach (var element in bothElements) // Re-adding these nodes afterwards is safer 
            {
                foreach (var node in element.Nodes) bodyNodes.Add(node);
            }

            ReportTipElements(tipElements);
        }

        private void FindSignedAreasOfElement(XContinuumElement2D element,
            out double positiveArea, out double negativeArea)
        {
            SortedSet<ICartesianPoint2D> triangleVertices = FindTriangleVertices(element);
            IReadOnlyList<TriangleCartesian2D> triangles = triangulator.CreateMesh(triangleVertices);
            ReportTriangulation(element, triangleVertices, triangles);

            positiveArea = 0.0;
            negativeArea = 0.0;
            foreach (var triangle in triangles)
            {
                ICartesianPoint2D v0 = triangle.Vertices[0];
                ICartesianPoint2D v1 = triangle.Vertices[1];
                ICartesianPoint2D v2 = triangle.Vertices[2];
                double area = 0.5 * Math.Abs(v0.X * (v1.Y - v2.Y) + v1.X * (v2.Y - v0.Y) + v2.X * (v0.Y - v1.Y));

                // The sign of the area can be derived from any node with signed distance != 0
                int sign = 0;
                foreach (var vertex in triangle.Vertices)
                {
                    sign = Math.Sign(SignedDistanceOfPoint(vertex));
                    if (sign != 0) break;
                }

                // If no node with non-zero signed distance is found, then find the signed distance of its centroid
                if (sign == 0)
                {
                    // Report this instance in DEBUG messages. It should not happen.
                    Console.WriteLine("--- DEBUG: Triangulation resulted in a triangle where all vertices are on the crack. ---");
                    var centroid = new CartesianPoint2D((v0.X + v1.X + v2.X) / 3.0, (v0.Y + v1.Y + v2.Y) / 3.0);
                    sign = Math.Sign(SignedDistanceOfPoint(centroid));
                }

                if (sign > 0) positiveArea += area;
                else if (sign < 0) negativeArea += area;
                else throw new Exception(
                    "Even after finding the signed distance of its centroid, the sign of the area is unidentified");
            }
        }

        private void ResolveHeavisideEnrichmentDependencies(HashSet<XNode2D> bodyNodes)
        {
            const double toleranceHeavisideEnrichmentArea = 1e-4;
            var processedElements = new Dictionary<XContinuumElement2D, Tuple<double, double>>();
            var nodesToRemove = new List<XNode2D>(); // Can't remove them while iterating the collection.
            foreach (var node in bodyNodes)
            {
                double nodePositiveArea = 0.0;
                double nodeNegativeArea = 0.0;

                foreach (var element in Mesh.FindElementsWithNode(node))
                {
                    Tuple<double, double> elementPosNegAreas;
                    bool alreadyProcessed = processedElements.TryGetValue(element, out elementPosNegAreas);
                    if (!alreadyProcessed)
                    {
                        double elementPosArea, elementNegArea;
                        FindSignedAreasOfElement(element, out elementPosArea, out elementNegArea);
                        elementPosNegAreas = new Tuple<double, double>(elementPosArea, elementNegArea);
                        processedElements[element] = elementPosNegAreas;
                    }
                    nodePositiveArea += elementPosNegAreas.Item1;
                    nodeNegativeArea += elementPosNegAreas.Item2;
                }

                if (SignedDistanceOfPoint(node) >= 0.0)
                {
                    double negativeAreaRatio = nodeNegativeArea / (nodePositiveArea + nodeNegativeArea);
                    if (negativeAreaRatio < toleranceHeavisideEnrichmentArea) nodesToRemove.Add(node);
                }
                else
                {
                    double positiveAreaRatio = nodePositiveArea / (nodePositiveArea + nodeNegativeArea);
                    if (positiveAreaRatio < toleranceHeavisideEnrichmentArea) nodesToRemove.Add(node);
                }
            }

            foreach (var node in nodesToRemove) bodyNodes.Remove(node);
        }

        private double SignedDistanceOfPoint(ICartesianPoint2D globalPoint)
        {
            if (Segments.Count == 1) return Segments[0].TransformGlobalToLocalPoint(globalPoint).Y;

            var distances = new List<double>();
            bool afterPreviousSegment = false;

            // First segment
            ICartesianPoint2D localPoint = Segments[0].TransformGlobalToLocalPoint(globalPoint);
            if (localPoint.X < Segments[0].Length) distances.Add(localPoint.Y);
            else afterPreviousSegment = true;

            // Subsequent segments
            for (int i = 1; i < Segments.Count - 1; ++i)
            {
                localPoint = Segments[i].TransformGlobalToLocalPoint(globalPoint);
                if (localPoint.X < 0.0)
                {
                    if (afterPreviousSegment)
                    {
                        // Compute the distance from the vertex between this segment and the previous
                        double dx = globalPoint.X - Vertices[i].X;
                        double dy = globalPoint.Y - Vertices[i].Y;
                        double distance = Math.Sqrt(dx * dx + dy * dy);
                        int sign = -Math.Sign(Angles[i-1]); // If growth angle > 0, the convex angle faces the positive area.
                        distances.Add(sign * distance);
                    }
                    afterPreviousSegment = false;
                }
                else if (localPoint.X <= Segments[i].Length)
                {
                    distances.Add(localPoint.Y);
                    afterPreviousSegment = false;
                }
                else afterPreviousSegment = true;
            }

            // Last segment
            int last = Segments.Count - 1;
            localPoint = Segments[last].TransformGlobalToLocalPoint(globalPoint);
            if (localPoint.X < 0.0)
            {
                if (afterPreviousSegment)
                {
                    // Compute the distance from the vertex between this segment and the previous
                    double dx = globalPoint.X - Vertices[last].X;
                    double dy = globalPoint.Y - Vertices[last].Y;
                    double distance = Math.Sqrt(dx * dx + dy * dy);
                    int sign = -Math.Sign(Angles[last - 1]); // If growth angle > 0, the convex angle faces the positive area.
                    distances.Add(sign * distance);
                }
                afterPreviousSegment = false;
            }
            else distances.Add(localPoint.Y);

            return distances[MathUtilities.IndexOfMinAbs(distances)];
        }

        [ConditionalAttribute("DEBUG")]
        private void ReportTipElements(IReadOnlyList<XContinuumElement2D> tipElements)
        {
            if (!reports) return;
            Console.WriteLine("------ DEBUG: TIP ELEMENTS/ ------");
            if (tipElements.Count < 1) throw new Exception("No tip element found");
            Console.WriteLine("Tip elements:");
            for (int e = 0; e < tipElements.Count; ++e)
            {
                Console.WriteLine("Tip element " + e + " with nodes: ");
                foreach (var node in tipElements[e].Nodes)
                {
                    Console.WriteLine(node);
                }
            }
            Console.WriteLine("------ /DEBUG: TIP ELEMENTS ------");
        }

        [ConditionalAttribute("DEBUG")]
        private void ReportTriangulation(XContinuumElement2D element, SortedSet<ICartesianPoint2D> triangleVertices, 
            IReadOnlyList<TriangleCartesian2D> triangles)
        {
            if (!reports) return;
            Console.WriteLine("------ DEBUG: TRIANGULATION/ ------");
            Console.WriteLine("Element with nodes: ");
            foreach (var node in element.Nodes) 
            {
                Console.Write(node);
                Console.Write(' ');
            }

            Console.WriteLine("\nVertices of triangular mesh: ");
            foreach (var vertex in triangleVertices)
            {
                Console.Write(vertex);
                Console.Write(' ');
            }

            Console.WriteLine("\nTriangles: ");
            foreach (var triangle in triangles)
            {
                Console.WriteLine(triangle);
            }

            Console.WriteLine("------ /DEBUG: TRIANGULATION ------");
        }

        public void Propagate(IDofOrderer dofOrderer, Vector totalFreeDisplacements, Vector totalConstrainedDisplacements)
        {
            throw new NotImplementedException();
        }

        public void InitializeGeometry(PolyLine2D initialCrack)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Represents the type of enrichment that will be applied to all nodes of the element. 
        /// </summary>
        private enum ElementEnrichmentType { Standard, Heaviside, Tip, Both }
    }
}
