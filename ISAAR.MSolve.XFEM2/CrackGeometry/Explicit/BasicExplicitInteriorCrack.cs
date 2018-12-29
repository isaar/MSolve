using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;


namespace ISAAR.MSolve.XFEM.CrackGeometry.Explicit
{
    class BasicExplicitInteriorCrack
    {
        private static readonly bool reports = false;
        private static readonly IComparer<ICartesianPoint2D> pointComparer = new Point2DComparerXMajor();

        private readonly double enrichmentRadiusOverElementSize;
        private readonly CartesianTriangulator triangulator;

        // TODO: Not too fond of the setters, but at least the enrichments are immutable. Perhaps I can pass their
        // parameters to a CrackDescription builder and construct them there, without involving the user 
        // (given how easy it is to forget the setters, it is a must).
        public IMesh2D<XNode2D, XContinuumElement2D> Mesh { get; set; }
        public CrackBodyEnrichment2D CrackBodyEnrichment { get; set; }
        public CrackTipEnrichments2D StartTipEnrichments { get; set; }
        public CrackTipEnrichments2D EndTipEnrichments { get; set; }

        private LinkedList<ICartesianPoint2D> Vertices { get; }
        private LinkedList<DirectedSegment2D> Segments { get; }
        // Angles[i-1] is the angle of segment i w.r.t segment i-1, aka the crack growth angle.
        private LinkedList<double> Angles { get; }

        private TipCoordinateSystem startTipSystem;
        private TipCoordinateSystem endTipSystem;
        private List<XContinuumElement2D> startTipElements;
        private List<XContinuumElement2D> endTipElements;

        public BasicExplicitInteriorCrack(double enrichmentRadiusOverElementSize = 0.0)
        {
            this.enrichmentRadiusOverElementSize = enrichmentRadiusOverElementSize;
            this.triangulator = new CartesianTriangulator();
            this.startTipElements = new List<XContinuumElement2D>();
            this.endTipElements = new List<XContinuumElement2D>();

            Vertices = new LinkedList<ICartesianPoint2D>();
            Segments = new LinkedList<DirectedSegment2D>();
            Angles = new LinkedList<double>();
        }

        public ICartesianPoint2D GetCrackTip(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Start) return Vertices.First.Value;
            else if (tipPosition == CrackTipPosition.End) return Vertices.Last.Value;
            else throw new ArgumentException("Invalid tip position");
        }

        public TipCoordinateSystem GetTipSystem(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Start) return startTipSystem;
            else if (tipPosition == CrackTipPosition.End) return endTipSystem;
            else throw new ArgumentException("Invalid tip position");
        }

        public IReadOnlyList<XContinuumElement2D> GetTipElements(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Start) return startTipElements;
            else if (tipPosition == CrackTipPosition.End) return endTipElements;
            else throw new ArgumentException("Invalid tip position");
        }

        public void InitializeGeometry(ICartesianPoint2D startTip, ICartesianPoint2D endTip)
        {
            double dx = endTip.X - startTip.X;
            double dy = endTip.Y - startTip.Y;
            double tangentSlope = Math.Atan2(dy, dx);
            endTipSystem = new TipCoordinateSystem(endTip, tangentSlope);
            startTipSystem = new TipCoordinateSystem(startTip, tangentSlope + Math.PI);

            Vertices.AddLast(startTip);
            Vertices.AddLast(endTip);
            Segments.AddLast(new DirectedSegment2D(startTip, endTip));
        }

        public void UpdateGeometry(double localGrowthAngleStart, double growthLengthStart,
            double localGrowthAngleEnd, double growthLengthEnd)
        {
            // End tip
            double globalAngleEnd = AngleUtilities.Wrap(localGrowthAngleEnd + endTipSystem.RotationAngle);
            double dxEnd = growthLengthEnd * Math.Cos(globalAngleEnd);
            double dyEnd = growthLengthEnd * Math.Sin(globalAngleEnd);
            var oldEndTip = Vertices.Last.Value; ;
            var newEndTip = new CartesianPoint2D(oldEndTip.X + dxEnd, oldEndTip.Y + dyEnd);
            Vertices.AddLast(newEndTip);
            Segments.AddLast(new DirectedSegment2D(oldEndTip, newEndTip));
            Angles.AddLast(localGrowthAngleEnd); // These are independent of the global coordinate system
            endTipSystem = new TipCoordinateSystem(newEndTip, globalAngleEnd);

            //StartTip
            double globalAngleStart = AngleUtilities.Wrap(localGrowthAngleStart + startTipSystem.RotationAngle);
            double dxStart = growthLengthStart * Math.Cos(globalAngleStart);
            double dyStart = growthLengthStart * Math.Sin(globalAngleStart);
            var oldStartTip = Vertices.First.Value;
            var newStartTip = new CartesianPoint2D(oldStartTip.X + dxEnd, oldStartTip.Y + dyEnd);
            Vertices.AddFirst(newStartTip);
            Segments.AddFirst(new DirectedSegment2D(newStartTip, oldStartTip));
            Angles.AddFirst(-localGrowthAngleEnd); // From new to old start segment, use the opposite angle
            endTipSystem = new TipCoordinateSystem(newStartTip, globalAngleStart);
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
            var startTipNodes = new HashSet<XNode2D>();
            var endTipNodes = new HashSet<XNode2D>();
            startTipElements.Clear();
            endTipElements.Clear();

            FindBodyAndTipNodesAndElements(bodyNodes, startTipNodes, endTipNodes);
            ApplyFixedEnrichmentArea(Vertices.First.Value, startTipElements[0], startTipNodes, StartTipEnrichments);
            ApplyFixedEnrichmentArea(Vertices.Last.Value, endTipElements[0], endTipNodes, EndTipEnrichments);
            ResolveHeavisideEnrichmentDependencies(bodyNodes);

            ApplyEnrichmentFunctions(bodyNodes, startTipNodes, endTipNodes);
        }

        private void ApplyEnrichmentFunctions(HashSet<XNode2D> bodyNodes, HashSet<XNode2D> startTipNodes, 
            HashSet<XNode2D> endTipNodes)
        {
            // O(n) operation. TODO: This could be sped up by tracking the tip enriched nodes of each step.
            foreach (var node in Mesh.Vertices)
            {
                node.EnrichmentItems.Remove(StartTipEnrichments);
                node.EnrichmentItems.Remove(EndTipEnrichments);
            }

            foreach (var node in startTipNodes)
            {
                double[] enrichmentValues = StartTipEnrichments.EvaluateFunctionsAt(node);
                node.EnrichmentItems[StartTipEnrichments] = enrichmentValues;
            }

            foreach (var node in endTipNodes)
            {
                double[] enrichmentValues = EndTipEnrichments.EvaluateFunctionsAt(node);
                node.EnrichmentItems[EndTipEnrichments] = enrichmentValues;
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
        private void ApplyFixedEnrichmentArea(ICartesianPoint2D crackTip, XContinuumElement2D tipElement, 
            HashSet<XNode2D> tipNodes, IEnrichmentItem2D tipEnrichments)
        {
            if (enrichmentRadiusOverElementSize > 0)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(tipElement.Nodes);
                double elementArea = outline.ComputeArea();
                double radius = enrichmentRadiusOverElementSize * Math.Sqrt(elementArea);
                var enrichmentArea = new Circle2D(crackTip, radius);

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
                    if (completelyInside) element.EnrichmentItems.Add(tipEnrichments);
                }
            }
        }

        private void FindBodyAndTipNodesAndElements(HashSet<XNode2D> bodyNodes, HashSet<XNode2D> startTipNodes, 
            HashSet<XNode2D> endTipNodes)
        {
            var bodyElements = new HashSet<XContinuumElement2D>();
            foreach (var element in Mesh.Cells)
            {
                element.EnrichmentItems.Clear();
                bool isCut = IsCutElement(element);
                bool containsStartTip = IsTipElement(element, CrackTipPosition.Start);
                bool containsEndTip = IsTipElement(element, CrackTipPosition.End);
                
                if (containsStartTip)
                {
                    startTipElements.Add(element);
                    foreach (var node in element.Nodes) startTipNodes.Add(node);
                    element.EnrichmentItems.Add(StartTipEnrichments);
                }
                if (containsEndTip)
                {
                    endTipElements.Add(element);
                    foreach (var node in element.Nodes) endTipNodes.Add(node);
                    element.EnrichmentItems.Add(EndTipEnrichments);
                }
                if (isCut)
                {
                    bodyElements.Add(element);
                    foreach (var node in element.Nodes) bodyNodes.Add(node);
                    element.EnrichmentItems.Add(CrackBodyEnrichment);
                }
            }

            // After all Heaviside nodes are aggregated remove the nodes of tip elements. 
            // TODO: Handle the case where the tip element contains 2 vertices. In this case the Heaviside enrichment should not be removed.
            foreach (var element in startTipElements)
            {
                foreach (var node in element.Nodes) bodyNodes.Remove(node);
            }
            foreach (var element in endTipElements)
            {
                foreach (var node in element.Nodes) bodyNodes.Remove(node);
            }
            //ReportTipElements(startTipElements, endTipElements);
        }

        private bool IsTipElement(XContinuumElement2D element, CrackTipPosition tipPosition)
        {
            ICartesianPoint2D crackTip, adjacentVertex;
            if (tipPosition == CrackTipPosition.Start)
            {
                crackTip = Vertices.First.Value;
                adjacentVertex = Vertices.First.Next.Value;
            }
            else if (tipPosition == CrackTipPosition.End)
            {
                crackTip = Vertices.Last.Value;
                adjacentVertex = Vertices.Last.Previous.Value;
            }
            else throw new ArgumentException("Tip position can be either start or end");

            var polygon = ConvexPolygon2D.CreateUnsafe(element.Nodes);
            PolygonPointPosition relativeTipPos = polygon.FindRelativePositionOfPoint(crackTip);
            if (relativeTipPos == PolygonPointPosition.Inside || relativeTipPos == PolygonPointPosition.OnEdge ||
                    relativeTipPos == PolygonPointPosition.OnVertex)
            {
                PolygonPointPosition previousVertexPos = polygon.FindRelativePositionOfPoint(adjacentVertex);
                if (previousVertexPos == PolygonPointPosition.Inside)
                {
                    throw new NotImplementedException("Problem with blending elements, if the tip element is also " +
                        "enriched with Heaviside. What happens after the crack tip? Based on the LSM, the signed " +
                        "distance of the blending element after the crack tip should have a positive and negative " +
                        "region, however that element is not split by the crack and  thus should not have " +
                        "discontinuity in the displacement field");
                    //return ElementEnrichmentType.Both;
                }
                return true;
            }
            return false;
        }

        private bool IsCutElement(XContinuumElement2D element)
        {
            var polygon = ConvexPolygon2D.CreateUnsafe(element.Nodes);

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
                        return true;
                    }
                    else if (position == LineSegment2D.SegmentSegmentPosition.Overlapping)
                    {
                        return true;
                    }
                }
            }

            // Look at the vertices 
            // (if a segment is entirely inside an element, it will not be caught by checking the segment itself)
            bool previousVertexOnEdge = false;
            LinkedListNode<ICartesianPoint2D> currentNode = Vertices.First.Next;
            LinkedListNode<ICartesianPoint2D> lastNode = Vertices.Last;
            while(currentNode != lastNode)
            {
                PolygonPointPosition position = polygon.FindRelativePositionOfPoint(currentNode.Value);
                if (position == PolygonPointPosition.Inside) return true;
                else if (position == PolygonPointPosition.OnEdge || position == PolygonPointPosition.OnVertex)
                {
                    if (previousVertexOnEdge) return true;
                    else previousVertexOnEdge = true;
                }
                else previousVertexOnEdge = false;
                currentNode = currentNode.Next;
            }

            return false;
        }

        private double SignedDistanceOfPoint(ICartesianPoint2D globalPoint)
        {
            if (Segments.Count == 1) return Segments.First.Value.TransformGlobalToLocalPoint(globalPoint).Y;

            var distances = new List<double>();
            bool afterPreviousSegment = false;
            LinkedList<DirectedSegment2D>.Enumerator segmentEnumerator = Segments.GetEnumerator();
            LinkedList<ICartesianPoint2D>.Enumerator vertexEnumerator = Vertices.GetEnumerator();
            LinkedList<double>.Enumerator angleEnumerator = Angles.GetEnumerator();

            // First segment
            segmentEnumerator.MoveNext(); // Do not advance angleEnumerator yet, since there is 1 less angle than segments.
            vertexEnumerator.MoveNext();
            ICartesianPoint2D localPoint = segmentEnumerator.Current.TransformGlobalToLocalPoint(globalPoint);
            if (localPoint.X < segmentEnumerator.Current.Length) distances.Add(localPoint.Y);
            else afterPreviousSegment = true;

            // Subsequent segments
            for (int i = 1; i < Segments.Count - 1; ++i)
            {
                segmentEnumerator.MoveNext();
                vertexEnumerator.MoveNext();
                angleEnumerator.MoveNext(); 
                localPoint = segmentEnumerator.Current.TransformGlobalToLocalPoint(globalPoint);
                if (localPoint.X < 0.0)
                {
                    if (afterPreviousSegment)
                    {
                        // Compute the distance from the vertex between this segment and the previous
                        double dx = globalPoint.X - vertexEnumerator.Current.X;
                        double dy = globalPoint.Y - vertexEnumerator.Current.Y;
                        double distance = Math.Sqrt(dx * dx + dy * dy);
                        int sign = -Math.Sign(angleEnumerator.Current); // If growth angle > 0, the convex angle faces the positive area.
                        distances.Add(sign * distance);
                    }
                    afterPreviousSegment = false;
                }
                else if (localPoint.X <= segmentEnumerator.Current.Length)
                {
                    distances.Add(localPoint.Y);
                    afterPreviousSegment = false;
                }
                else afterPreviousSegment = true;
            }

            // Last segment
            segmentEnumerator.MoveNext();
            vertexEnumerator.MoveNext();
            angleEnumerator.MoveNext();
            localPoint = segmentEnumerator.Current.TransformGlobalToLocalPoint(globalPoint);
            if (localPoint.X < 0.0)
            {
                if (afterPreviousSegment)
                {
                    // Compute the distance from the vertex between this segment and the previous
                    double dx = globalPoint.X - vertexEnumerator.Current.X;
                    double dy = globalPoint.Y - vertexEnumerator.Current.Y;
                    double distance = Math.Sqrt(dx * dx + dy * dy);
                    int sign = -Math.Sign(angleEnumerator.Current); // If growth angle > 0, the convex angle faces the positive area.
                    distances.Add(sign * distance);
                }
                afterPreviousSegment = false;
            }
            else distances.Add(localPoint.Y);

            return distances[MathUtilities.IndexOfMinAbs(distances)];
        }

        private void FindSignedAreasOfElement(XContinuumElement2D element,
            out double positiveArea, out double negativeArea)
        {
            SortedSet<ICartesianPoint2D> triangleVertices = FindTriangleVertices(element);
            IReadOnlyList<TriangleCartesian2D> triangles = triangulator.CreateMesh(triangleVertices);

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

        [ConditionalAttribute("DEBUG")]
        private void ReportTipElements(IReadOnlyList<XContinuumElement2D> startTipElements, 
            IReadOnlyList<XContinuumElement2D> endTipElements)
        {
            if (!reports) return;
            Console.WriteLine("------ DEBUG/ ------");

            // Start tip
            if (startTipElements.Count < 1) throw new Exception("No start tip element found");
            Console.WriteLine("Tip elements:");
            for (int e = 0; e < startTipElements.Count; ++e)
            {
                Console.WriteLine("Tip element " + e + " with nodes: ");
                foreach (var node in startTipElements[e].Nodes)
                {
                    Console.WriteLine(node);
                }
            }

            // End tip
            if (endTipElements.Count < 1) throw new Exception("No end tip element found");
            Console.WriteLine("End tip elements:");
            for (int e = 0; e < endTipElements.Count; ++e)
            {
                Console.WriteLine("Tip element " + e + " with nodes: ");
                foreach (var node in endTipElements[e].Nodes)
                {
                    Console.WriteLine(node);
                }
            }

            Console.WriteLine("------ /DEBUG ------");
        }

        public IReadOnlyList<ICartesianPoint2D> GetCrackTips()
        {
            throw new NotImplementedException();
        }

        public void Propagate(IDofOrderer dofOrderer, Vector totalFreeDisplacements, Vector totalConstrainedDisplacements)
        {
            throw new NotImplementedException();
        }
    }
}
