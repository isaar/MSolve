using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Commons;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Geometry.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

// TODO: Consider removing the bookkeeping of enrichment items in elements. It creates a lot of opportunities for mistakes.
//       Could the enrichment type of an element be infered by just looking at its nodes, without storing state. Could it be
//       cached for efficiency and how would the cache be updated safely?
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit
{
    // For mouth cracks only (single tip). Warning: may misclassify elements as tip elements, causing gross errors.
    class BasicCrackLSM
    {
        private static readonly bool reports = false;
        private static readonly IComparer<CartesianPoint> pointComparer = new Point2DComparerXMajor();

        private readonly double tipEnrichmentAreaRadius;
        private readonly Triangulator2D<CartesianPoint> triangulator;

        public CartesianPoint CrackMouth { get; private set; }

        // TODO: Not too fond of the setters, but at least the enrichments are immutable. Perhaps I can pass their
        // parameters to a CrackDescription builder and construct them there, without involving the user.
        public IMesh2D<XNode, XContinuumElement2D> Mesh { get; set; }
        public CrackBodyEnrichment2D CrackBodyEnrichment { get; set; }
        public CrackTipEnrichments2D CrackTipEnrichments { get; set; }

        public ISet<XNode> CrackBodyNodesAll => throw new NotImplementedException();
        public ISet<XNode> CrackTipNodesNew => throw new NotImplementedException();
        public ISet<XNode> CrackBodyNodesNew => throw new NotImplementedException();
        public ISet<XNode> CrackTipNodesOld => throw new NotImplementedException();

        private readonly Dictionary<XNode, double> levelSetsBody;
        private readonly Dictionary<XNode, double> levelSetsTip;

        private CartesianPoint crackTip;
        private TipCoordinateSystem tipSystem;
        private List<XContinuumElement2D> tipElements;

        public BasicCrackLSM(double tipEnrichmentAreaRadius = 0.0)
        {
            this.tipEnrichmentAreaRadius = tipEnrichmentAreaRadius;
            this.triangulator = new Triangulator2D<CartesianPoint>((x, y) => new CartesianPoint(x, y));
            this.tipElements = new List<XContinuumElement2D>();
            levelSetsBody = new Dictionary<XNode, double>();
            levelSetsTip = new Dictionary<XNode, double>();
        }

        public CartesianPoint GetCrackTip(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Single) return crackTip;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        public TipCoordinateSystem GetTipSystem(CrackTipPosition tip)
        {
            if (tip == CrackTipPosition.Single) return tipSystem;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        public IReadOnlyList<XContinuumElement2D> GetTipElements(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Single) return tipElements;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        public void InitializeGeometry(CartesianPoint crackMouth, CartesianPoint crackTip)
        {
            CrackMouth = crackMouth;
            var segment = new DirectedSegment2D(crackMouth, crackTip);

            double tangentX = crackTip.X - crackMouth.X;
            double tangentY = crackTip.Y - crackMouth.Y;
            double length = Math.Sqrt(tangentX * tangentX + tangentY * tangentY);
            double tangentSlope = Math.Atan2(tangentY, tangentX);
            this.crackTip = crackTip;
            tipSystem = new TipCoordinateSystem(crackTip, tangentSlope);

            tangentX /= length;
            tangentY /= length;

            foreach (XNode node in Mesh.Nodes)
            {
                levelSetsBody[node] = segment.SignedDistanceOf(node);
                levelSetsTip[node] = (node.X - crackTip.X) * tangentX + (node.Y - crackTip.Y) * tangentY;
            }
        }

        public void UpdateGeometry(double localGrowthAngle, double growthLength)
        {
            double globalGrowthAngle = MathUtilities.WrapAngle(localGrowthAngle + tipSystem.RotationAngle);
            double dx = growthLength * Math.Cos(globalGrowthAngle);
            double dy = growthLength * Math.Sin(globalGrowthAngle);
            double unitDx = dx / growthLength;
            double unitDy = dy / growthLength;

            var oldTip = crackTip;
            var newTip = new CartesianPoint(oldTip.X + dx, oldTip.Y + dy);
            crackTip = newTip;
            tipSystem = new TipCoordinateSystem(newTip, globalGrowthAngle);

            var newSegment = new DirectedSegment2D(oldTip, newTip);
            foreach (XNode node in Mesh.Nodes)
            {
                // Rotate the ALL tip level sets towards the new tip and then advance them
                double rotatedTipLevelSet = (node.X - crackTip.X) * unitDx + (node.Y - crackTip.Y) * unitDy;
                levelSetsTip[node] = rotatedTipLevelSet - newSegment.Length;
                 
                if (rotatedTipLevelSet > 0.0) // Only some body level sets are updated (See Stolarska 2001) 
                {
                    levelSetsBody[node] = newSegment.SignedDistanceOf(node);
                }
            }
        }

        /// <summary>
        /// Warning: with narrow band this should throw an exception if the node is not tracked.
        /// </summary>
        /// <param name="node"></param>
        /// <returns></returns>
        public double SignedDistanceOf(XNode node)
        {
            return levelSetsBody[node];
        }

        /// <summary>
        /// Warning: with narrow band this should throw an exception if the element/nodes are not tracked.
        /// </summary>
        /// <param name="point"></param>
        /// <param name="elementNodes"></param>
        /// <param name="interpolation"></param>
        /// <returns></returns>
        public double SignedDistanceOf(NaturalPoint point, XContinuumElement2D element,
             EvalInterpolation2D interpolation)
        {
            double signedDistance = 0.0;
            for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
            {
                signedDistance += interpolation.ShapeFunctions[nodeIdx] * levelSetsBody[element.Nodes[nodeIdx]];
            }
            return signedDistance;
        }

        public Tuple<double, double> SignedDistanceGradientThrough(NaturalPoint point,
            IReadOnlyList<XNode> elementNodes, EvalInterpolation2D interpolation)
        {
            double gradientX = 0.0;
            double gradientY = 0.0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                double levelSet = levelSetsBody[elementNodes[nodeIdx]];
                gradientX += interpolation.ShapeGradientsCartesian[nodeIdx, 0] * levelSet;
                gradientY += interpolation.ShapeGradientsCartesian[nodeIdx, 1] * levelSet;
            }
            return new Tuple<double, double>(gradientX, gradientY);
        }

        public SortedSet<CartesianPoint> FindTriangleVertices(XContinuumElement2D element)
        {
            var triangleVertices = new SortedSet<CartesianPoint>(element.Nodes, pointComparer);
            int nodesCount = element.Nodes.Count;
            ElementEnrichmentType type = CharacterizeElementEnrichment(element);
            if (type != ElementEnrichmentType.Standard)
            {
                // Find the intersections between element edges and the crack. TODO: See Serafeim's Msc Thesis for a correct procedure.
                for (int i = 0; i < nodesCount; ++i)
                {
                    XNode node1 = element.Nodes[i];
                    XNode node2 = element.Nodes[(i + 1) % nodesCount];
                    double levelSet1 = levelSetsBody[node1];
                    double levelSet2 = levelSetsBody[node2];

                    if (levelSet1 * levelSet2 < 0.0) 
                    {
                        // The intersection point between these nodes can be found using the linear interpolation, see 
                        // Sukumar 2001
                        double k = -levelSet1 / (levelSet2 - levelSet1);
                        double x = node1.X + k * (node2.X - node1.X);
                        double y = node1.Y + k * (node2.Y - node1.Y);

                        // TODO: For the tip element one intersection point is on the crack extension and does not  
                        // need to be added. It is not wrong though.
                        triangleVertices.Add(new CartesianPoint(x, y)); 
                    }
                    else if (levelSet1 == 0.0) triangleVertices.Add(node1); // TODO: perhaps some tolerance is needed.
                    else if (levelSet2 == 0.0) triangleVertices.Add(node2);
                }

                if (type == ElementEnrichmentType.Tip) triangleVertices.Add(crackTip);
            }
            return triangleVertices;
        }

        // The tip enrichments are cleared and reapplied at each call. In constrast, nodes previously enriched with Heavise will 
        // continue to do so, with the addition of newly Heaviside enriched nodes. If the caller needs the nodes to be cleared of 
        // Heaviside enrichments he must do so himself.
        public void UpdateEnrichments() 
        {
            var bodyNodes = new HashSet<XNode>();
            var tipNodes = new HashSet<XNode>();
            tipElements.Clear();

            FindBodyAndTipNodesAndElements(bodyNodes, tipNodes, tipElements);
            ApplyFixedEnrichmentArea(tipNodes, tipElements[0]);
            ResolveHeavisideEnrichmentDependencies(bodyNodes);

            ApplyEnrichmentFunctions(bodyNodes, tipNodes); 
        }

        private void ApplyEnrichmentFunctions(HashSet<XNode> bodyNodes, HashSet<XNode> tipNodes)
        {
            // O(n) operation. TODO: This could be sped up by tracking the tip enriched nodes of each step.
            foreach (var node in Mesh.Nodes) node.EnrichmentItems.Remove(CrackTipEnrichments);
            foreach (var node in tipNodes)
            {
                double[] enrichmentValues = CrackTipEnrichments.EvaluateFunctionsAt(node);
                node.EnrichmentItems[CrackTipEnrichments] = enrichmentValues;
            }

            // Heaviside enrichment is never removed (unless the crack curves towards itself, but that creates a lot of
            // problems and cannot be modeled with LSM accurately). Thus there is no need to process each mesh node. 
            // TODO: It could be sped up by only updating the Heaviside enrichments of nodes that have updated body  
            //      level sets, which requires tracking them.
            // TODO: should I also clear and reapply all Heaviside enrichments? It is safer and might be useful for e.g. 
            //      reanalysis. Certainly I must not clear all node enrichments, as they may include material interfaces etc.
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
        private void ApplyFixedEnrichmentArea(HashSet<XNode> tipNodes, XContinuumElement2D tipElement)
        {
            if (tipEnrichmentAreaRadius > 0)
            {
                var enrichmentArea = new Circle2D(crackTip, tipEnrichmentAreaRadius);
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
            double minBodyLevelSet = double.MaxValue;
            double maxBodyLevelSet = double.MinValue;
            double minTipLevelSet = double.MaxValue;
            double maxTipLevelSet = double.MinValue;

            foreach (XNode node in element.Nodes)
            {
                double bodyLevelSet = levelSetsBody[node];
                double tipLevelSet = levelSetsTip[node];
                if (bodyLevelSet < minBodyLevelSet) minBodyLevelSet = bodyLevelSet;
                if (bodyLevelSet > maxBodyLevelSet) maxBodyLevelSet = bodyLevelSet;
                if (tipLevelSet < minTipLevelSet) minTipLevelSet = tipLevelSet;
                if (tipLevelSet > maxTipLevelSet) maxTipLevelSet = tipLevelSet;
            }

            // Warning: This criterion might give false positives for tip elements (see Serafeim's thesis for details)
            if (minBodyLevelSet * maxBodyLevelSet <= 0.0)
            {
                if (minTipLevelSet * maxTipLevelSet <= 0) return ElementEnrichmentType.Tip;
                else if (maxTipLevelSet < 0) return ElementEnrichmentType.Heaviside;
            }
            return ElementEnrichmentType.Standard;
        }

        private void FindBodyAndTipNodesAndElements(HashSet<XNode> bodyNodes, HashSet<XNode> tipNodes, 
            List<XContinuumElement2D> tipElements)
        {
            foreach (var element in Mesh.Elements)
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
                    // Cut elements next to tip elements, they will be enriched with both Heaviside and tip functions. If all 
                    // nodes were enriched with tip functions, the element would not be enriched with Heaviside, but then it 
                    // would be a tip element and not fall under this case.
                    element.EnrichmentItems.Add(CrackBodyEnrichment);
                }
            }
            foreach (var node in tipNodes) bodyNodes.Remove(node); // tip element's nodes are not enriched with Heaviside

            ReportTipElements(tipElements);
        }

        private void FindSignedAreasOfElement(XContinuumElement2D element, 
            out double positiveArea, out double negativeArea)
        {
            SortedSet<CartesianPoint> triangleVertices = FindTriangleVertices(element);
            IReadOnlyList<Triangle2D<CartesianPoint>> triangles = triangulator.CreateMesh(triangleVertices);

            positiveArea = 0.0;
            negativeArea = 0.0;
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
                    if (vertex is XNode)
                    {
                        sign = Math.Sign(levelSetsBody[(XNode)vertex]);
                        if (sign != 0) break;
                    }
                }

                // If no node with non-zero body level set is found, then find the body level set of its centroid
                if (sign == 0)
                {
                    // Report this instance in DEBUG messages. It should not happen with linear level sets and only 1 crack.
                    if (reports)
                    {
                        Console.WriteLine("--- DEBUG: Triangulation resulted in a triangle where no vertex is an element node. ---");
                    }
                   
                    
                    var centroid = new CartesianPoint((v0.X + v1.X + v2.X) / 3.0, (v0.Y + v1.Y + v2.Y) / 3.0);
                    NaturalPoint centroidNatural = element.Interpolation.
                        CreateInverseMappingFor(element.Nodes).TransformPointCartesianToNatural(centroid);
                    EvalInterpolation2D centroidInterpolation = 
                        element.Interpolation.EvaluateAllAt(element.Nodes, centroidNatural);
                    sign = Math.Sign(SignedDistanceOf(centroidNatural, element, centroidInterpolation));
                }

                if (sign > 0) positiveArea += area;
                else if (sign < 0) negativeArea += area;
                else throw new Exception(
                    "Even after finding the signed distance of its centroid, the sign of the area is unidentified");
            }
        }

        private void ResolveHeavisideEnrichmentDependencies(HashSet<XNode> bodyNodes)
        {
            const double toleranceHeavisideEnrichmentArea = 1e-4;
            var processedElements = new Dictionary<XContinuumElement2D, Tuple<double, double>>();
            var nodesToRemove = new List<XNode>(); // Can't remove them while iterating the collection.
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

                if (levelSetsBody[node] >= 0.0)
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
        private void ReportTipElements(IReadOnlyList<XContinuumElement2D> tipElements)
        {
            if (!reports) return;
            Console.WriteLine("------ DEBUG/ ------");
            if (tipElements.Count < 1) throw new Exception("No tip element found");
            Console.WriteLine("Tip elements:");
            for (int e = 0; e < tipElements.Count; ++e)
            {
                Console.WriteLine("Tip element " + e + " with nodes: ");
                foreach (var node in tipElements[e].Nodes)
                {
                    Console.Write(node);
                    Console.Write(" - body level set = " + levelSetsBody[node]);
                    Console.WriteLine(" - tip level set = " + levelSetsTip[node]);
                }
            }
            Console.WriteLine("------ /DEBUG ------");
        }

        public IReadOnlyList<CartesianPoint> GetCrackTips()
        {
            throw new NotImplementedException();
        }

        public void Propagate(IDofOrderer dofOrderer, Vector totalFreeDisplacements, Vector totalConstrainedDisplacements)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Represents the type of enrichment that will be applied to all nodes of the element. In LSM with linear 
        /// interpolations, an element enriched with tip functions does not need to be enriched with Heaviside too. 
        /// This is because, even if there are kinks inside the element, the linear interpolation cannot reproduce them.
        /// </summary>
        private enum ElementEnrichmentType { Standard, Heaviside, Tip } 
    }
}
