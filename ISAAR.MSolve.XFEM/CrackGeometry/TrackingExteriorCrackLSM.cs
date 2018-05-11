using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

//TODO: replace BasicCrackLSM with this one after testing is done
//TODO: Consider removing the bookkeeping of enrichment items in elements. It creates a lot of opportunities for mistakes.
//      Could the enrichment type of an element be infered by just looking at its nodes, without storing state? Should it be
//      cached for efficiency and how would the cache be updated safely?
//TODO: using narrow band, a superset of tipNodes and bodyNodes would include those sets. Need core class SuperSet.
//TODO: perhaps the bookkeeping of nodes and elements can be done by a dedicated class. Narrow banding would then be implemented
//      there. In general, this is a god class and should be broken down to smaller ones.
//TODO: Crack tips should be handled differently than using enums. Interior and exterior cracks should compose their common 
//      dedicated strategy classes with their common functionality and expose appropriate properties for the crack tip data.
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    /// <summary>
    /// Warning: may misclassify elements as tip elements, causing gross errors.
    /// </summary>
    class TrackingExteriorCrackLSM: IExteriorCrack
    {
        private static readonly TipDetectionScheme tipDetection = TipDetectionScheme.Explicit;
        private static readonly bool reports = false;
        private static readonly IComparer<ICartesianPoint2D> pointComparer = new Point2DComparerXMajor();
        private readonly Dictionary<XNode2D, double> levelSetsBody;
        private readonly Dictionary<XNode2D, double> levelSetsTip;
        private readonly double tipEnrichmentAreaRadius;
        private readonly List<XContinuumElement2D> tipElements; // Ideally there is only 1, but what if the tip falls on the edge bewteen elements?
        
        private readonly CartesianTriangulator triangulator;
        private ICartesianPoint2D crackTip;
        private TipCoordinateSystem tipSystem;

        public TrackingExteriorCrackLSM(double tipEnrichmentAreaRadius = 0.0)
        {
            this.levelSetsBody = new Dictionary<XNode2D, double>();
            this.levelSetsTip = new Dictionary<XNode2D, double>();
            this.tipEnrichmentAreaRadius = tipEnrichmentAreaRadius;
            this.tipElements = new List<XContinuumElement2D>();
            this.CrackBodyNodesAll = new HashSet<XNode2D>();
            this.CrackBodyNodesNew = new HashSet<XNode2D>();
            this.CrackTipNodesNew = new HashSet<XNode2D>();
            this.CrackTipNodesOld = new HashSet<XNode2D>();
            this.triangulator = new CartesianTriangulator();
        }

        // TODO: Not too fond of the setters, but at least the enrichments are immutable. Perhaps I can pass their
        // parameters to a CrackDescription builder and construct them there, without involving the user.
        public CrackBodyEnrichment2D CrackBodyEnrichment { get; set; }
        public CrackTipEnrichments2D CrackTipEnrichments { get; set; }
        public ISet<XNode2D> CrackBodyNodesAll { get; private set; } // TODO: a TreeSet might be better if set intersections are applied
        public ISet<XNode2D> CrackBodyNodesNew { get; private set; }
        public ISet<XNode2D> CrackTipNodesNew { get; private set; }
        public ISet<XNode2D> CrackTipNodesOld { get; private set; }
        public ICartesianPoint2D CrackMouth { get; private set; }
        public IReadOnlyDictionary<XNode2D, double> LevelSetsBody { get { return levelSetsBody; } }
        public IReadOnlyDictionary<XNode2D, double> LevelSetsTip { get { return levelSetsTip; } }
        public LsmLogger Logger { get; set; }
        public IMesh2D<XNode2D, XContinuumElement2D> Mesh { get; set; }

        public ICartesianPoint2D GetCrackTip(CrackTipPosition tipPosition) 
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

        public void InitializeGeometry(ICartesianPoint2D crackMouth, ICartesianPoint2D crackTip)
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

            foreach (XNode2D node in Mesh.Vertices)
            {
                levelSetsBody[node] = segment.SignedDistanceOf(node);
                levelSetsTip[node] = (node.X - crackTip.X) * tangentX + (node.Y - crackTip.Y) * tangentY;
            }

            if (Logger != null) Logger.InitialLog(); //TODO: handle this with a NullLogger.
        }

        public void UpdateGeometry(double localGrowthAngle, double growthLength)
        {
            double globalGrowthAngle = AngleUtilities.Wrap(localGrowthAngle + tipSystem.RotationAngle);
            double dx = growthLength * Math.Cos(globalGrowthAngle);
            double dy = growthLength * Math.Sin(globalGrowthAngle);
            double unitDx = dx / growthLength;
            double unitDy = dy / growthLength;

            var oldTip = crackTip;
            var newTip = new CartesianPoint2D(oldTip.X + dx, oldTip.Y + dy);
            crackTip = newTip;
            tipSystem = new TipCoordinateSystem(newTip, globalGrowthAngle);

            var newSegment = new DirectedSegment2D(oldTip, newTip);
            foreach (XNode2D node in Mesh.Vertices)
            {
                // Rotate the ALL tip level sets towards the new tip and then advance them
                double rotatedTipLevelSet = (node.X - crackTip.X) * unitDx + (node.Y - crackTip.Y) * unitDy;
                levelSetsTip[node] = rotatedTipLevelSet - newSegment.Length;

                if (rotatedTipLevelSet > 0.0) // Only some body level sets are updated (See Stolarska 2001) 
                {
                    levelSetsBody[node] = newSegment.SignedDistanceOf(node);
                }
            }

            if (Logger != null) Logger.Log(); //TODO: handle this with a NullLogger.
        }

        /// <summary>
        /// Warning: with narrow band this should throw an exception if the node is not tracked.
        /// </summary>
        /// <param name="node"></param>
        /// <returns></returns>
        public double SignedDistanceOf(XNode2D node)
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
        public double SignedDistanceOf(INaturalPoint2D point, XContinuumElement2D element,
             EvaluatedInterpolation2D interpolation)
        {
            double signedDistance = 0.0;
            foreach (XNode2D node in element.Nodes)
            {
                signedDistance += interpolation.GetValueOf(node) * levelSetsBody[node];
            }
            return signedDistance;
        }

        public Tuple<double, double> SignedDistanceGradientThrough(INaturalPoint2D point,
            IReadOnlyList<XNode2D> elementNodes, EvaluatedInterpolation2D interpolation)
        {
            double gradientX = 0.0;
            double gradientY = 0.0;
            foreach (XNode2D node in elementNodes)
            {
                double levelSet = levelSetsBody[node];
                var shapeFunctionGradient = interpolation.GetGlobalCartesianDerivativesOf(node);
                gradientX += shapeFunctionGradient[0] * levelSet;
                gradientY += shapeFunctionGradient[1] * levelSet;
            }
            return new Tuple<double, double>(gradientX, gradientY);
        }

        public SortedSet<ICartesianPoint2D> FindTriangleVertices(XContinuumElement2D element)
        {
            var triangleVertices = new SortedSet<ICartesianPoint2D>(element.Nodes, pointComparer);
            int nodesCount = element.Nodes.Count;
            ElementEnrichmentType type = CharacterizeElementEnrichment(element);
            if (type != ElementEnrichmentType.Standard)
            {
                // Find the intersections between element edges and the crack. TODO: See Serafeim's Msc Thesis for a correct procedure.
                for (int i = 0; i < nodesCount; ++i)
                {
                    XNode2D node1 = element.Nodes[i];
                    XNode2D node2 = element.Nodes[(i + 1) % nodesCount];
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
                        triangleVertices.Add(new CartesianPoint2D(x, y));
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
            // The order these are called is important, as they mess with state stored in the fields
            ClearPreviousEnrichments();
            FindBodyAndTipNodesAndElements();
            ApplyFixedEnrichmentArea(tipElements[0]);
            ResolveHeavisideEnrichmentDependencies();
            ApplyEnrichmentFunctions();
        }

        private void ApplyEnrichmentFunctions() 
        {
            foreach (var node in CrackTipNodesNew)
            {
                double[] enrichmentValues = CrackTipEnrichments.EvaluateFunctionsAt(node);
                node.EnrichmentItems[CrackTipEnrichments] = enrichmentValues;
            }

            // There is no need to process each mesh node. Once a node is enriched with Heaviside it will stay that way until the 
            // end. Even if the crack curves towards itself and a crack tip comes near the node, the original discontinuity must 
            // be represented by the original Heaviside (this case creates a lot of problems and cannot be modeled with LSM 
            // accurately anyway). 
            // I am not sure if the value of the Heaviside enrichment doesn't change for elements around the crack tip, once the 
            // crack propagates. In first order LSM this is improbable, since there cannot be kinks inside the element, but what 
            // about explicit cracks and higher order LSM?
            //TODO: It could be sped up by only updating the Heaviside enrichments of nodes that have updated body  
            //      level sets, which requires tracking them.
            //      - Done. Tracking newly enriched nodes is useful for many reasons.
            //TODO: should I also clear and reapply all Heaviside enrichments? It is safer and might be useful for e.g. 
            //      reanalysis. Certainly I must not clear all node enrichments, as they may include material interfaces etc.
            //      - Ans: nope in reanalysis old Heaviside enrichments are assumed to stay the same (not sure if that is correct 
            //          though)
            foreach (var node in CrackBodyNodesNew)
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
        /// <param name="tipElement"></param>
        private void ApplyFixedEnrichmentArea(XContinuumElement2D tipElement)
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
                            CrackTipNodesNew.Add(node);
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

            foreach (XNode2D node in element.Nodes)
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
                if (IsTipElement(element, minTipLevelSet, maxTipLevelSet)) return ElementEnrichmentType.Tip;
                else if (maxTipLevelSet < 0) return ElementEnrichmentType.Heaviside;
            }
            return ElementEnrichmentType.Standard;
        }

        private void ClearPreviousEnrichments()
        {
            CrackBodyNodesNew = new HashSet<XNode2D>();
            CrackTipNodesOld = CrackTipNodesNew;
            CrackTipNodesNew = new HashSet<XNode2D>();
            foreach (var node in CrackTipNodesOld) node.EnrichmentItems.Remove(CrackTipEnrichments);
            tipElements.Clear();
        }

        private void FindBodyAndTipNodesAndElements()
        {
            foreach (var element in Mesh.Cells)
            {
                element.EnrichmentItems.Clear(); //TODO: not too fond of saving enrichment state in elements. It should be confined in nodes
                ElementEnrichmentType type = CharacterizeElementEnrichment(element);
                if (type == ElementEnrichmentType.Tip)
                {
                    tipElements.Add(element);
                    foreach (var node in element.Nodes) CrackTipNodesNew.Add(node);
                    element.EnrichmentItems.Add(CrackTipEnrichments);
                }
                else if (type == ElementEnrichmentType.Heaviside)
                {
                    foreach (var node in element.Nodes)
                    {
                        bool isNew = CrackBodyNodesAll.Add(node);
                        if (isNew) CrackBodyNodesNew.Add(node);
                    }
                    // Cut elements next to tip elements will be enriched with both Heaviside and tip functions. If all 
                    // nodes were enriched with tip functions, the element would not be enriched with Heaviside, but then it 
                    // would be a tip element and not fall under this case.
                    element.EnrichmentItems.Add(CrackBodyEnrichment);
                }
            }
            foreach (var node in CrackTipNodesNew) // tip element's nodes are not enriched with Heaviside
            {
                CrackBodyNodesAll.Remove(node);
                CrackBodyNodesNew.Remove(node);
            }

            ReportTipElements(tipElements);
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

                // The sign of the area can be derived from any node with body level set != 0
                int sign = 0;
                foreach (var vertex in triangle.Vertices)
                {
                    if (vertex is XNode2D)
                    {
                        sign = Math.Sign(levelSetsBody[(XNode2D)vertex]);
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
                    sign = Math.Sign(SignedDistanceOf(centroidNatural, element, centroidInterpolation));
                }

                if (sign > 0) positiveArea += area;
                else if (sign < 0) negativeArea += area;
                else throw new Exception(
                    "Even after finding the signed distance of its centroid, the sign of the area is unidentified");
            }
        }

        private bool IsTipElement(XContinuumElement2D elements, double minTipLevelSet, double maxTipLevelSet)
        {
            if (tipDetection == TipDetectionScheme.Implicit) //Stolarska criterion
            {
                if (minTipLevelSet * maxTipLevelSet <= 0) return true;
                else return false;
            }
            else if (tipDetection == TipDetectionScheme.Explicit)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(elements.Nodes);
                if (outline.IsPointInsidePolygon(crackTip)) return true;
                else return false;
            }
            else if (tipDetection == TipDetectionScheme.Hybrid)
            {
                if (minTipLevelSet * maxTipLevelSet <= 0)
                {
                    var outline = ConvexPolygon2D.CreateUnsafe(elements.Nodes);
                    if (outline.IsPointInsidePolygon(crackTip)) return true;
                }
                return false;
            }
            else throw new Exception("Unreachable code");
        }

        //TODO: I suspect this is method is responsible for singular global matrices popping up randomly for many configurations.
        //      Either this doesn't work as good as advertised or I messed it up. Try checking if there is at least  
        //      least 1 GP on either side of the crack instead. Bonus points if the GPs are cached somehow (e.g. all std elements
        //      have the same GPs in natural system, many enriched elements may also have the same active integration rule) for
        //      when the integration actually happens.
        private void ResolveHeavisideEnrichmentDependencies()
        {
            const double toleranceHeavisideEnrichmentArea = 1e-4;
            var processedElements = new Dictionary<XContinuumElement2D, Tuple<double, double>>();
            var nodesToRemove = new List<XNode2D>(); // Can't remove them while iterating the collection.
            foreach (var node in CrackBodyNodesAll) //TODO: Is it safe to only search the newly enriched nodes?
            {
                double nodePositiveArea = 0.0;
                double nodeNegativeArea = 0.0;

                foreach (var element in Mesh.FindElementsWithNode(node))
                {
                    bool alreadyProcessed = processedElements.TryGetValue(element, out Tuple<double, double> elementPosNegAreas);
                    if (!alreadyProcessed)
                    {
                        FindSignedAreasOfElement(element, out double elementPosArea, out double elementNegArea);
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

            foreach (var node in nodesToRemove) // using sets and set operations might be better and faster
            {
                //Console.WriteLine("Removing Heaviside enrichment from node: " + node);
                CrackBodyNodesAll.Remove(node);
                CrackBodyNodesNew.Remove(node);
            }
        }

        [ConditionalAttribute("DEBUG")]
        private void ReportTipElements(IReadOnlyList<XContinuumElement2D> tipElements)
        {
            if (!reports) return;
            Console.WriteLine("------ DEBUG/ ------");
            Console.WriteLine("Crack tip: " + crackTip);
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

        /// <summary>
        /// Represents the type of enrichment that will be applied to all nodes of the element. In LSM with linear 
        /// interpolations, an element enriched with tip functions does not need to be enriched with Heaviside too. 
        /// This is because, even if there are kinks inside the element, the linear interpolation cannot reproduce them.
        /// </summary>
        private enum ElementEnrichmentType { Standard, Heaviside, Tip }

        private enum TipDetectionScheme { Implicit, Explicit, Hybrid }
    }
}
