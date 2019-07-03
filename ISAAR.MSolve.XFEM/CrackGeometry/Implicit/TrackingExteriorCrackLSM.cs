using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Commons;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.LevelSetUpdating;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: replace BasicCrackLSM with this one after testing is done
//TODO: Consider removing the bookkeeping of enrichment items in elements. It creates a lot of opportunities for mistakes.
//      Could the enrichment type of an element be infered by just looking at its nodes, without storing state? Should it be
//      cached for efficiency and how would the cache be updated safely?
//TODO: using narrow band, a superset of tipNodes and bodyNodes would include those sets. Need core class SuperSet.
//TODO: perhaps the bookkeeping of nodes and elements can be done by a dedicated class. Narrow banding would then be implemented
//      there. In general, this is a god class and should be broken down to smaller ones.
//TODO: Crack tips should be handled differently than using enums. Interior and exterior cracks should compose their common 
//      dedicated strategy classes with their common functionality and expose appropriate properties for the crack tip data.
//TODO: Perhaps all loggers can be grouped and called together
//TODO: a lot of the tracking is just wasted memory and cpu time for most cases. It would be better to use observers to do it.
//      However, syncing the observers with the LSM is error prone and needs to be done at well defined points, without changing
//      the LSM itself and without too much memory duplication.
//TODO: A lot of functionality should be delegated to strategy classes. This can be done by having the strategy classes share
//      the fields of LSM and then mutate them when called upon. Each strategy classes gets injected with the fields it needs  
//      during construction. Alternatively I could have a readonly LSM interface that only exposes readonly properties, and a 
//      mutable one for the various strategy classes to mutate LSM data that they pull.
//TODO: If I do delegate a lot of functionality to strategy classes, how can the observers be updated correctly and efficiently,
//      namely without a lot of memory copying?
//TODO: Use a builder. It deserves one.
//TODO: a lot of bookkeeping could be done by an abstract base class.
namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit
{
    /// <summary>
    /// Warning: may misclassify elements as tip elements, causing gross errors.
    /// </summary>
    public class TrackingExteriorCrackLSM : IExteriorCrack
    {
        private static readonly bool reports = false;
        private static readonly IComparer<CartesianPoint> pointComparer = new Point2DComparerXMajor();
        private readonly Dictionary<XNode, double> levelSetsBody;
        private readonly Dictionary<XNode, double> levelSetsTip;
        private readonly double tipEnrichmentAreaRadius;
        private readonly List<XContinuumElement2D> tipElements; // Ideally there is only 1, but what if the tip falls on the edge bewteen elements?
        private readonly ILevelSetUpdater levelSetUpdater;
        private readonly IMeshInteraction meshInteraction;
        private readonly IPropagator propagator;

        private readonly List<CartesianPoint> crackPath;
        private CartesianPoint crackMouth;
        private CartesianPoint crackTip;
        private TipCoordinateSystem tipSystem;
        private ISet<XNode> crackBodyNodesAll; // TODO: a TreeSet might be better if set intersections are applied
        private ISet<XNode> crackBodyNodesModified;
        private ISet<XNode> crackBodyNodesNearModified;
        private ISet<XNode> crackBodyNodesNew;
        private ISet<XNode> crackBodyNodesRejected;
        private ISet<XNode> crackTipNodesNew;
        private ISet<XNode> crackTipNodesOld;

        public TrackingExteriorCrackLSM(IPropagator propagator, double tipEnrichmentAreaRadius,
            IHeavisideSingularityResolver singularityResolver)
        {
            this.propagator = propagator;
            this.tipEnrichmentAreaRadius = tipEnrichmentAreaRadius;

            this.crackPath = new List<CartesianPoint>();
            this.levelSetsBody = new Dictionary<XNode, double>();
            this.levelSetsTip = new Dictionary<XNode, double>();
            this.tipElements = new List<XContinuumElement2D>();

            this.crackBodyNodesAll = new HashSet<XNode>();
            this.crackBodyNodesModified = new HashSet<XNode>();
            this.crackBodyNodesNearModified = new HashSet<XNode>();
            this.crackBodyNodesNew = new HashSet<XNode>();
            this.crackBodyNodesRejected = new HashSet<XNode>();
            this.crackTipNodesNew = new HashSet<XNode>();
            this.crackTipNodesOld = new HashSet<XNode>();
            this.ElementsModified = new HashSet<XContinuumElement2D>();

            //this.levelSetUpdater = new LevelSetUpdaterOLD();
            this.levelSetUpdater = new LevelSetUpdaterStolarska();
            //this.meshInteraction = new StolarskaMeshInteraction(this);
            //this.meshInteraction = new HybridMeshInteraction(this);
            this.meshInteraction = new SerafeimMeshInteraction(this);
            this.SingularityResolver = singularityResolver;
        }

        public TrackingExteriorCrackLSM(IPropagator propagator, double tipEnrichmentAreaRadius = 0.0) :
            this(propagator, tipEnrichmentAreaRadius, null)
        {
            this.SingularityResolver = new RelativeAreaResolver(1E-4);
        }

        // TODO: Not too fond of the setters, but at least the enrichments are immutable. Perhaps I can pass their
        // parameters to a CrackDescription builder and construct them there, without involving the user.
        public CrackBodyEnrichment2D CrackBodyEnrichment { get; set; }
        public CrackTipEnrichments2D CrackTipEnrichments { get; set; }
        public IReadOnlyList<CartesianPoint> CrackPath { get { return crackPath; } }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesAll
        {
            get
            {
                return new Dictionary<CrackBodyEnrichment2D, ISet<XNode>> { [CrackBodyEnrichment] = crackBodyNodesAll };
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesModified
        {
            get
            {
                return new Dictionary<CrackBodyEnrichment2D, ISet<XNode>> { [CrackBodyEnrichment] = crackBodyNodesModified };
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNearModified
        {
            get
            {
                return new Dictionary<CrackBodyEnrichment2D, ISet<XNode>> { [CrackBodyEnrichment] = crackBodyNodesNearModified };
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNew
        {
            get
            {
                return new Dictionary<CrackBodyEnrichment2D, ISet<XNode>> { [CrackBodyEnrichment] = crackBodyNodesNew };
            }
        }

        public IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesRejected
        {
            get
            {
                return new Dictionary<CrackBodyEnrichment2D, ISet<XNode>> { [CrackBodyEnrichment] = crackBodyNodesRejected };
            }
        }

        public IReadOnlyList<CartesianPoint> CrackTips { get { return new CartesianPoint[] { crackTip }; } }

        public IReadOnlyDictionary<CartesianPoint, IReadOnlyList<XContinuumElement2D>> CrackTipElements
        {
            get
            {
                var crackTipElements = new Dictionary<CartesianPoint, IReadOnlyList<XContinuumElement2D>>();
                crackTipElements.Add(this.crackTip, this.tipElements);
                return crackTipElements;
            }
        }

        public IReadOnlyDictionary<CartesianPoint, IPropagator> CrackTipPropagators
        {
            get
            {
                var tipPropagators = new Dictionary<CartesianPoint, IPropagator>();
                tipPropagators.Add(this.crackTip, this.propagator);
                return tipPropagators;
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesNew
        {
            get
            {
                return new Dictionary<CrackTipEnrichments2D, ISet<XNode>> { [CrackTipEnrichments] = crackTipNodesNew };
            }
        }

        public IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesOld
        {
            get
            {
                return new Dictionary<CrackTipEnrichments2D, ISet<XNode>> { [CrackTipEnrichments] = crackTipNodesOld };
            }
        }

        public ISet<XContinuumElement2D> ElementsModified { get; private set; }

        public IReadOnlyList<IEnrichmentItem2D> Enrichments
        {
            get { return new IEnrichmentItem2D[] { CrackBodyEnrichment, CrackTipEnrichments }; }
        }

        public Dictionary<XNode, double> LevelSetsBody { get { return levelSetsBody; } }
        public Dictionary<XNode, double> LevelSetsTip { get { return levelSetsTip; } }
        public EnrichmentLogger EnrichmentLogger { get; set; }
        public LevelSetLogger LevelSetLogger { get; set; }
        public PreviousLevelSetComparer LevelSetComparer { get; set; }
        public BidirectionalMesh2D<XNode, XContinuumElement2D> Mesh { get; set; }
        public IHeavisideSingularityResolver SingularityResolver { get; }
        public IReadOnlyList<ISingleCrack> SingleCracks { get { return new ISingleCrack[] { this }; } }

        // TODO: remove this
        public CartesianPoint GetCrackTip(CrackTipPosition tipPosition) 
        {
            if (tipPosition == CrackTipPosition.Single) return crackTip;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        // TODO: remove this
        public TipCoordinateSystem GetTipSystem(CrackTipPosition tip)
        {
            if (tip == CrackTipPosition.Single) return tipSystem;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        // TODO: remove this
        public IReadOnlyList<XContinuumElement2D> GetTipElements(CrackTipPosition tipPosition)
        {
            if (tipPosition == CrackTipPosition.Single) return tipElements;
            else throw new ArgumentException("Only works for single tip cracks.");
        }

        public void InitializeGeometry(CartesianPoint crackMouth, CartesianPoint crackTip)
        {
            this.crackMouth = crackMouth;
            crackPath.Add(crackMouth);
            var segment = new DirectedSegment2D(crackMouth, crackTip);

            double tangentX = crackTip.X - crackMouth.X;
            double tangentY = crackTip.Y - crackMouth.Y;
            double length = Math.Sqrt(tangentX * tangentX + tangentY * tangentY);
            double tangentSlope = Math.Atan2(tangentY, tangentX);
            this.crackTip = crackTip;
            crackPath.Add(crackTip);
            tipSystem = new TipCoordinateSystem(crackTip, tangentSlope);
            CrackTipEnrichments.TipSystem = tipSystem;

            tangentX /= length;
            tangentY /= length;

            foreach (XNode node in Mesh.Nodes)
            {
                levelSetsBody[node] = segment.SignedDistanceOf(node);
                levelSetsTip[node] = (node.X - crackTip.X) * tangentX + (node.Y - crackTip.Y) * tangentY;
            }

            if (LevelSetLogger != null) LevelSetLogger.InitialLog(); //TODO: handle this with a NullLogger.
        }

        //TODO: This should work for any IOpenCurve2D. Same for all ICrackGeometryDescriptions.
        //TODO: The tangent stuff should be done by the initial curve.
        public void InitializeGeometry(PolyLine2D initialCrack)
        {
            foreach (var vertex in initialCrack.Vertices) crackPath.Add(vertex);

            crackMouth = initialCrack.Start;
            var lastSegment = initialCrack.Segments[initialCrack.Segments.Count - 1];
            
            double tangentX = lastSegment.End.X - lastSegment.Start.X;
            double tangentY = lastSegment.End.Y - lastSegment.Start.Y;
            double length = Math.Sqrt(tangentX * tangentX + tangentY * tangentY);
            double tangentSlope = Math.Atan2(tangentY, tangentX);
            this.crackTip = initialCrack.End;
            tipSystem = new TipCoordinateSystem(crackTip, tangentSlope);
            CrackTipEnrichments.TipSystem = tipSystem;

            tangentX /= length;
            tangentY /= length;

            foreach (XNode node in Mesh.Nodes)
            {
                levelSetsBody[node] = initialCrack.SignedDistanceOf(node);
                levelSetsTip[node] = (node.X - crackTip.X) * tangentX + (node.Y - crackTip.Y) * tangentY;
            }

            if (LevelSetLogger != null) LevelSetLogger.InitialLog(); //TODO: handle this with a NullLogger.
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
                double dNdx = interpolation.ShapeGradientsCartesian[nodeIdx, 0];
                double dNdy = interpolation.ShapeGradientsCartesian[nodeIdx, 1];

                double levelSet = levelSetsBody[elementNodes[nodeIdx]];
                gradientX += dNdx * levelSet;
                gradientY += dNdy * levelSet;
            }
            return new Tuple<double, double>(gradientX, gradientY);
        }

        public SortedSet<CartesianPoint> FindTriangleVertices(XContinuumElement2D element)
        {
            var triangleVertices = new SortedSet<CartesianPoint>(element.Nodes, pointComparer);
            int nodesCount = element.Nodes.Count;
            CrackElementPosition relativePosition = meshInteraction.FindRelativePositionOf(element);
            if (relativePosition != CrackElementPosition.Irrelevant)
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

                if (relativePosition == CrackElementPosition.ContainsTip) triangleVertices.Add(crackTip);
            }
            return triangleVertices;
        }

        public void Propagate(Dictionary<int, Vector> totalFreeDisplacements)
        {
            (double growthAngle, double growthLength) = propagator.Propagate(totalFreeDisplacements,
                crackTip, tipSystem, tipElements);
            UpdateGeometry(growthAngle, growthLength);
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

            // Modified elements
            ElementsModified.Clear();
            var modifiedNodes = new HashSet<XNode>(crackBodyNodesNew);
            modifiedNodes.UnionWith(crackTipNodesNew);
            modifiedNodes.UnionWith(crackTipNodesOld);
            modifiedNodes.UnionWith(crackBodyNodesModified);
            foreach (var node in modifiedNodes)
            {
                foreach (var element in Mesh.FindElementsWithNode(node)) ElementsModified.Add(element);
            }

            // Unmodified body nodes of modified elements
            crackBodyNodesNearModified.Clear();
            foreach (var element in ElementsModified)
            {
                foreach (var node in element.Nodes)
                {
                    // Only Heaviside enriched nodes of these elements are of concern.
                    // TODO: what about tip enriched nodes?
                    bool isEnrichedHeaviside = node.EnrichmentItems.ContainsKey(CrackBodyEnrichment);
                    if (!modifiedNodes.Contains(node) && isEnrichedHeaviside) crackBodyNodesNearModified.Add(node);

                    //WARNING: The next would also mark std nodes which is incorrect. Code left over as a reminder.
                    //if (!modifiedNodes.Contains(node)) CrackBodyNodesNearModified.Add(node);
                }
            }

            if (EnrichmentLogger != null) EnrichmentLogger.Log(); //TODO: handle this with a NullLogger.
            if (LevelSetComparer != null) LevelSetComparer.Log();
        }

        //TODO: make this private
        public void UpdateGeometry(double localGrowthAngle, double growthLength)
        {
            double globalGrowthAngle = MathUtilities.WrapAngle(localGrowthAngle + tipSystem.RotationAngle);
            double dx = growthLength * Math.Cos(globalGrowthAngle);
            double dy = growthLength * Math.Sin(globalGrowthAngle);
            var oldTip = crackTip;
            var newTip = new CartesianPoint(oldTip.X + dx, oldTip.Y + dy);
            crackTip = newTip;
            crackPath.Add(newTip);
            tipSystem = new TipCoordinateSystem(newTip, globalGrowthAngle);
            CrackTipEnrichments.TipSystem = tipSystem;

            //TODO: it is inconsistent that the modified body nodes are updated here, while the other in UpdateEnrichments(); 
            crackBodyNodesModified = levelSetUpdater.Update(oldTip, localGrowthAngle, growthLength, dx, dy, Mesh.Nodes, 
                crackBodyNodesAll, levelSetsBody, levelSetsTip);

            if (LevelSetLogger != null) LevelSetLogger.Log(); //TODO: handle this with a NullLogger.
        }

        private void ApplyEnrichmentFunctions() 
        {
            foreach (var node in crackTipNodesNew)
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
            foreach (var node in crackBodyNodesNew)
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
                            crackTipNodesNew.Add(node);
                        }
                        else completelyInside = false;
                    }

                    //TODO: I tried an alternative approach, ie elements access their enrichments from their nodes. 
                    //      My original thought that this approach (storing enrichments in elements, unless they are standard /
                    //      blending) wouldn't work for blending elements, was incorrect, as elements with 0 enrichments
                    //      were then examined and separated into standard / blending.
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

        private void ClearPreviousEnrichments()
        {
            //TODO: this method should not exist. Clearing these sets should be done in the methods where they are updated.
            //WARNING: Do not clear CrackBodyNodesModified here, since they are updated during UpdateGeometry() which is called first. 
            crackBodyNodesNew = new HashSet<XNode>();
            crackTipNodesOld = crackTipNodesNew;
            crackTipNodesNew = new HashSet<XNode>();
            foreach (var node in crackTipNodesOld) node.EnrichmentItems.Remove(CrackTipEnrichments);
            tipElements.Clear();
        }

        private void FindBodyAndTipNodesAndElements() //TODO: perhaps the singularity resolution should be done inside this method
        {
            foreach (var element in Mesh.Elements)
            {
                //TODO: I tried an alternative approach, ie elements access their enrichments from their nodes. 
                //      My original thought that this approach (storing enrichments in elements, unless they are standard /
                //      blending) wouldn't work for blending elements, was incorrect, as elements with 0 enrichments
                //      were then examined and separated into standard / blending.
                element.EnrichmentItems.Clear(); //TODO: not too fond of saving enrichment state in elements. It should be confined in nodes

                CrackElementPosition relativePosition = meshInteraction.FindRelativePositionOf(element);
                if (relativePosition == CrackElementPosition.ContainsTip)
                {
                    tipElements.Add(element);
                    foreach (var node in element.Nodes) crackTipNodesNew.Add(node);

                    element.EnrichmentItems.Add(CrackTipEnrichments);
                }
                else if (relativePosition == CrackElementPosition.Intersected)
                {
                    foreach (var node in element.Nodes)
                    {
                        bool isNew = crackBodyNodesAll.Add(node);
                        if (isNew) crackBodyNodesNew.Add(node);
                    }

                    // Cut elements next to tip elements will be enriched with both Heaviside and tip functions. If all 
                    // nodes were enriched with tip functions, the element would not be enriched with Heaviside, but then it 
                    // would be a tip element and not fall under this case.
                    element.EnrichmentItems.Add(CrackBodyEnrichment);
                }
            }
            foreach (var node in crackTipNodesNew) // tip element's nodes are not enriched with Heaviside
            {
                crackBodyNodesAll.Remove(node);
                crackBodyNodesNew.Remove(node);
            }

            ReportTipElements(tipElements);
        }

        private void ResolveHeavisideEnrichmentDependencies()
        {
            //TODO: Is it safe to only search the newly enriched nodes? Update the rejected set appropriately
            crackBodyNodesRejected.Clear();
            crackBodyNodesRejected = SingularityResolver.FindHeavisideNodesToRemove(this, Mesh, crackBodyNodesAll); 
            foreach (var node in crackBodyNodesRejected) // using set operations might be better and faster
            {
                //Console.WriteLine("Removing Heaviside enrichment from node: " + node);
                crackBodyNodesAll.Remove(node);
                crackBodyNodesNew.Remove(node);
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
    }
}
