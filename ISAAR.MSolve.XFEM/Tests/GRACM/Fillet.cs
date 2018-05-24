using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit.Logging;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;


namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class Fillet: IBenchmark
    {
        public static Builder SetupBenchmark()
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Meshes\fillet.msh";
            string plotPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Plots";
            string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Timing";
            //string meshPath = @"C:\Users\seraf\Desktop\GRACM\Fillet\Meshes\fillet.msh";
            //string plotPath = @"C:\Users\seraf\Desktop\GRACM\Fillet\Plots";
            //string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Fillet\Timing";

            double growthLength = 6; // mm. Must be sufficiently larger than the element size.
            var builder = new Builder(meshPath, growthLength, timingPath);
            builder.LsmOutputDirectory = plotPath;
            builder.MaxIterations = 20;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0; 

            // TODO: enter the fixed propagator here, perhaps by solving the benchmark once.
            // These are for mesh: ...
            //var propagationLogger = new PropagationLogger();
            //propagationLogger.GrowthAngles.Add(-0.0349434289780521);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.0729848767244629);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.125892740180586);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.200116860828933);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.258299391791769);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.264803465603906);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.201411670680886);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.123234163953279);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //propagationLogger.GrowthAngles.Add(-0.0816346096256186);
            //propagationLogger.GrowthLengths.Add(growthLength);
            //builder.KnownPropagation = propagationLogger;

            return builder;
        }


        #region constants
        /// <summary>
        /// Thickness of the whole domain
        /// </summary>
        private const double t = 1.0; // mm

        /// <summary>
        /// Young's modulus
        /// </summary>
        private const double E = 2.1e12; // kN/mm^2

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3;

        /// <summary>
        /// The material used for the J-integral computation. It msut be stored separately from individual element materials.
        /// </summary>
        private static readonly HomogeneousElasticMaterial2D globalHomogeneousMaterial =
            HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);

        /// <summary>
        /// The maximum value that the effective SIF can reach before collapse occurs.
        /// </summary>
        private const double fractureToughness = double.MaxValue;

        /// <summary>
        /// Magnitude in kN of "opening" load
        /// </summary>
        private const double load = 20.0; // kN


        // Geometry (m)
        private const double bottomWidth = 375.0, topWidth = 75.0, radius = 20.0;// mm
        private const double flangeHeight = 75.0, totalHeight = 150.0; // mm
        private const double crackHeight = flangeHeight + radius, crackLength = 5, webLeft = 0.5 * (bottomWidth - topWidth); //mm
        private const double infCrackHeight = 90.0, supCrackHeight = 105.0; //mm

        #endregion
        private readonly bool constrainBottomEdge;

        /// <summary>
        /// The length by which the crack grows in each iteration.
        /// </summary>
        private readonly double growthLength;

        /// <summary>
        /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
        /// radius is larger than the length of the crack segments.
        /// </summary>
        private readonly double jIntegralRadiusOverElementSize;

        /// <summary>
        /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
        /// rather than the using the J-integral method to predict it.
        /// </summary>
        private readonly PropagationLogger knownPropagation;

        private readonly string meshPath;
        private readonly string lsmOutputDirectory;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;

        private BiMesh2D mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private Fillet(string meshPath, double growthLength, double jIntegralRadiusOverElementSize,
             PropagationLogger knownPropagation, string lsmOutputDirectory, int maxIterations, bool constrainBottomEdge)
        {
            this.meshPath = meshPath;
            this.growthLength = growthLength;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.knownPropagation = knownPropagation;
            this.lsmOutputDirectory = lsmOutputDirectory;
            this.maxIterations = maxIterations;
            this.constrainBottomEdge = constrainBottomEdge;
        }

        /// <summary>
        /// The crack geometry description
        /// </summary>
        public TrackingExteriorCrackLSM Crack { get; private set; }

        public IDecomposer Decomposer { get; private set; }

        public IReadOnlyList<XNode2D> EnrichedArea { get; private set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public Model2D Model { get; private set; }


        public IReadOnlyList<ICartesianPoint2D> Analyze(ISolver solver)
        {
            //var crackPath = new List<ICartesianPoint2D>();
            //crackPath.Add(new CartesianPoint2D(webLeft, crackHeight));
            //crackPath.Add(new CartesianPoint2D(webLeft + crackLength, crackHeight));

            var actualPropagator = new Propagator(mesh, Crack, CrackTipPosition.Single, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            IPropagator propagator;
            if (knownPropagation != null) propagator = new FixedPropagator(knownPropagation, null);
            else propagator = actualPropagator;
            var analysis = new QuasiStaticAnalysis(Model, mesh, Crack, solver, propagator, fractureToughness, maxIterations);
            analysis.Analyze();
            //crackPath.AddRange();
            return Crack.CrackPath;
        }

        public void InitializeModel()
        {
            Model = new Model2D();
            CreateMesh();
            ApplyBoundaryConditions();
            InitializeCrack();
            LimitEnrichedArea();
            DomainDecomposition();
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(Model, 1e-6);

            // Constraints
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            Model.AddConstraint(bottomLeftNode, DisplacementDof.X, 0.0);
            if (constrainBottomEdge)
            {
                foreach (var node in finder.FindNodesWithY(0.0))
                {
                    Model.AddConstraint(node, DisplacementDof.Y, 0.0);
                }
            }
            else
            {
                XNode2D bottomRightNode = finder.FindNodeWith(bottomWidth, 0.0);
                Model.AddConstraint(bottomRightNode, DisplacementDof.Y, 0.0);
                Model.AddConstraint(bottomLeftNode, DisplacementDof.Y, 0.0);
            }


            // Loads
            double distributedLoad = load / topWidth;
            var distrubutor = new LoadDistributor();
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(totalHeight);
            double[,] topLoads = distrubutor.DistributeLoad(topNodes, 0.0, distributedLoad);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                Model.AddNodalLoad(topNodes[i], DisplacementDof.X, topLoads[i, 0]);
                Model.AddNodalLoad(topNodes[i], DisplacementDof.Y, topLoads[i, 1]);
            }
        }

        private void CreateMesh()
        {
            // Mesh generation
            var reader = new GmshReader(meshPath, true);
            Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> meshEntities = reader.ReadMesh();
            var nodes = new XNode2D[meshEntities.Item1.Count];
            var elementConnectivity = new List<XNode2D[]>();
            nodes = meshEntities.Item1.ToArray();
            foreach (var element in meshEntities.Item2)
            {
                elementConnectivity.Add(element.Nodes.ToArray());
            }

            // Nodes
            foreach (XNode2D node in nodes) Model.AddNode(node);

            // Integration rules
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            // Elements
            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
                Model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }

            // Mesh usable for crack-mesh interaction
            
            mesh = new BiMesh2D(Model.Nodes, Model.Elements, new FilletBoundary());
        }

        private void DomainDecomposition() //TODO: this should not be hardcoded, but provided by the caller of the solver
        {
            //var regions = new PolygonalRegion[4];

            //var vertices1 = new CartesianPoint2D[4];
            //vertices1[0] = new CartesianPoint2D(0.0, h);
            //vertices1[1] = new CartesianPoint2D(0.0, 0.0);
            //vertices1[2] = new CartesianPoint2D(L / 5.0, 0.0);
            //vertices1[3] = new CartesianPoint2D(L / 3.0, h);
            //var boundaries1 = new HashSet<LineSegment2D>();
            //boundaries1.Add(new LineSegment2D(vertices1[2], vertices1[3]));
            //regions[0] = new PolygonalRegion(vertices1, boundaries1);

            //var vertices2 = new CartesianPoint2D[4];
            //vertices2[0] = new CartesianPoint2D(L / 3.0, h);
            //vertices2[1] = new CartesianPoint2D(L / 5.0, 0.0);
            //vertices2[2] = new CartesianPoint2D(L / 2.0, 0.0);
            //vertices2[3] = new CartesianPoint2D(L / 2.0, h);
            //var boundaries2 = new HashSet<LineSegment2D>();
            //boundaries2.Add(new LineSegment2D(vertices2[0], vertices2[1]));
            //boundaries2.Add(new LineSegment2D(vertices2[2], vertices2[3]));
            //regions[1] = new PolygonalRegion(vertices2, boundaries2);

            //var vertices3 = new CartesianPoint2D[4];
            //vertices3[0] = new CartesianPoint2D(L / 2.0, h);
            //vertices3[1] = new CartesianPoint2D(L / 2.0, 0.0);
            //vertices3[2] = new CartesianPoint2D(3.0 * L / 4.0, 0.0);
            //vertices3[3] = new CartesianPoint2D(2.0 * L / 3.0, h);
            //var boundaries3 = new HashSet<LineSegment2D>();
            //boundaries3.Add(new LineSegment2D(vertices3[0], vertices3[1]));
            //boundaries3.Add(new LineSegment2D(vertices3[2], vertices3[3]));
            //regions[2] = new PolygonalRegion(vertices3, boundaries3);

            //var vertices4 = new CartesianPoint2D[4];
            //vertices4[0] = new CartesianPoint2D(2.0 * L / 3.0, h);
            //vertices4[1] = new CartesianPoint2D(3.0 * L / 4.0, 0.0);
            //vertices4[2] = new CartesianPoint2D(L, 0.0);
            //vertices4[3] = new CartesianPoint2D(L, h);
            //var boundaries4 = new HashSet<LineSegment2D>();
            //boundaries4.Add(new LineSegment2D(vertices4[0], vertices4[1]));
            //regions[3] = new PolygonalRegion(vertices4, boundaries4);

            //Decomposer = new GuideDecomposer(regions, mesh);
        }

        private void InitializeCrack()
        {
            var crackMouth = new CartesianPoint2D(webLeft, crackHeight);
            var crackTip = new CartesianPoint2D(webLeft + crackLength, crackHeight);
            var initialCrack = new PolyLine2D(crackMouth, crackTip);
            var lsmCrack = new TrackingExteriorCrackLSM();
            lsmCrack.Mesh = mesh;

            // Create enrichments          
            lsmCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(lsmCrack);
            lsmCrack.CrackTipEnrichments = new CrackTipEnrichments2D(lsmCrack, CrackTipPosition.Single);
            if (lsmOutputDirectory != null)
            {
                lsmCrack.EnrichmentLogger = new EnrichmentLogger(Model, lsmCrack, lsmOutputDirectory);
                lsmCrack.LevelSetLogger = new LevelSetLogger(Model, lsmCrack, lsmOutputDirectory);
                lsmCrack.LevelSetComparer = new PreviousLevelSetComparer(lsmCrack, lsmOutputDirectory);
            }

            // Mesh geometry interaction
            lsmCrack.InitializeGeometry(initialCrack);
            this.Crack = lsmCrack;
        }

        private void LimitEnrichedArea()
        {
            EnrichedArea = Model.Nodes.Where(node => (node.Y >= infCrackHeight) && (node.Y <= supCrackHeight)).ToList();
        }

        public class Builder: IBenchmarkBuilder
        {
            private readonly double growthLength;
            private readonly string meshPath;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="meshPath">The absolute path of the mesh file.</param>
            /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
            /// <param name="timingPath">The absolute path of the file where slover timing will be written.</param>
            public Builder(string meshPath, double growthLength, string timingPath)
            {
                this.growthLength = growthLength;
                this.meshPath = meshPath;
                this.TimingPath = timingPath;
            }

            /// <summary>
            /// If set to true, the global y dof will be constrained for all nodes along the bottom edge, in order to simulate 
            /// a very thick rigid I-beam. If set to false, the global � degree of freedom is constrained only for the corner
            ///  nodes of the bottom edge, in order to simulate a very thin, flexible I-beam
            /// </summary>
            public bool ConstrainBottomEdge { get; set; } = true;

            /// <summary>
            /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
            /// radius is larger than the length of the crack segments.
            /// </summary>
            public double JintegralRadiusOverElementSize { get; set; } = 2.0;

            /// <summary>
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            public PropagationLogger KnownPropagation { get; set; } = null;

            /// <summary>
            /// The absolute path of the directory where output vtk files with the crack path and the level set functions at 
            /// each iteration will be written. Leave it null to avoid the performance cost it will introduce.
            /// </summary>
            public string LsmOutputDirectory { get; set; } = null;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = int.MaxValue;

            /// <summary>
            /// The absolute path of the file where slover timing will be written.
            /// </summary>
            public string TimingPath { get; }

            public IBenchmark BuildBenchmark()
            {
                return new Fillet(meshPath, growthLength, JintegralRadiusOverElementSize, KnownPropagation, LsmOutputDirectory,
                    MaxIterations, ConstrainBottomEdge);
            }
        }

        private class FilletBoundary : IDomainBoundary
        {
            private readonly double voidRectWidth, centerY, leftCenterX, rightCenterX;

            public FilletBoundary()
            {
                voidRectWidth = 0.5 * (bottomWidth - topWidth);
                centerY = flangeHeight + radius;
                leftCenterX = voidRectWidth - radius;
                rightCenterX = bottomWidth - voidRectWidth + radius;
            }

            public bool IsInside(ICartesianPoint2D point)
            {
                // Shapes
                var rectHull = new RectangularBoundary(0.0, bottomWidth, 0.0, totalHeight);
                var leftVoid = new RectangularBoundary(0.0, voidRectWidth, flangeHeight, totalHeight);
                var rightVoid = new RectangularBoundary(bottomWidth - voidRectWidth, bottomWidth, flangeHeight, totalHeight);
                var leftCircle = new Circle2D(new CartesianPoint2D(leftCenterX, centerY), radius);
                var rightCircle = new Circle2D(new CartesianPoint2D(rightCenterX, centerY), radius);

                if (rectHull.IsInside(point))
                {
                    if (leftVoid.IsInside(point)) // Over flange, left of web
                    {
                        if ((point.X > leftCenterX) && (point.Y < centerY))
                        {
                            if (leftCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside)
                            {
                                return true; // Inside left fillet
                            }
                            else return false;
                        }
                        else return false;
                    }
                    else if (rightVoid.IsInside(point)) // Over flange, right of web
                    {
                        if ((point.X < leftCenterX) && (point.Y < centerY))
                        {
                            if (leftCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside)
                            {
                                return true; // Inside right fillet
                            }
                            else return false;
                        }
                        else return false;
                    }
                    else return true; // Inside the flange or the web
                }
                else return false;
            }
        }
    }
}
