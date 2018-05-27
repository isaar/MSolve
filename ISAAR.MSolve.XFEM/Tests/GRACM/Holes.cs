using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
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
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;

//TODO: use fixed enrichment area
//TODO: limit enriched nodes in the areas around the cracks
namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class Holes: IBenchmark
    {
        public static Builder SetupBenchmark()
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Meshes\holes.msh";
            string plotPathLeft = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Plots\Left";
            string plotPathRight = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Plots\Right";
            string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Holes\Timing\results.txt";
            //string meshPath = @"C:\Users\seraf\Desktop\GRACM\Holes\Meshes\holes.msh";
            //string plotPathLeft = @"C:\Users\seraf\Desktop\GRACM\Benchmark_Holes\Plots\Left";
            //string plotPathRight = @"C:\Users\seraf\Desktop\GRACM\Benchmark_Holes\Plots\Right";
            //string timingPath = @"C:\Users\seraf\Desktop\GRACM\Holes\Timing\results.txt";

            double growthLength = 1.0; // mm. Must be sufficiently larger than the element size. 0.8 works for left crack, 0.9 is so-so for right crack
            var builder = new Builder(meshPath, growthLength, timingPath);
            builder.HeavisideEnrichmentTolerance = 0.03;
            builder.LeftLsmPlotDirectory = plotPathLeft;
            builder.RightLsmPlotDirectory = plotPathRight;
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
        private const double E = 2e5; // N/mm^2

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3;

        /// <summary>
        /// Fracture toughness
        /// </summary>
        private const double KIc = double.MaxValue; // N.mm^(3/2) actually 1500

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
        /// Prescribed displacement
        /// </summary>
        private const double displacement = 0.1; //mm

        // Geometry
        private const double minX = 0.0, minY = 0.0, maxX = 20.0, maxY = 10.0; //mm
        private const double holeRadius = 2.0, leftHoleX = 3.0, leftHoleY = 7.0, rightHoleX = 17.0, rightHoleY = 3.0; //mm
        private const double initialCrackLength = 1.0;
        private const double leftCrackMouthX = 0.0, leftCrackMouthY = 2.85, leftCrackTipX = minX + initialCrackLength, leftCrackTipY = 2.85; //mm
        private const double rightCrackMouthX = 20.0, rightCrackMouthY = 7.15, rightCrackTipX = maxX - initialCrackLength, rightCrackTipY = 7.15; //mm

        private const bool fixCornersInsteadOfBottom = false;
        #endregion

        /// <summary>
        /// The length by which the crack grows in each iteration.
        /// </summary>
        private readonly double growthLength;

        private readonly double heavisideTol;

        /// <summary>
        /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
        /// radius is larger than the length of the crack segments.
        /// </summary>
        private readonly double jIntegralRadiusOverElementSize;

        /// <summary>
        /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
        /// rather than the using the J-integral method to predict it.
        /// </summary>
        private readonly PropagationLogger knownLeftPropagation;

        /// <summary>
        /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
        /// rather than the using the J-integral method to predict it.
        /// </summary>
        private readonly PropagationLogger knownRightPropagation;

        private readonly string meshPath;
        private readonly string leftLsmPlotDirectory;
        private readonly string rightLsmPlotDirectory;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;

        private TrackingExteriorCrackLSM leftCrack;
        private TrackingExteriorCrackLSM rightCrack;
        private BiMesh2D mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private Holes(string meshPath, double growthLength, double jIntegralRadiusOverElementSize,
             PropagationLogger knownLeftPropagation, PropagationLogger knownRightPropagation, 
             string leftLsmPlotDirectory, string rightLsmPlotDirectory, int maxIterations, double heavisideTol)
        {
            this.meshPath = meshPath;
            this.growthLength = growthLength;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.knownLeftPropagation = knownLeftPropagation;
            this.knownRightPropagation = knownRightPropagation;
            this.leftLsmPlotDirectory = leftLsmPlotDirectory;
            this.rightLsmPlotDirectory = rightLsmPlotDirectory;
            this.maxIterations = maxIterations;
            this.heavisideTol = heavisideTol;
        }

        /// <summary>
        /// The crack geometry description
        /// </summary>
        public ICrackDescription Crack { get; private set; }

        public IDecomposer Decomposer { get; private set; }

        public IReadOnlyList<XNode2D> EnrichedArea { get; private set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public Model2D Model { get; private set; }

        public Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> PossibleEnrichments { get; private set; }

        public void Analyze(ISolver solver)
        {
            var analysis = new QuasiStaticAnalysis(Model, mesh, Crack, solver, fractureToughness, maxIterations);
            analysis.Analyze();

            Console.WriteLine("Left crack path:");
            foreach (var point in leftCrack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();
            Console.WriteLine("Left crack growth angles:");
            foreach (var angle in leftCrack.GetCrackTipPropagators()[0].Logger.GrowthAngles)
            {
                Console.WriteLine(angle);
            }
            Console.WriteLine();
            Console.WriteLine("Right crack path:");
            foreach (var point in rightCrack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();
            Console.WriteLine("Right crack growth angles:");
            foreach (var angle in rightCrack.GetCrackTipPropagators()[0].Logger.GrowthAngles)
            {
                Console.WriteLine(angle);
            }
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

            // Dirichlet
            foreach (var node in finder.FindNodesWithY(minY))
            {
                Model.AddConstraint(node, DisplacementDof.X, 0.0);
                Model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }
            foreach (var node in finder.FindNodesWithY(maxY))
            {
                Model.AddConstraint(node, DisplacementDof.Y, displacement);
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

            mesh = new BiMesh2D(Model.Nodes, Model.Elements, new HolesBoundary());
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
            // Left crack
            // Propagator
            IPropagator leftPropagator;
            var actualPropagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            if (knownLeftPropagation != null) leftPropagator = new FixedPropagator(knownLeftPropagation, null);
            else leftPropagator = actualPropagator;

            var initialLeftCrack = new PolyLine2D(new CartesianPoint2D(leftCrackMouthX, leftCrackMouthY),
                new CartesianPoint2D(leftCrackTipX, leftCrackTipY));
            leftCrack = new TrackingExteriorCrackLSM(leftPropagator, 0.0, new RelativeAreaResolver(heavisideTol));
            leftCrack.Mesh = mesh;

            // Create enrichments          
            leftCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(leftCrack);
            leftCrack.CrackTipEnrichments = new CrackTipEnrichments2D(leftCrack, CrackTipPosition.Single);
            if (leftLsmPlotDirectory != null)
            {
                leftCrack.EnrichmentLogger = new EnrichmentLogger(Model, leftCrack, leftLsmPlotDirectory);
                leftCrack.LevelSetLogger = new LevelSetLogger(Model, leftCrack, leftLsmPlotDirectory);
                leftCrack.LevelSetComparer = new PreviousLevelSetComparer(leftCrack, leftLsmPlotDirectory);
            }

            // Mesh geometry interaction
            leftCrack.InitializeGeometry(initialLeftCrack);


            // Right crack
            // Propagator
            IPropagator rightPropagator;
            actualPropagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            if (knownLeftPropagation != null) rightPropagator = new FixedPropagator(knownRightPropagation, null);
            else rightPropagator = actualPropagator;

            var initialRightCrack = new PolyLine2D(new CartesianPoint2D(rightCrackMouthX, rightCrackMouthY),
                new CartesianPoint2D(rightCrackTipX, rightCrackTipY));
            rightCrack = new TrackingExteriorCrackLSM(rightPropagator, 0.0, new RelativeAreaResolver(heavisideTol));
            rightCrack.Mesh = mesh;

            // Create enrichments          
            rightCrack.CrackBodyEnrichment = new CrackBodyEnrichment2D(rightCrack);
            rightCrack.CrackTipEnrichments = new CrackTipEnrichments2D(rightCrack, CrackTipPosition.Single);
            if (rightLsmPlotDirectory != null)
            {
                rightCrack.EnrichmentLogger = new EnrichmentLogger(Model, rightCrack, rightLsmPlotDirectory);
                rightCrack.LevelSetLogger = new LevelSetLogger(Model, rightCrack, rightLsmPlotDirectory);
                rightCrack.LevelSetComparer = new PreviousLevelSetComparer(rightCrack, rightLsmPlotDirectory);
            }

            // Mesh geometry interaction
            rightCrack.InitializeGeometry(initialRightCrack);


            // Container for both cracks
            Crack = new MultipleCracksDisjoint(new TrackingExteriorCrackLSM[] { leftCrack, rightCrack });
        }

        private void LimitEnrichedArea()
        {
            //TODO: do that
            //EnrichedArea = Model.Nodes.Where(node => (node.Y >= infCrackHeight) && (node.Y <= supCrackHeight)).ToList();
            PossibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();

            var leftEnrichedNodes = Model.Nodes;
            PossibleEnrichments.Add(leftCrack.CrackBodyEnrichment, leftEnrichedNodes);
            PossibleEnrichments.Add(leftCrack.CrackTipEnrichments, leftEnrichedNodes);

            var rightEnrichedNodes = Model.Nodes;
            PossibleEnrichments.Add(rightCrack.CrackBodyEnrichment, rightEnrichedNodes);
            PossibleEnrichments.Add(rightCrack.CrackTipEnrichments, rightEnrichedNodes);
        }

        public class Builder : IBenchmarkBuilder
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
            /// A node that lies in the positive halfplane defined by the body level set of a crack, will be enriched with 
            /// heaviside enrichment if Apos/(Apos+Aneg) &gt; <see cref="HeavisideEnrichmentTolerance"/> where Apos, Aneg 
            /// are the subareas of its nodal support  in the positive and negative halfplanes respectively. Similarly a
            /// node in the negative halfplane will be enriched if Aneg/(Apos+Aneg) &gt; 
            /// <see cref="HeavisideEnrichmentTolerance"/>.
            /// </summary>
            public double HeavisideEnrichmentTolerance { get; set; } = 0.0001;

            /// <summary>
            /// Controls how large will the radius of the J-integral contour be. WARNING: errors are introduced if the J-integral 
            /// radius is larger than the length of the crack segments.
            /// </summary>
            public double JintegralRadiusOverElementSize { get; set; } = 2.0;

            /// <summary>
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            public PropagationLogger KnownLeftPropagation { get; set; } = null;

            /// <summary>
            /// If it isn't null, the crack propagtion described by the passed <see cref="PropagationLogger"/> will be enforced, 
            /// rather than the using the J-integral method to predict it.
            /// </summary>
            public PropagationLogger KnownRightPropagation { get; set; } = null;

            /// <summary>
            /// The absolute path of the directory where output vtk files with the crack path and the level set functions at 
            /// each iteration will be written. Leave it null to avoid the performance cost it will introduce.
            /// </summary>
            public string LeftLsmPlotDirectory { get; set; } = null;

            /// <summary>
            /// The absolute path of the directory where output vtk files with the crack path and the level set functions at 
            /// each iteration will be written. Leave it null to avoid the performance cost it will introduce.
            /// </summary>
            public string RightLsmPlotDirectory { get; set; } = null;

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
                return new Holes(meshPath, growthLength, JintegralRadiusOverElementSize, 
                    KnownLeftPropagation, KnownRightPropagation, LeftLsmPlotDirectory, RightLsmPlotDirectory, MaxIterations,
                    HeavisideEnrichmentTolerance);
            }
        }

        private class HolesBoundary : IDomainBoundary
        {
            public HolesBoundary()
            {
            }

            public bool IsInside(ICartesianPoint2D point)
            {
                // Shapes
                var rectHull = new RectangularBoundary(minX, maxX, minY, maxY);
                var leftCircle = new Circle2D(new CartesianPoint2D(leftHoleX, leftHoleY), holeRadius);
                var rightCircle = new Circle2D(new CartesianPoint2D(rightHoleX, rightHoleY), holeRadius);

                // Intrnal points lie inside the rectangle, but outside the circular holes.
                if (rectHull.IsInside(point))
                {
                    if (leftCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside) return true;
                    if (rightCircle.FindRelativePositionOfPoint(point) == CirclePointPosition.Outside) return true;
                }
                return false;
            }
        }
    }
}
