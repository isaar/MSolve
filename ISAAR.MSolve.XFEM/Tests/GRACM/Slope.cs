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
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class Slope: IBenchmark
    {
        public static Builder SetupBenchmark()
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Slope\Meshes\fillet.msh";
            string plotPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Slope\Plots";
            string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Slope\Timing\results.txt";
            //string meshPath = @"C:\Users\seraf\Desktop\GRACM\Slope\Meshes\fillet.msh";
            //string plotPath = @"C:\Users\seraf\Desktop\GRACM\Slope\Plots";
            //string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Slope\Timing";

            double growthLength = 2.0; // mm. Must be sufficiently larger than the element size.
            var builder = new Builder(meshPath, growthLength, timingPath);
            builder.LsmPlotDirectory = plotPath;
            builder.MaxIterations = 10;

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
        private const double t = 1.0;

        /// <summary>
        /// Young's modulus in Pa
        /// </summary>
        private const double E = 10e6;

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3;

        /// <summary>
        /// The material used for the J-integral computation. It msut be stored separately from individual element materials.
        /// </summary>
        private static readonly HomogeneousElasticMaterial2D globalHomogeneousMaterial =
            HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

        /// <summary>
        /// The maximum value that the effective SIF can reach before collapse occurs.
        /// </summary>
        private const double fractureToughness = double.MaxValue;

        /// <summary>
        /// Magnitude in lbs of "opening" point loads applied at the upper left and lower left corners of the domain.
        /// </summary>
        private const double load = 4e5; // N


        // Geometry (m)
        private const double yTop = 10.0;
        private const double yExtra = 10.4;
        private const double xRight = 20.0;
        private const double crackMouthX = 15.5;
        private const double crackMouthY = 10.0;
        private const double crackTipX = 15.5;
        private const double crackTipY = 9.5;

        #endregion


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

        private TrackingExteriorCrackLSM crack;
        private BiMesh2D mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private Slope(string meshPath, double growthLength, double jIntegralRadiusOverElementSize,
             PropagationLogger knownPropagation, string lsmOutputDirectory, int maxIterations)
        {
            this.meshPath = meshPath;
            this.growthLength = growthLength;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.knownPropagation = knownPropagation;
            this.lsmOutputDirectory = lsmOutputDirectory;
            this.maxIterations = maxIterations;
        }
        
        /// <summary>
        /// The crack geometry description
        /// </summary>
        public ICrackDescription Crack { get { return crack; } }

        public IDomainDecomposer Decomposer { get; private set; }

        public string Name { get { return "GRACM Slope"; } }

        public string PlotDirectory { get { return lsmOutputDirectory; } }

        public Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>> PossibleEnrichments { get; private set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public Model2D Model { get; private set; }

        public void Analyze(ISolver solver)
        {
            var analysis = new QuasiStaticAnalysis(Model, mesh, crack, solver, fractureToughness, maxIterations);
            analysis.Analyze();

            Console.WriteLine("Crack path:");
            foreach (var point in crack.CrackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
            Console.WriteLine();
            Console.WriteLine("Crack growth angles:");
            foreach (var angle in crack.CrackTipPropagators[crack.CrackTips[0]].Logger.GrowthAngles)
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

            // Fixed dofs
            foreach (var node in finder.FindNodesWithY(0.0))
            {
                Model.AddConstraint(node, DisplacementDof.X, 0.0);
                Model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }
            foreach (var node in finder.FindNodesWithX(xRight))
            {
                if (node.Y != 0.0) Model.AddConstraint(node, DisplacementDof.X, 0.0);
            }

            // Loads
            IReadOnlyList<XNode2D> loadedNodes = finder.FindNodesWithY(yTop).Where(node => node.X <= 14.0).ToArray();
            double pressure = load / loadedNodes.Count;
            foreach (var node in loadedNodes)
            {
                Model.AddNodalLoad(node, DisplacementDof.Y, -pressure);
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
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                Model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }

            // Mesh usable for crack-mesh interaction
            var vertices = new CartesianPoint2D[4];
            vertices[0] = new CartesianPoint2D(0.0, 0.0);
            vertices[1] = new CartesianPoint2D(xRight, 0.0);
            vertices[2] = new CartesianPoint2D(xRight, 10.0);
            vertices[3] = new CartesianPoint2D(10.0, 10.0);
            mesh = new BiMesh2D(Model.Nodes, Model.Elements, new PolygonalBoundary(vertices));
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
            IPropagator propagator;
            var actualPropagator = new Propagator(mesh,jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            if (knownPropagation != null) propagator = new FixedPropagator(knownPropagation, null);
            else propagator = actualPropagator;

            var crackVertex0 = new CartesianPoint2D(crackMouthX, crackMouthY);
            var crackVertex1 = new CartesianPoint2D(crackTipX, crackTipY);
            var initialCrack = new PolyLine2D(crackVertex0, crackVertex1);
            var lsmCrack = new TrackingExteriorCrackLSM(propagator);
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
            this.crack = lsmCrack;
        }

        private void LimitEnrichedArea()
        {
            var enrichedNodes = Model.Nodes;
            PossibleEnrichments = new Dictionary<IEnrichmentItem2D, IReadOnlyList<XNode2D>>();
            PossibleEnrichments.Add(crack.CrackBodyEnrichment, enrichedNodes);
            PossibleEnrichments.Add(crack.CrackTipEnrichments, enrichedNodes);
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
                this.TimingOutputDirectory = timingPath;
            }

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
            public string LsmPlotDirectory { get; set; } = null;

            /// <summary>
            /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
            /// boundary or if the fracture toughness is exceeded.
            /// </summary>
            public int MaxIterations { get; set; } = int.MaxValue;

            public IBenchmark BuildBenchmark()
            {
                return new Slope(meshPath, growthLength, JintegralRadiusOverElementSize, KnownPropagation, LsmPlotDirectory,
                    MaxIterations);
            }

            /// <summary>
            /// The absolute path of the file where slover timing will be written.
            /// </summary>
            public string TimingOutputDirectory { get; }
        }
    }
}
