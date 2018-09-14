using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
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
using ISAAR.MSolve.XFEM.Geometry.Polyline;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests.GRACM;
using ISAAR.MSolve.XFEM.Utilities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.XFEM.Tests.Visualization
{
    class IntersectedMeshCantilever
    {
        public static void Run()
        {
            string outputDirectory = @"C:\Users\Serafeim\Desktop\Presentation\XFEM\Plots";

            double growthLength = 0.3; //0.3, 0.4 is good
            double elementSize = 0.08; //0.08, 0.2 is good
            IMeshProvider meshProvider = new DCBUniformMeshProvider(elementSize);
            //double fineElementSize = 0.11;
            //IMeshProvider meshProvider = new DCBRefinedMeshProvider(fineElementSize, 10 * fineElementSize);
            //IMeshProvider meshProvider = new GmshMeshProvider(@"C: \Users\Serafeim\Desktop\GMSH\dcb.msh");
            int maxIterations = 10;
            double jIntegralRadiusOverElementSize = 2.0;

            var benchmark = new IntersectedMeshCantilever(growthLength, jIntegralRadiusOverElementSize, null,
                outputDirectory, outputDirectory, maxIterations, meshProvider);
            benchmark.InitializeModel();
            ISolver solver = new CholeskyAMDSolver(benchmark.Model);

            benchmark.Analyze(solver);
        }

        #region dimensions and material constants
        /// <summary>
        /// Horizontal dimension of the whole domain in inches.
        /// </summary>
        public const double L = 11.82;

        /// <summary>
        /// Vertical dimension of the whole domain in inches.
        /// </summary>
        public const double h = 3.94;

        /// <summary>
        /// Length of first crack segment (horizontal) in inches.
        /// </summary>
        public const double a = 3.95;

        /// <summary>
        /// Length of second crack segment (inclined) in inches.
        /// </summary>
        public const double da = 0.5;

        /// <summary>
        /// Initial crack angle in radians.
        /// </summary>
        public const double dTheta = 5.71 * Math.PI / 180;

        /// <summary>
        /// Thickness of the whole domain in inches.
        /// </summary>
        private const double t = 1.0;

        /// <summary>
        /// Young's modulus in psi=lbs/in^2
        /// </summary>
        private const double E = 3e7;

        /// <summary>
        /// Poisson's ratio
        /// </summary>
        private const double v = 0.3; // psi=lbs/in^2

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
        private const double load = 197; // lbs
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

        private readonly string lsmOutputDirectory;
        private readonly string fieldOutputDirectory;

        /// <summary>
        /// The maximum number of crack propagation steps. The analysis may stop earlier if the crack has reached the domain 
        /// boundary or if the fracture toughness is exceeded.
        /// </summary>
        private readonly int maxIterations;

        private readonly IMeshProvider meshProvider;

        private BiMesh2D mesh;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="growthLength">The length by which the crack grows in each iteration.</param>
        private IntersectedMeshCantilever(double growthLength, double jIntegralRadiusOverElementSize, PropagationLogger knownPropagation,
            string lsmOutputDirectory, string fieldOutputDirectory, int maxIterations, IMeshProvider meshProvider)
        {
            this.growthLength = growthLength;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.knownPropagation = knownPropagation;
            this.lsmOutputDirectory = lsmOutputDirectory;
            this.fieldOutputDirectory = fieldOutputDirectory;
            this.maxIterations = maxIterations;
            this.meshProvider = meshProvider;
        }

        /// <summary>
        /// The crack geometry description
        /// </summary>
        public TrackingExteriorCrackLSM Crack { get; private set; }

        /// <summary>
        /// Before accessing it, make sure <see cref="InitializeModel"/> has been called.
        /// </summary>
        public Model2D Model { get; private set; }


        public IReadOnlyList<ICartesianPoint2D> Analyze(ISolver solver)
        {
            var crackPath = new List<ICartesianPoint2D>();
            crackPath.Add(new CartesianPoint2D(0.0, h / 2));
            crackPath.Add(new CartesianPoint2D(a, h / 2));
            double dx = da * Math.Cos(dTheta);
            double dy = da * Math.Sin(dTheta);
            crackPath.Add(new CartesianPoint2D(crackPath.Last().X + dx, crackPath.Last().Y - dy));

            string fieldOuputPath = fieldOutputDirectory + "\\field_output";
            var fieldOutput = new IntersectedMeshOutput(Model, Crack, fieldOuputPath);
            var analysis = new QuasiStaticAnalysis(Model, mesh, Crack, solver, fractureToughness, maxIterations, fieldOutput);
            analysis.Analyze();
            crackPath.AddRange(Crack.CrackPath);
            return crackPath;
        }

        public void InitializeModel()
        {
            Model = new Model2D();
            CreateMesh();
            ApplyBoundaryConditions();
            InitializeCrack();
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(Model, 1e-6);

            // Fixed dofs
            foreach (var node in finder.FindNodesWithX(L))
            {
                Model.AddConstraint(node, DisplacementDof.X, 0.0);
                Model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }

            // Loads
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            XNode2D topLeftNode = finder.FindNodeWith(0.0, h);
            Model.AddNodalLoad(bottomLeftNode, DisplacementDof.Y, -load);
            Model.AddNodalLoad(topLeftNode, DisplacementDof.Y, load);
        }

        private void CreateMesh()
        {
            // Mesh generation
            (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) = meshProvider.CreateMesh();

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
            var boundary = new RectangularBoundary(0.0, L, 0.0, h);
            mesh = new BiMesh2D(Model.Nodes, Model.Elements, boundary);
        }

        private void InitializeCrack()
        {
            var crackVertex0 = new CartesianPoint2D(0.0, h / 2);
            var crackVertex1 = new CartesianPoint2D(a, h / 2);
            var initialCrack = new PolyLine2D(crackVertex0, crackVertex1);
            initialCrack.UpdateGeometry(-dTheta, da);

            var actualPropagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            IPropagator propagator;
            if (knownPropagation != null) propagator = new FixedPropagator(knownPropagation, null);
            else propagator = actualPropagator;

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

            this.Crack = lsmCrack;
        }
    }
}
