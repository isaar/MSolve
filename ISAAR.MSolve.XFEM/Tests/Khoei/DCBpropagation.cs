using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class DCBpropagation
    {
        private static readonly double DIM_X = 60, DIM_Y = 20;
        private static readonly double E = 2e6, v = 0.3;
        private static readonly ICartesianPoint2D CRACK_MOUTH = new CartesianPoint2D(DIM_X, 0.5 * DIM_Y);
        private static readonly bool structuredMesh = true;
        private static readonly string triMeshFile = "dcb_tri.msh";
        private static readonly string quadMeshFile = "dcb_quad.msh";
        private static readonly bool integrationWithTriangles = false;

        public static void Main()
        {
            // Parameters
            int elementsPerY = 15;
            double jIntegralRadiusOverElementSize = 2.0;
            double fractureToughness = double.MaxValue;
            int maxIterations = 10;
            double growthLength = 2; // cm

            var benchmark = new DCBpropagation(elementsPerY, jIntegralRadiusOverElementSize, fractureToughness,
                maxIterations, growthLength);
            benchmark.Analyze();
        }

        
        private Model2D model;
        private TrackingExteriorCrackLSM crack;
        private IIntegrationStrategy2D<XContinuumElement2D> integration;
        private IIntegrationStrategy2D<XContinuumElement2D> jIntegration;
        private readonly double jIntegralRadiusOverElementSize;
        private Propagator propagator;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private readonly double fractureToughness;
        private readonly int maxIterations;
        private readonly double growthLength;

        public DCBpropagation(int elementsPerY, double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations, double growthLength)
        {
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
            this.growthLength = growthLength;

            CreateModel(elementsPerY);
        }

        private void CreateModel(int elementsPerY)
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
            model = new Model2D();
            HandleIntegrations();
            if (structuredMesh) CreateStructuredMesh(elementsPerY);
            else CreateUnstructuredMesh();
            ApplyBoundaryConditions();
            HandleCrack();
        }

        private void HandleIntegrations()
        {
            if (integrationWithTriangles)
            {
                ITriangulator2D triangulator = new IncrementalTriangulator();
                integration = new IntegrationForCrackPropagation2D(
                    new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order2Points3, crack, triangulator),
                    new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order2Points3, crack, triangulator));
                jIntegration = new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order3Points4, crack,
                    triangulator);
            }
            else
            {
                integration = new IntegrationForCrackPropagation2D(
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
                jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);
                //var jIntegration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
                //    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4));
            }
        }

        private void CreateStructuredMesh(int elementsPerY)
        {
            var meshGenerator = new UniformRectilinearMeshGenerator(DIM_X, DIM_Y, 3 * elementsPerY, elementsPerY);
            (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) = meshGenerator.CreateMesh();
            foreach (XNode2D node in nodes) model.AddNode(node);
            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }
        }

        private void CreateUnstructuredMesh()
        {
            //var meshReader = new GmshReader(triMeshFile);
            var meshReader = new GmshReader(quadMeshFile);
            Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> meshEntities = meshReader.ReadMesh();
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);
            foreach (var element in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                model.AddElement(new XContinuumElement2D(element.ElementType, element.Nodes, materialField,
                    integration, jIntegration));
            }
        }

        private void ApplyBoundaryConditions()
        {
            var finder = new EntityFinder(model, 1e-6);

            // Fixed dofs
            foreach (var node in finder.FindNodesWithX(0.0))
            {
                model.AddConstraint(node, DisplacementDof.X, 0.0);
                model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }

            // Prescribed displacements
            XNode2D bottomRightNode = finder.FindNodeWith(DIM_X, 0.0);
            XNode2D topRightNode = finder.FindNodeWith(DIM_X, DIM_Y);
            model.AddConstraint(bottomRightNode, DisplacementDof.Y, -0.05);
            model.AddConstraint(topRightNode, DisplacementDof.Y, 0.05);
        }

        private void HandleCrack()
        {
            var boundary = new RectangularBoundary(0.0, DIM_X, 0.0, DIM_Y);
            var mesh = new BiMesh2D(model.Nodes, model.Elements, boundary);
            propagator = new Propagator(mesh, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));
            crack = new TrackingExteriorCrackLSM(propagator);
            crack.Mesh = mesh;

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack, new SignFunctionOpposite2D());
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);

            // Mesh geometry interaction
            var crackTip = new CartesianPoint2D(0.5 * DIM_X, 0.5 * DIM_Y);
            crack.InitializeGeometry(CRACK_MOUTH, crackTip);
        }

        private void Analyze()
        {
            var solver = new SkylineSolver(model);
            QuasiStaticAnalysis analysis = new QuasiStaticAnalysis(model, crack.Mesh, crack, solver, 
                fractureToughness, maxIterations);
            analysis.Analyze();
        }
    }
}
