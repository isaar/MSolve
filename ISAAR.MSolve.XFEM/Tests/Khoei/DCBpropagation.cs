using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Analysis;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.LinearAlgebra;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class DCBpropagation
    {
        private static readonly double DIM_X = 60, DIM_Y = 20;
        private static readonly double E = 2e6, v = 0.3;
        private static readonly ICartesianPoint2D CRACK_MOUTH = new CartesianPoint2D(DIM_X, 0.5 * DIM_Y);
        private static readonly SubmatrixChecker checker = new SubmatrixChecker(1e-4);

        public static void Main()
        {
            // Parameters
            int elementsPerY = 15;
            double jIntegralRadiusOverElementSize = 2.0;
            double fractureTaoughness = double.MaxValue;
            int maxIterations = 10;
            double growthLength = 5; // cm

            var benchmark = new DCBpropagation(elementsPerY, jIntegralRadiusOverElementSize, fractureTaoughness,
                maxIterations, growthLength);
            benchmark.Analyze();
        }

        
        private Model2D model;
        private BasicExplicitCrack2D crack;
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

            CreateModel(3 * elementsPerY, elementsPerY);
            HandleCrack();
        }

        private void CreateModel(int elementsPerX, int elementsPerY)
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
            model = new Model2D();

            // Mesh
            var meshGenerator = new UniformRectilinearMeshGenerator(DIM_X, DIM_Y, elementsPerX, elementsPerY);
            Tuple<XNode2D[], List<XNode2D[]>> meshEntities = meshGenerator.CreateMesh();

            // Nodes
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);

            // Elements
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));

            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration, jIntegration));
            }

            // Boundary conditions
            var finder = new EntityFinder(model, 1e-6);
            var bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            var topLeftNode = finder.FindNodeWith(0.0, DIM_Y);
            var bottomRightNode = finder.FindNodeWith(DIM_X, 0.0);
            var topRightNode = finder.FindNodeWith(DIM_X, DIM_Y);

            model.AddConstraint(bottomLeftNode, StandardDOFType.X, 0.0);
            model.AddConstraint(bottomLeftNode, StandardDOFType.Y, 0.0);
            model.AddConstraint(topLeftNode, StandardDOFType.X, 0.0);
            model.AddConstraint(topLeftNode, StandardDOFType.Y, 0.0);
            model.AddConstraint(bottomRightNode, StandardDOFType.Y, -0.05);
            model.AddConstraint(topRightNode, StandardDOFType.Y, 0.05);
        }

        private void HandleCrack()
        {
            crack = new BasicExplicitCrack2D();
            var boundary = new RectangularBoundary(0.0, DIM_X, 0.0, DIM_Y);
            crack.Mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements, boundary);

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack, new SignFunctionOpposite2D());
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);

            // Mesh geometry interaction
            var crackTip = new CartesianPoint2D(0.5 * DIM_X, 0.5 * DIM_Y);
            crack.InitializeGeometry(CRACK_MOUTH, crackTip);
            crack.UpdateEnrichments();
        }

        private void Analyze()
        {
            propagator = new Propagator(crack.Mesh, crack, CrackTipPosition.Single, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(growthLength));

            QuasiStaticAnalysis analysis = new QuasiStaticAnalysis(model, crack.Mesh, crack, 
                propagator, fractureToughness, maxIterations);
            analysis.Analyze();
        }
    }
}
