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
    //TODO: Fix the bugs when a/s=4.
    class InfinitePlateHorizontalCrack
    {
        private static readonly double E = 2.1e6; // kg/cm^2
        private static readonly double v = 0.3;
        private static readonly double tension = 2000; // kg/cm (actually kg/cm^2, but the extra cm is thickness)
        private static readonly double width = 100;
        private static readonly double crackLength = width / 25;
        private static readonly int elementsPerCoarseRegion = 10;


        public static void Main()
        {
            var test = new InfinitePlateHorizontalCrack(0.5*crackLength, 2);
            IVector solution = test.Solve();
            //test.CheckSolution(solution);
            //test.Propagate(solution);
        }

        private readonly SubmatrixChecker checker;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private Model2D model;
        private BasicExplicitInteriorCrack crack;
        private readonly double finiteElementSize;
        private readonly double jIntegralRadiusOverElementSize;

        public InfinitePlateHorizontalCrack(double finiteElementSize, double jIntegralRadiusOverElementSize)
        {
            this.finiteElementSize = finiteElementSize;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            checker = new SubmatrixChecker(1e-4);
            CreateModel();
            HandleCrack();
        }

        private void CreateModel()
        {
            model = new Model2D();
            //Tuple<XNode2D[], List<XNode2D[]>> meshEntities = CreateUniformMesh(45);
            Tuple<XNode2D[], List<XNode2D[]>> meshEntities = CreateRectilinearMesh(finiteElementSize);

            // Nodes
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);

            // Elements
            var integration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration, jIntegration));
            }

            ApplyBoundaryConditions(model);
        }

        private Tuple<XNode2D[], List<XNode2D[]>> CreateUniformMesh(int elementsPerAxis)
        {
            var meshGenerator = new UniformRectilinearMeshGenerator(width, width, elementsPerAxis, elementsPerAxis);
            return meshGenerator.CreateMesh();
        }

        private Tuple<XNode2D[], List<XNode2D[]>> CreateRectilinearMesh(double finiteElementSize)
        {
            var meshGenerator = new FineAtCenterRectilinearMeshGenerator();
            meshGenerator.domainLowerBounds = new double[] { 0, 0 };
            meshGenerator.domainUpperBounds = new double[] { width, width };
            meshGenerator.interestAreaDimensions = new double[] { crackLength, crackLength };
            meshGenerator.coarseElementCountPerRegion = elementsPerCoarseRegion;
            meshGenerator.fineElementSize = finiteElementSize;
            return meshGenerator.CreateMesh();
        }

        private void ApplyBoundaryConditions(Model2D model)
        {
            var finder = new EntityFinder(model, 1e-6);

            // Supports: Constrain x, y at all bottom nodes
            IReadOnlyList<XNode2D> bottomNodes = finder.FindNodesWithY(0.0);
            foreach (var node in bottomNodes)
            {
                model.AddConstraint(node, StandardDOFType.X, 0.0);
                model.AddConstraint(node, StandardDOFType.Y, 0.0);
            }

            // Loads: Top and bottom sides are subject to tension
            var distrubutor = new LoadDistributor();
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(width);
            double[,] topLoads = distrubutor.DistributeLoad(topNodes, 0.0, tension);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                model.AddNodalLoad(topNodes[i], StandardDOFType.X, topLoads[i, 0]);
                model.AddNodalLoad(topNodes[i], StandardDOFType.Y, topLoads[i, 1]);
            }
        }

        private void ApplyBoundaryConditions2(Model2D model)
        {
            var finder = new EntityFinder(model, 1e-6);
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            XNode2D topLeftNode = finder.FindNodeWith(0.0, width);
            XNode2D bottomRightNode = finder.FindNodeWith(width, 0.0);
            XNode2D topRightNode = finder.FindNodeWith(width, width);

            // Supports: Constrain x, y at the corner nodes
            model.AddConstraint(bottomLeftNode, StandardDOFType.X, 0.0);
            model.AddConstraint(bottomLeftNode, StandardDOFType.Y, 0.0);
            model.AddConstraint(topLeftNode, StandardDOFType.X, 0.0);
            model.AddConstraint(topLeftNode, StandardDOFType.Y, 0.0);
            model.AddConstraint(bottomRightNode, StandardDOFType.X, 0.0);
            model.AddConstraint(bottomRightNode, StandardDOFType.Y, 0.0);
            model.AddConstraint(topRightNode, StandardDOFType.X, 0.0);
            model.AddConstraint(topRightNode, StandardDOFType.Y, 0.0);

            // Loads: Top and bottom sides are subject to tension
            IEnumerable<XNode2D> bottomNodes =
                finder.FindNodesWithY(0).Except(new XNode2D[] { bottomLeftNode, bottomRightNode });
            IEnumerable<XNode2D> topNodes = 
                finder.FindNodesWithY(width).Except(new XNode2D[] { topLeftNode, topRightNode });
            throw new NotImplementedException("Distribute the load to the nodes");
        }

        private void HandleCrack()
        {
            crack = new BasicExplicitInteriorCrack();
            crack.Mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements);

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            crack.StartTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Start);
            crack.EndTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.End);

            // Mesh geometry interaction
            var startTip = new CartesianPoint2D(0.5 * width - 0.5 * crackLength, 0.5 * width);
            var endTip = new CartesianPoint2D(0.5 * width + 0.5 * crackLength, 0.5 * width);
            crack.InitializeGeometry(startTip, endTip);
            crack.UpdateEnrichments();
        }

        private IVector Solve()
        {
            model.EnumerateDofs();
            var analysis = new LinearStaticAnalysis(model);
            analysis.Solve();
            return analysis.Solution;
        }
    }
}
