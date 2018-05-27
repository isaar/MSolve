using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class InfinitePlate
    {
        private static readonly string outputFile = "InfinitePlate";
        private static readonly double E = 2.1e6; // kg/cm^2
        private static readonly double v = 0.3;
        private static readonly double tension = 2000; // kg/cm (actually kg/cm^2, but the extra cm is thickness)
        private static readonly double crackLength = 4;
        private static readonly double width = 200;
        private static readonly int elementsPerCoarseRegion = 10;
        private static readonly double fineAreaExtentsOverCrackLength = 2.0;
        private static readonly bool uniformMesh = false;
        private static readonly bool constrainTopLeft = false, constrainTopRight = true;
        private static readonly bool constrainBottomLeft = false, constrainBottomRight = true;
        private static readonly bool clampedBottom = false; // overwrites the flags for corner nodes above

        public static void Main()
        {
            double defaultCrackLengthOverElementSize = 4.001;
            double defaultJIntegralRadusOverElementSize = 3.0;
            double defaultEnrichmentRadiusOverElementSize = 0.0;
            double defaultCrackAngleDegrees = 0.0;

            if (true)
            {
                var test = new InfinitePlate(defaultCrackAngleDegrees, crackLength / 2.01,
                    defaultJIntegralRadusOverElementSize, defaultEnrichmentRadiusOverElementSize);
                test.PrintMesh();
            }

            Console.WriteLine("Results:");
            if (uniformMesh) Console.WriteLine("Mesh: uniform");
            else Console.WriteLine("Mesh: rectilinear");
            if (clampedBottom) Console.WriteLine("Constraints: x,y on all bottom edge nodes");
            else
            {
                Console.Write("Constraints: ");
                Console.WriteLine("Top left node: x = {0}, y = {0}",(constrainTopLeft ? "constrained" : "free"));
                Console.WriteLine("Top right node: x = {0}, y = {0}", (constrainTopRight ? "constrained" : "free"));
                Console.WriteLine("Bottom left node: x = {0}, y = {0}", (constrainBottomLeft ? "constrained" : "free"));
                Console.WriteLine("Bottom right node: x = {0}, y = {0}", (constrainBottomRight ? "constrained" : "free"));
            }

            // Parametric mesh size
            //double[] crackLengthsOverElementSize = new double[] { 1.001, 2.001, 3.001, 4.001, 5.001, 6.001 };
            //double[] crackLengthsOverElementSize = new double[] { 1.001, 2.001, 3.001, 4.001, 5.001, 6.001, 7.001, 8.001, 9.001, 10.001 };
            double[] crackLengthsOverElementSize = new double[] { 1.001, 1.501, 2.001, 2.501, 3.001, 3.501, 4.001, 4.501, 5.001, 5.501, 6.001};
            Console.WriteLine();
            Console.WriteLine(" -------------- Parametric mesh size --------------");
            Console.Write("Crack angle (degrees) = " + defaultCrackAngleDegrees);
            Console.Write(", J-integral radius / element size = " + defaultJIntegralRadusOverElementSize);
            Console.Write(", Fixed enrichment area radius / element size = " + defaultEnrichmentRadiusOverElementSize);
            Console.WriteLine();
            Console.WriteLine("Crack length / element size | analytic J | analytic K1 | analytic K2 | "
                    + "start J | start K1 | start K2 | end J | end K1 | end K2");
            for (int i = 0; i < crackLengthsOverElementSize.Length; ++i)
            {
                var test = new InfinitePlate(defaultCrackAngleDegrees, crackLength / crackLengthsOverElementSize[i],
                    defaultJIntegralRadusOverElementSize, defaultEnrichmentRadiusOverElementSize);
                Vector solution = test.Solve();
                Report report = test.Propagate(solution);
                Console.Write(crackLength / report.tipElementSize);
                Console.Write(" " + report.analyticJIntegral);
                Console.Write(" " + report.analyticSIF1);
                Console.Write(" " + report.analyticSIF2);
                Console.Write(" " + report.startJIntegral);
                Console.Write(" " + report.startSIF1);
                Console.Write(" " + report.startSIF2);
                Console.Write(" " + report.endJIntegral);
                Console.Write(" " + report.endSIF1);
                Console.WriteLine(" " + report.endSIF2);
            }

            // Parametric J-integral radius
            double[] jIntegralRadiiOverElementSize = new double[] { 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5.0, 5.5, 6.0};
            Console.WriteLine();
            Console.WriteLine(" -------------- Parametric J-integral radius --------------");
            Console.Write("Crack angle (degrees) = " + defaultCrackAngleDegrees);
            Console.Write(", Crack length / element size = " + defaultCrackLengthOverElementSize);
            Console.Write(", Fixed enrichment area radius / element size = " + defaultEnrichmentRadiusOverElementSize);
            Console.WriteLine();
            Console.WriteLine("J-integral radius / element size | analytic J | analytic K1 | analytic K2 | "
                    + "start J | start K1 | start K2 | end J | end K1 | end K2");
            for (int i = 0; i < jIntegralRadiiOverElementSize.Length; ++i)
            {
                var test = new InfinitePlate(defaultCrackAngleDegrees, crackLength / defaultCrackLengthOverElementSize,
                    jIntegralRadiiOverElementSize[i], defaultEnrichmentRadiusOverElementSize);
                Vector solution = test.Solve();
                Report report = test.Propagate(solution);
                Console.Write(report.jIntegralRadiusOverElementSize);
                Console.Write(" " + report.analyticJIntegral);
                Console.Write(" " + report.analyticSIF1);
                Console.Write(" " + report.analyticSIF2);
                Console.Write(" " + report.startJIntegral);
                Console.Write(" " + report.startSIF1);
                Console.Write(" " + report.startSIF2);
                Console.Write(" " + report.endJIntegral);
                Console.Write(" " + report.endSIF1);
                Console.WriteLine(" " + report.endSIF2);
            }

            // Parametric fixed enrichment area radius
            double[] enrichmentRadiiOverElementSize = new double[] { };
            //double[] enrichmentRadiiOverElementSize = new double[] { 0, 1, 2, 3};
            Console.WriteLine();
            Console.WriteLine(" -------------- Parametric fixed enrichment area radius --------------");
            Console.Write("Crack angle (degrees) = " + defaultCrackAngleDegrees);
            Console.Write(", Crack length / element size = " + defaultCrackLengthOverElementSize);
            Console.Write(", J-integral radius / element size = " + defaultJIntegralRadusOverElementSize);
            Console.WriteLine();
            Console.WriteLine("Fixed enrichment area radius / element size | analytic J | analytic K1 | analytic K2 | "
                    + "start J | start K1 | start K2 | end J | end K1 | end K2");
            for (int i = 0; i < enrichmentRadiiOverElementSize.Length; ++i)
            {
                var test = new InfinitePlate(defaultCrackAngleDegrees, crackLength / defaultCrackLengthOverElementSize,
                    defaultJIntegralRadusOverElementSize, enrichmentRadiiOverElementSize[i]);
                Vector solution = test.Solve();
                Report report = test.Propagate(solution);
                Console.Write(report.enrichmentRadiusOverElementSize);
                Console.Write(" " + report.analyticJIntegral);
                Console.Write(" " + report.analyticSIF1);
                Console.Write(" " + report.analyticSIF2);
                Console.Write(" " + report.startJIntegral);
                Console.Write(" " + report.startSIF1);
                Console.Write(" " + report.startSIF2);
                Console.Write(" " + report.endJIntegral);
                Console.Write(" " + report.endSIF1);
                Console.WriteLine(" " + report.endSIF2);
            }

            // Parametric crack angle
            double[] crackAnglesDegrees = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
            Console.WriteLine();
            Console.WriteLine(" -------------- Parametric crack angle --------------");
            Console.Write("Crack length / element size = " + defaultCrackLengthOverElementSize);
            Console.Write(", J-integral radius / element size = " + defaultJIntegralRadusOverElementSize);
            Console.Write(", Fixed enrichment area radius / element size = " + defaultEnrichmentRadiusOverElementSize);
            Console.WriteLine();
            Console.WriteLine("Crack angle (degrees) / element size | analytic J | analytic K1 | analytic K2 | "
                    + "start J | start K1 | start K2 | end J | end K1 | end K2");
            for (int i = 0; i < crackAnglesDegrees.Length; ++i)
            {
                var test = new InfinitePlate(crackAnglesDegrees[i], crackLength / defaultCrackLengthOverElementSize,
                    defaultJIntegralRadusOverElementSize, defaultEnrichmentRadiusOverElementSize);
                Vector solution = test.Solve();
                Report report = test.Propagate(solution);
                Console.Write(report.crackAngleDegrees);
                Console.Write(" " + report.analyticJIntegral);
                Console.Write(" " + report.analyticSIF1);
                Console.Write(" " + report.analyticSIF2);
                Console.Write(" " + report.startJIntegral);
                Console.Write(" " + report.startSIF1);
                Console.Write(" " + report.startSIF2);
                Console.Write(" " + report.endJIntegral);
                Console.Write(" " + report.endSIF1);
                Console.WriteLine(" " + report.endSIF2);
            }
        }


        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private Model2D model;
        private IDofOrderer dofOrderer;
        private BasicExplicitInteriorCrack crack;
        private readonly Report report;
        private readonly double fineElementSize;

        public InfinitePlate(double crackAngleDegrees, double fineElementSize,
            double jIntegralRadiusOverElementSize, double enrichmentRadiusOverElementSize)
        {
            report = new Report();
            report.crackAngleDegrees = crackAngleDegrees;
            report.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            report.enrichmentRadiusOverElementSize = enrichmentRadiusOverElementSize;
            this.fineElementSize = fineElementSize;

            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, 1.0);
            CalculateAnalytic();

            CreateModel();
            HandleCrack();
        }

        private void CalculateAnalytic()
        {
            double equivalentE = globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;
            double angleRad = report.crackAngleDegrees * Math.PI / 180;
            double cosBeta = Math.Cos(angleRad);
            double sinBeta = Math.Sin(angleRad);
            double sif1 = tension * Math.Sqrt(Math.PI * crackLength / 2) * cosBeta * cosBeta;
            double sif2 = tension * Math.Sqrt(Math.PI * crackLength / 2) * sinBeta * cosBeta;
            report.analyticSIF1 = sif1;
            report.analyticSIF2 = sif2;
            report.analyticJIntegral = (sif1 * sif1 + sif2 * sif2) / equivalentE;
        }

        private void CreateModel()
        {
            model = new Model2D();
            XNode2D[] nodes;
            List<XNode2D[]> elementConnectivity;
            if (uniformMesh) (nodes, elementConnectivity) = CreateUniformMesh(fineElementSize);
            else (nodes, elementConnectivity) = CreateRectilinearMesh(fineElementSize);

            // Nodes
            foreach (XNode2D node in nodes) model.AddNode(node);

            // Elements
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4));
            //var integration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //        new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, 1.0);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration, jIntegration));
            }

            if (clampedBottom) ApplyBoundaryConditions1(model);
            else ApplyBoundaryConditions2(model);
        }

        private (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateUniformMesh(double elementSize)
        {
            int elementsPerAxis = (int)(width / elementSize) + 1;
            var meshGenerator = new UniformRectilinearMeshGenerator(width, width, elementsPerAxis, elementsPerAxis);
            return meshGenerator.CreateMesh();
        }

        private (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateRectilinearMesh(double fineElementSize)
        {
            var meshGenerator = new FineAtCenterRectilinearMeshGenerator();
            meshGenerator.domainLowerBounds = new double[] { 0, 0 };
            meshGenerator.domainUpperBounds = new double[] { width, width };
            meshGenerator.interestAreaDimensions = new double[] { fineAreaExtentsOverCrackLength * crackLength,
                    fineAreaExtentsOverCrackLength * crackLength };
            meshGenerator.coarseElementCountPerRegion = elementsPerCoarseRegion;
            meshGenerator.fineElementSize = fineElementSize;
            return meshGenerator.CreateMesh();
        }

        // Clamp the bottom edge
        private void ApplyBoundaryConditions1(Model2D model)
        {
            var finder = new EntityFinder(model, 1e-6);

            // Supports: Constrain x, y at all bottom nodes
            IReadOnlyList<XNode2D> bottomNodes = finder.FindNodesWithY(0.0);
            foreach (var node in bottomNodes)
            {
                model.AddConstraint(node, DisplacementDof.X, 0.0);
                model.AddConstraint(node, DisplacementDof.Y, 0.0);
            }

            // Loads: Top and bottom sides are subject to tension
            var distrubutor = new LoadDistributor();
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(width);
            double[,] topLoads = distrubutor.DistributeLoad(topNodes, 0.0, tension);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                model.AddNodalLoad(topNodes[i], DisplacementDof.X, topLoads[i, 0]);
                model.AddNodalLoad(topNodes[i], DisplacementDof.Y, topLoads[i, 1]);
            }
        }

        // Isostatic supports at top & bottom left corner nodes
        private void ApplyBoundaryConditions2(Model2D model)
        {
            var finder = new EntityFinder(model, 1e-6);
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            XNode2D bottomRightNode = finder.FindNodeWith(width, 0.0);
            XNode2D topLeftNode = finder.FindNodeWith(0.0, width);
            XNode2D topRightNode = finder.FindNodeWith(width, width);
            IReadOnlyList<XNode2D> bottomNodes = finder.FindNodesWithY(0.0);
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(width);

            // Supports: Constrain x, y at the corner nodes
            if (constrainBottomLeft)
            {
                model.AddConstraint(bottomLeftNode, DisplacementDof.X, 0.0);
                model.AddConstraint(bottomLeftNode, DisplacementDof.Y, 0.0);
            }
            if (constrainBottomRight)
            {
                model.AddConstraint(bottomRightNode, DisplacementDof.X, 0.0);
                model.AddConstraint(bottomRightNode, DisplacementDof.Y, 0.0);
            }
            if (constrainTopLeft)
            {
                model.AddConstraint(topLeftNode, DisplacementDof.X, 0.0);
                model.AddConstraint(topLeftNode, DisplacementDof.Y, 0.0);
            }
            if (constrainTopRight)
            {
                model.AddConstraint(topRightNode, DisplacementDof.X, 0.0);
                model.AddConstraint(topRightNode, DisplacementDof.Y, 0.0);
            }

            // Loads: Top and bottom sides are subject to tension
            var distrubutor = new LoadDistributor();
            double[,] topLoads = distrubutor.DistributeLoad(topNodes, 0.0, tension);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                if ((topNodes[i] == topLeftNode) && constrainTopLeft) continue;
                else if ((topNodes[i] == topRightNode) && constrainTopRight) continue;
                model.AddNodalLoad(topNodes[i], DisplacementDof.X, topLoads[i, 0]);
                model.AddNodalLoad(topNodes[i], DisplacementDof.Y, topLoads[i, 1]);
            }
            double[,] bottomLoads = distrubutor.DistributeLoad(bottomNodes, 0.0, -tension);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                if ((bottomNodes[i] == bottomLeftNode) && constrainBottomLeft) continue;
                else if ((bottomNodes[i] == bottomRightNode) && constrainBottomRight) continue;
                model.AddNodalLoad(bottomNodes[i], DisplacementDof.X, bottomLoads[i, 0]);
                model.AddNodalLoad(bottomNodes[i], DisplacementDof.Y, bottomLoads[i, 1]);
            }
        }

        private void HandleCrack()
        {
            throw new NotImplementedException();

            //crack = new BasicExplicitInteriorCrack(report.enrichmentRadiusOverElementSize);
            //var boundary = new RectangularBoundary(0.0, width, 0.0, width);
            //crack.Mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements, boundary);

            //// Create enrichments          
            //crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            //crack.StartTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Start);
            //crack.EndTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.End);

            //// Mesh geometry interaction
            //double dx = crackLength / 2 * Math.Cos(report.crackAngleDegrees * Math.PI / 180);
            //double dy = crackLength / 2 * Math.Sin(report.crackAngleDegrees * Math.PI / 180);
            //var startTip = new CartesianPoint2D(0.5 * width - dx, 0.5 * width - dy);
            //var endTip = new CartesianPoint2D(0.5 * width + dx, 0.5 * width + dy);
            //crack.InitializeGeometry(startTip, endTip);
            //crack.UpdateEnrichments();
        }

        private Vector Solve()
        {
            var solver = new SkylineSolver(model);
            solver.Initialize();
            solver.Solve();
            dofOrderer = solver.DofOrderer;
            return solver.Solution;
        }

        private Report Propagate(Vector solution)
        {
            throw new NotImplementedException();
            //Vector totalConstrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);

            //// Start tip propagation
            //var startPropagator = new Propagator(crack.Mesh, crack, CrackTipPosition.Start, 
            //    report.jIntegralRadiusOverElementSize,
            //    new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //    new HomogeneousSIFCalculator(globalHomogeneousMaterial),
            //    new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.5 * crackLength));
            //(double startGrowthAngle, double startGrowthIncrement) = startPropagator.Propagate(dofOrderer, solution,
            //    totalConstrainedDisplacements);
            //double startSIF1 = startPropagator.Logger.SIFsMode1[0];
            //double startSIF2 = startPropagator.Logger.SIFsMode2[0];

            //// End tip propagation
            //var endPropagator = new Propagator(crack.Mesh, crack, CrackTipPosition.End, 
            //    report.jIntegralRadiusOverElementSize,
            //    new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //    new HomogeneousSIFCalculator(globalHomogeneousMaterial),
            //    new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.5 * crackLength));
            //(double endGrowthAngle, double endGrowthIncrement) = endPropagator.Propagate(dofOrderer, solution,
            //    totalConstrainedDisplacements);
            //double endSIF1 = endPropagator.Logger.SIFsMode1[0];
            //double endSIF2 = endPropagator.Logger.SIFsMode2[0];

            //// Propagation results
            //double equivalentE = globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;
            //report.startSIF1 = startSIF1;
            //report.startSIF2 = startSIF2;
            //report.startJIntegral = (startSIF1 * startSIF1 + startSIF2 * startSIF2) / equivalentE;
            //report.endSIF1 = endSIF1;
            //report.endSIF2 = endSIF2;
            //report.endJIntegral = (endSIF1 * endSIF1 + endSIF2 * endSIF2) / equivalentE;

            //// Actual mesh size
            //var startTipElement = crack.GetTipElements(CrackTipPosition.Start)[0];
            //double startTipElementArea = ConvexPolygon2D.CreateUnsafe(startTipElement.Nodes).ComputeArea();
            //var endTipElement = crack.GetTipElements(CrackTipPosition.End)[0];
            //double endTipElementArea = ConvexPolygon2D.CreateUnsafe(endTipElement.Nodes).ComputeArea();
            //report.tipElementSize = Math.Sqrt(0.5 * (startTipElementArea + endTipElementArea));

            //return report;
        }

        private void PrintMesh()
        {
            var writer = new DiscontinuousMeshVTKWriter(model);
            writer.InitializeFile(outputFile);
            writer.CloseCurrentFile();
        }

        class Report
        {
            public double crackAngleDegrees;
            public double tipElementSize;
            public double jIntegralRadiusOverElementSize;
            public double enrichmentRadiusOverElementSize;
            public double? startJIntegral;
            public double? startSIF1;
            public double? startSIF2;
            public double endJIntegral;
            public double endSIF1;
            public double endSIF2;
            public double analyticJIntegral;
            public double analyticSIF1;
            public double analyticSIF2;
        }
    }
}
