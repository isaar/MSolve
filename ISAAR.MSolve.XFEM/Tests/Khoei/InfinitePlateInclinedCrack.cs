using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.XFEM.Analysis;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class InfinitePlateInclinedCrack
    {
        private static readonly double E = 2.1e6; // kg/cm^2
        private static readonly double v = 0.3;
        private static readonly double tension = 2000; // kg/cm (actually kg/cm^2, but the extra cm is thickness)
        private static readonly double width = 100;
        private static readonly double crackLength = width / 25;
        private static readonly int elementsPerCoarseRegion = 10;


        public static void Main()
        {
            double jIntegralRadiusOverElementSize = 3;
            double crackLengthOverElementSize = 7;
            double[] crackAnglesDegrees = new double[] { 0, 10, 20, 30, 40, 50, 60, 70};
            int entries = crackAnglesDegrees.Length;

            double[] analyticJintegrals = new double[entries];
            double[] analyticSIFs1 = new double[entries];
            double[] analyticSIFs2 = new double[entries];
            double[] startJintegrals = new double[entries];
            double[] startSIFs1 = new double[entries];
            double[] startSIFs2 = new double[entries];
            double[] endJintegrals = new double[entries];
            double[] endSIFs1 = new double[entries];
            double[] endSIFs2 = new double[entries];

            for (int i = 0; i < entries; ++i)
            {
                var test = new InfinitePlateInclinedCrack(crackAnglesDegrees[i], 
                    crackLength / crackLengthOverElementSize,
                    jIntegralRadiusOverElementSize);
                IVectorOLD solution = test.Solve();
                PropagationResults results = test.Propagate(solution);
                analyticJintegrals[i] = (double)results.referenceJIntegral;
                analyticSIFs1[i] = (double)results.referenceSIF1;
                analyticSIFs2[i] = (double)results.referenceSIF2;
                startJintegrals[i] = (double)results.startJIntegral;
                startSIFs1[i] = (double)results.startSIF1;
                startSIFs2[i] = (double)results.startSIF2;
                endJintegrals[i] = (double)results.endJIntegral;
                endSIFs1[i] = (double)results.endSIF1;
                endSIFs2[i] = (double)results.endSIF2;
            }

            Console.WriteLine("******************************* RESULTS *******************************");
            Console.Write("Rows = crack angle (degrees) ");
            ArrayPrinter.PrintArrayAsColumn(crackAnglesDegrees);
            Console.WriteLine();

            Console.WriteLine("\nAnalytic J-integral:");
            ArrayPrinter.PrintArrayAsColumn(analyticJintegrals);
            Console.WriteLine();

            Console.WriteLine("\nAnalytic SIF mode I:");
            ArrayPrinter.PrintArrayAsColumn(analyticSIFs1);
            Console.WriteLine();

            Console.WriteLine("\nAnalytic SIF mode II:");
            ArrayPrinter.PrintArrayAsColumn(analyticSIFs2);
            Console.WriteLine();

            Console.WriteLine("\nStart tip J-integral:");
            ArrayPrinter.PrintArrayAsColumn(startJintegrals);
            Console.WriteLine();

            Console.WriteLine("\nStart tip SIF mode I:");
            ArrayPrinter.PrintArrayAsColumn(startSIFs1);
            Console.WriteLine();

            Console.WriteLine("\nStart tip SIF mode II:");
            ArrayPrinter.PrintArrayAsColumn(startSIFs2);
            Console.WriteLine();

            Console.WriteLine("\nEnd tip J-integral:");
            ArrayPrinter.PrintArrayAsColumn(endJintegrals);
            Console.WriteLine();

            Console.WriteLine("\nEnd tip SIF mode I:");
            ArrayPrinter.PrintArrayAsColumn(endSIFs1);
            Console.WriteLine();

            Console.WriteLine("\nEnd tip SIF mode II:");
            ArrayPrinter.PrintArrayAsColumn(endSIFs2);
            Console.WriteLine();

            //Console.WriteLine("Enter any key to exit");
            //Console.ReadLine();
        }

        private readonly double crackAngleDegrees;
        private readonly double fineElementSize;
        private readonly double jIntegralRadiusOverElementSize;

        private readonly SubmatrixChecker checker;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private Model2D model;
        private BasicExplicitInteriorCrack crack;

        private readonly PropagationResults results;

        public InfinitePlateInclinedCrack(double crackAngleDegrees, double fineElementSize, 
            double jIntegralRadiusOverElementSize)
        {
            this.crackAngleDegrees = crackAngleDegrees;
            this.fineElementSize = fineElementSize;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            checker = new SubmatrixChecker(1e-4);

            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlainStress(E, v, 1.0);
            results = new PropagationResults();
            CalculateAnalytic();
            

            CreateModel();
            HandleCrack();
        }

        private void CalculateAnalytic()
        {
            double equivalentE = globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;
            double angleRad = crackAngleDegrees * Math.PI / 180;
            double cosBeta = Math.Cos(angleRad);
            double sinBeta = Math.Sin(angleRad);
            double sif1 = tension * Math.Sqrt(Math.PI * crackLength / 2) * cosBeta * cosBeta;
            double sif2 = tension * Math.Sqrt(Math.PI * crackLength / 2) * sinBeta * cosBeta;
            results.referenceSIF1 = sif1;
            results.referenceSIF2 = sif2;
            results.referenceJIntegral = (sif1 * sif1 + sif2 * sif2) / equivalentE;
        }

        private void CreateModel()
        {
            model = new Model2D();
            Tuple<XNode2D[], List<XNode2D[]>> meshEntities = CreateUniformMesh(fineElementSize);
            //Tuple<XNode2D[], List<XNode2D[]>> meshEntities = CreateRectilinearMesh(fineElementSize);

            // Nodes
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);

            // Elements
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new XSimpleIntegration2D());
            //var integration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStress(E, v, 1.0);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration, jIntegration));
            }

            ApplyBoundaryConditions2(model);
        }

        private Tuple<XNode2D[], List<XNode2D[]>> CreateUniformMesh(double elementSize)
        {
            int elementsPerAxis = (int)(width / elementSize) + 1;
            var meshGenerator = new UniformRectilinearMeshGenerator(width, width, elementsPerAxis, elementsPerAxis);
            return meshGenerator.CreateMesh();
        }

        private Tuple<XNode2D[], List<XNode2D[]>> CreateRectilinearMesh(double fineElementSize)
        {
            var meshGenerator = new FineAtCenterRectilinearMeshGenerator();
            meshGenerator.domainLowerBounds = new double[] { 0, 0 };
            meshGenerator.domainUpperBounds = new double[] { width, width };
            meshGenerator.interestAreaDimensions = new double[] { 1.5 * crackLength, 1.5 * crackLength };
            meshGenerator.coarseElementCountPerRegion = elementsPerCoarseRegion;
            meshGenerator.fineElementSize = fineElementSize;
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
            //model.AddConstraint(bottomRightNode, StandardDOFType.X, 0.0);
            //model.AddConstraint(bottomRightNode, StandardDOFType.Y, 0.0);
            //model.AddConstraint(topRightNode, StandardDOFType.X, 0.0);
            //model.AddConstraint(topRightNode, StandardDOFType.Y, 0.0);

            // Loads: Top and bottom sides are subject to tension
            var distrubutor = new LoadDistributor();
            IReadOnlyList<XNode2D> bottomNodes = finder.FindNodesWithY(0.0);
            IReadOnlyList<XNode2D> topNodes = finder.FindNodesWithY(width);
            double[,] topLoads = distrubutor.DistributeLoad(topNodes, 0.0, tension);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                if (topNodes[i] != topLeftNode)
                //if ((topNodes[i] != topLeftNode) && (topNodes[i] != topRightNode))
                {
                    model.AddNodalLoad(topNodes[i], StandardDOFType.X, topLoads[i, 0]);
                    model.AddNodalLoad(topNodes[i], StandardDOFType.Y, topLoads[i, 1]);
                }
            }
            double[,] bottomLoads = distrubutor.DistributeLoad(bottomNodes, 0.0, -tension);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                if (bottomNodes[i] != bottomLeftNode)
                //if ((bottomNodes[i] != bottomLeftNode) && (bottomNodes[i] != bottomRightNode))
                {
                    model.AddNodalLoad(bottomNodes[i], StandardDOFType.X, bottomLoads[i, 0]);
                    model.AddNodalLoad(bottomNodes[i], StandardDOFType.Y, bottomLoads[i, 1]);
                }
            }
        }

        private void HandleCrack()
        {
            crack = new BasicExplicitInteriorCrack();
            var boundary = new RectangularBoundary(0.0, width, 0.0, width);
            crack.Mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements, boundary);

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            crack.StartTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Start);
            crack.EndTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.End);

            // Mesh geometry interaction
            double dx = crackLength / 2 * Math.Cos(crackAngleDegrees * Math.PI / 180);
            double dy = crackLength / 2 * Math.Sin(crackAngleDegrees * Math.PI / 180);
            var startTip = new CartesianPoint2D(0.5 * width - dx, 0.5 * width - dy);
            var endTip = new CartesianPoint2D(0.5 * width + dx, 0.5 * width + dy);
            crack.InitializeGeometry(startTip, endTip);
            crack.UpdateEnrichments();
        }

        private IVectorOLD Solve()
        {
            model.EnumerateDofs();
            var analysis = new LinearStaticAnalysisSkyline(model);
            analysis.Solve();
            return analysis.Solution;
        }

        private PropagationResults Propagate(IVectorOLD solution)
        {
            double[] totalConstrainedDisplacements = model.CalculateConstrainedDisplacements();
            double[] totalFreeDisplacements = new double[solution.Length];
            solution.CopyTo(totalFreeDisplacements, 0);
           

            // Start tip propagation
            var startPropagator = new Propagator(crack.Mesh, crack, CrackTipPosition.Start, jIntegralRadiusOverElementSize,
                new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
                new HomogeneousSIFCalculator(globalHomogeneousMaterial),
                new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.5 * crackLength));
            double startGrowthAngle, startGrowthIncrement;
            startPropagator.Propagate(model, totalFreeDisplacements, totalConstrainedDisplacements,
                out startGrowthAngle, out startGrowthIncrement);
            double startSIF1 = startPropagator.Logger.SIFsMode1[0];
            double startSIF2 = startPropagator.Logger.SIFsMode2[0];

            // End tip propagation
            var endPropagator = new Propagator(crack.Mesh, crack, CrackTipPosition.End, jIntegralRadiusOverElementSize,
               new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
               new HomogeneousSIFCalculator(globalHomogeneousMaterial),
               new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.5 * crackLength));
            double endGrowthAngle, endGrowthIncrement;
            endPropagator.Propagate(model, totalFreeDisplacements, totalConstrainedDisplacements,
                out endGrowthAngle, out endGrowthIncrement);
            double endSIF1 = endPropagator.Logger.SIFsMode1[0];
            double endSIF2 = endPropagator.Logger.SIFsMode2[0];

            // Return results
            double equivalentE = globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;
            results.startSIF1 = startSIF1;
            results.startSIF2 = startSIF2;
            results.startJIntegral = (startSIF1 * startSIF1 + startSIF2 * startSIF2) / equivalentE;
            results.endSIF1 = endSIF1;
            results.endSIF2 = endSIF2;
            results.endJIntegral = (endSIF1 * endSIF1 + endSIF2 * endSIF2) / equivalentE;

            return results;
        }
    }
}
