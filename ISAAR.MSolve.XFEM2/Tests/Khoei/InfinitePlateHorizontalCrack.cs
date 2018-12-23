using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.Mesh.SuiteSparse;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests.Tools;
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
        //private static readonly double crackLength = width / 50;
        private static readonly int elementsPerCoarseRegion = 10;


        public static void Main()
        {
            double[] crackLengthsOverElementSize = new double[] { 1.05, 2.05, 3, 4.05, 5, 6.05, 7, 8.05, 9, 10.05 };
            //double[] crackLengthsOverElementSize = new double[] { 1.05, 2, 3, 4.05, 5, 6};
            //double[] jIntegralRadiiOverElementSize = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double[] jIntegralRadiiOverElementSize = new double[] { 3 };
            int rows = crackLengthsOverElementSize.Length;
            int columns = jIntegralRadiiOverElementSize.Length;

            double[,] analyticJintegrals = new double[rows, columns];
            double[,] analyticSIFs1 = new double[rows, columns];
            double[,] analyticSIFs2 = new double[rows, columns];
            double[,] startJintegrals = new double[rows, columns];
            double[,] startSIFs1 = new double[rows, columns];
            double[,] startSIFs2 = new double[rows, columns];
            double[,] endJintegrals = new double[rows, columns];
            double[,] endSIFs1 = new double[rows, columns];
            double[,] endSIFs2 = new double[rows, columns];

            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < columns; ++j)
                {
                    var test = new InfinitePlateHorizontalCrack(crackLength / crackLengthsOverElementSize[i],
                        jIntegralRadiiOverElementSize[j]);
                    Vector solution = test.Solve();
                    PropagationResults results = test.Propagate(solution);
                    analyticJintegrals[i, j] = (double) results.referenceJIntegral;
                    analyticSIFs1[i, j] = (double)results.referenceSIF1;
                    analyticSIFs2[i, j] = (double)results.referenceSIF2;
                    startJintegrals[i, j] = (double)results.startJIntegral;
                    startSIFs1[i, j] = (double)results.startSIF1;
                    startSIFs2[i, j] = (double)results.startSIF2;
                    endJintegrals[i, j] = (double)results.endJIntegral;
                    endSIFs1[i, j] = (double)results.endSIF1;
                    endSIFs2[i, j] = (double)results.endSIF2;
                }
            }
            
            Console.WriteLine("******************************* RESULTS *******************************");

            Console.Write("Rows = crack length over element size: ");
            ArrayPrinter.PrintArrayAsColumn(crackLengthsOverElementSize);
            Console.WriteLine();
            Console.Write("Columns = J-integral radius over element size: ");
            ArrayPrinter.PrintArrayAsRow(jIntegralRadiiOverElementSize);
            Console.WriteLine();

            Console.WriteLine("\nAnalytic J-integral:");
            ArrayPrinter.PrintArray2D(analyticJintegrals);
            Console.WriteLine();

            Console.WriteLine("\nAnalytic SIF mode I:");
            ArrayPrinter.PrintArray2D(analyticSIFs1);
            Console.WriteLine();

            Console.WriteLine("\nAnalytic SIF mode II:");
            ArrayPrinter.PrintArray2D(analyticSIFs2);
            Console.WriteLine();

            Console.WriteLine("\nStart tip J-integral:");
            ArrayPrinter.PrintArray2D(startJintegrals);
            Console.WriteLine();

            Console.WriteLine("\nStart tip SIF mode I:");
            ArrayPrinter.PrintArray2D(startSIFs1);
            Console.WriteLine();

            Console.WriteLine("\nStart tip SIF mode II:");
            ArrayPrinter.PrintArray2D(startSIFs2);
            Console.WriteLine();

            Console.WriteLine("\nEnd tip J-integral:");
            ArrayPrinter.PrintArray2D(endJintegrals);
            Console.WriteLine();

            Console.WriteLine("\nEnd tip SIF mode I:");
            ArrayPrinter.PrintArray2D(endSIFs1);
            Console.WriteLine();

            Console.WriteLine("\nEnd tip SIF mode II:");
            ArrayPrinter.PrintArray2D(endSIFs2);
            Console.WriteLine();

            //Console.WriteLine("Enter any key to exit");
            //Console.ReadLine();
        }

        private readonly SubmatrixChecker checker;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private Model2D model;
        private IDofOrderer dofOrderer;
        private BasicExplicitInteriorCrack crack;
        private readonly double fineElementSize;
        private readonly double jIntegralRadiusOverElementSize;

        private readonly PropagationResults results;

        public InfinitePlateHorizontalCrack(double fineElementSize, double jIntegralRadiusOverElementSize)
        {
            this.fineElementSize = fineElementSize;
            this.jIntegralRadiusOverElementSize = jIntegralRadiusOverElementSize;
            checker = new SubmatrixChecker(1e-4);

            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, 1.0);
            results = new PropagationResults();
            CalculateAnalytic();

            CreateModel();
            HandleCrack();
        }

        private void CalculateAnalytic()
        {
            double equivalentE = globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;
            double sif1 = tension * Math.Sqrt(Math.PI * crackLength / 2);
            double sif2 = 0.0;
            results.referenceSIF1 = sif1;
            results.referenceSIF2 = sif2;
            results.referenceJIntegral = (sif1 * sif1 + sif2 * sif2) / equivalentE;
        }

        private void CreateModel()
        {
            model = new Model2D();
            (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) = CreateUniformMesh(fineElementSize);
            //Tuple<XNode2D[], List<XNode2D[]>> meshEntities = CreateRectilinearMesh(fineElementSize);

            // Nodes
            foreach (XNode2D node in nodes) model.AddNode(node);

            // Elements
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new XSimpleIntegration2D());
            //var integration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //        new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, 1.0);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration, jIntegration));
            }

            //ApplyBoundaryConditions(model);
            ApplyBoundaryConditions2(model);
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

        private void ApplyBoundaryConditions2(Model2D model)
        {
            var finder = new EntityFinder(model, 1e-6);
            XNode2D bottomLeftNode = finder.FindNodeWith(0.0, 0.0);
            XNode2D topLeftNode = finder.FindNodeWith(0.0, width);
            XNode2D bottomRightNode = finder.FindNodeWith(width, 0.0);
            XNode2D topRightNode = finder.FindNodeWith(width, width);

            // Supports: Constrain x, y at the corner nodes
            model.AddConstraint(bottomLeftNode, DisplacementDof.X, 0.0);
            model.AddConstraint(bottomLeftNode, DisplacementDof.Y, 0.0);
            model.AddConstraint(topLeftNode, DisplacementDof.X, 0.0);
            model.AddConstraint(topLeftNode, DisplacementDof.Y, 0.0);
            //model.AddConstraint(bottomRightNode, DisplacementDof.X, 0.0);
            //model.AddConstraint(bottomRightNode, DisplacementDof.Y, 0.0);
            //model.AddConstraint(topRightNode, DisplacementDof.X, 0.0);
            //model.AddConstraint(topRightNode, DisplacementDof.Y, 0.0);

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
                    model.AddNodalLoad(topNodes[i], DisplacementDof.X, topLoads[i, 0]);
                    model.AddNodalLoad(topNodes[i], DisplacementDof.Y, topLoads[i, 1]);
                }
            }
            double[,] bottomLoads = distrubutor.DistributeLoad(bottomNodes, 0.0, -tension);
            for (int i = 0; i < topNodes.Count; ++i)
            {
                if (bottomNodes[i] != bottomLeftNode)
                //if ((bottomNodes[i] != bottomLeftNode) && (bottomNodes[i] != bottomRightNode))
                {
                    model.AddNodalLoad(bottomNodes[i], DisplacementDof.X, bottomLoads[i, 0]);
                    model.AddNodalLoad(bottomNodes[i], DisplacementDof.Y, bottomLoads[i, 1]);
                }
            }
        }

        private void HandleCrack()
        {
            throw new NotImplementedException();

            //crack = new BasicExplicitInteriorCrack(3.0 * fineElementSize);
            //var boundary = new RectangularBoundary(0.0, width, 0.0, width);
            //crack.Mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements, boundary);

            //// Create enrichments          
            //crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            //crack.StartTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Start);
            //crack.EndTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.End);

            //// Mesh geometry interaction
            //var startTip = new CartesianPoint2D(0.5 * width - 0.5 * crackLength, 0.5 * width);
            //var endTip = new CartesianPoint2D(0.5 * width + 0.5 * crackLength, 0.5 * width);
            //crack.InitializeGeometry(startTip, endTip);
            //crack.UpdateEnrichments();
        }

        private Vector Solve()
        {
            var solver = new SkylineSolver(model);
            solver.Initialize();
            solver.Solve();
            dofOrderer = solver.DofOrderer;
            return (Vector)solver.Solution;
        }

        private PropagationResults Propagate(Vector solution)
        {
            throw new NotImplementedException();

            //Vector totalConstrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);

            //// Start tip propagation
            //var startPropagator = new Propagator(crack.Mesh, crack, CrackTipPosition.Start, jIntegralRadiusOverElementSize,
            //    new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //    new HomogeneousSIFCalculator(globalHomogeneousMaterial),
            //    new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.5 * crackLength));
            //(double startGrowthAngle, double startGrowthIncrement) = startPropagator.Propagate(dofOrderer, solution, 
            //    totalConstrainedDisplacements);
            //double startSIF1 = startPropagator.Logger.SIFsMode1[0];
            //double startSIF2 = startPropagator.Logger.SIFsMode2[0];            

            //// End tip propagation
            //var endPropagator = new Propagator(crack.Mesh, crack, CrackTipPosition.End, jIntegralRadiusOverElementSize,
            //   new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //   new HomogeneousSIFCalculator(globalHomogeneousMaterial),
            //   new MaximumCircumferentialTensileStressCriterion(), new ConstantIncrement2D(0.5 * crackLength));
            //(double endGrowthAngle, double endGrowthIncrement) = endPropagator.Propagate(dofOrderer, solution,
            //    totalConstrainedDisplacements);
            //double endSIF1 = endPropagator.Logger.SIFsMode1[0];
            //double endSIF2 = endPropagator.Logger.SIFsMode2[0];

            //// Return results
            //double equivalentE = globalHomogeneousMaterial.HomogeneousEquivalentYoungModulus;
            //results.startSIF1 = startSIF1;
            //results.startSIF2 = startSIF2;
            //results.startJIntegral = (startSIF1 * startSIF1 + startSIF2 * startSIF2) / equivalentE;
            //results.endSIF1 = endSIF1;
            //results.endSIF2 = endSIF2;
            //results.endJIntegral = (endSIF1 * endSIF1 + endSIF2 * endSIF2) / equivalentE;
            //return results;
        }
    }
}
