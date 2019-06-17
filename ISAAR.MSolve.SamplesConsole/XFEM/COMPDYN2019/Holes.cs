using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests;
using static ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019.Utilities;

namespace ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019
{
    public class Holes
    {
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_4442.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_13738.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_22666.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_29052.msh";
        private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_189740.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_357324.msh";
        private const string leftCrackPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Plots\Left";
        private const string rightCrackPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Plots\Right";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\solver_log.txt";

        public static void Run()
        {
            HolesBenchmark benchmarkSub10;

            // Skyline
            //HolesBenchmark benchmarkSub1 = CreateSingleSubdomainBenchmark();
            //ISolver skylineSolver = DefineSolver(benchmarkSub1, SolverType.Skyline);
            //Console.WriteLine("Uncracked analysis, 1 subdomain, Skyline   : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub1.Model, skylineSolver));
            //Console.WriteLine("Cracked analysis only 1 step, 1 subdomain, Skyline   : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub1.Model, benchmarkSub1.Crack, skylineSolver));
            //Console.WriteLine("Skyline solver, 1 subdomain: ");
            //RunCrackPropagationAnalysis(benchmarkSub1, skylineSolver);

            // FETI-1 10 subdomains
            //benchmarkSub10 = CreateMultiSubdomainBenchmark(10);
            //ISolver solverFeti1 = DefineSolver(benchmarkSub10, SolverType.Feti1);
            //PlotSubdomains(subdomainPlotPath, benchmarkSub10.Model);
            //Console.WriteLine("Uncracked analysis, 10 subdomains, FETI-1  : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub10.Model, solverFeti1));
            //Console.WriteLine("Cracked analysis only 1 step, 10 subdomains, FETI-1  : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub10.Model, benchmarkSub10.Crack, solverFeti1));
            //Console.WriteLine("FETI-1, 10 subdomains: ");
            //RunCrackPropagationAnalysis(benchmarkSub10, solverFeti1);

            // FETI-DP 10 subdomains
            benchmarkSub10 = CreateMultiSubdomainBenchmark(10);
            ISolver solverFetiDP = DefineSolver(benchmarkSub10, SolverType.FetiDP);
            PlotSubdomains(benchmarkSub10, solverFetiDP);
            //Console.WriteLine("Uncracked analysis, 10 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub10.Model, solverFetiDP));
            //Console.WriteLine("Cracked analysis only 1 step, 10 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub10.Model, benchmarkSub10.Crack, solverFetiDP));
            //Console.WriteLine("FETI-DP, 10 subdomains: ");
            //RunCrackPropagationAnalysis(benchmarkSub10, solverFetiDP);

            Console.Write("\nEnd");
        }

        private static HolesBenchmark CreateMultiSubdomainBenchmark(int numSubdomains)
        {
            // Define subdomain boundaries
            double tol = 1E-13;
            var regions = new Dictionary<int, RectangularRegion2D>();
            if (numSubdomains == 10)
            {
                double boundaryY = 5.0;
                double minX = 0.0, minY = 0.0, maxX = 20.0, maxY = 10.0;
                double leftBoundary1X = 3.5, leftBoundary2X = 7.0, leftBoundary3X = 10.5, leftBoundary4X = 15.1;
                double rightBoundary1X = 4.9, rightBoundary2X = 9.5, rightBoundary3X = 13.0, rightBoundary4X = 16.5;

                // Left crack regions: 
                regions[0] = new RectangularRegion2D(minX, minY, leftBoundary1X, boundaryY, tol);
                regions[0].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[0].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[1] = new RectangularRegion2D(leftBoundary1X, minY, leftBoundary2X, boundaryY, tol);
                regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[2] = new RectangularRegion2D(leftBoundary2X, minY, leftBoundary3X, boundaryY, tol);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[3] = new RectangularRegion2D(leftBoundary3X, minY, leftBoundary4X, boundaryY, tol);
                regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[4] = new RectangularRegion2D(leftBoundary4X, minY, maxX, boundaryY, tol);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);

                // Right crack regions:
                regions[5] = new RectangularRegion2D(minX, boundaryY, rightBoundary1X, maxY, tol);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[6] = new RectangularRegion2D(rightBoundary1X, boundaryY, rightBoundary2X, maxY, tol);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[7] = new RectangularRegion2D(rightBoundary2X, boundaryY, rightBoundary3X, maxY, tol);
                regions[7].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[7].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[7].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[8] = new RectangularRegion2D(rightBoundary3X, boundaryY, rightBoundary4X, maxY, tol);
                regions[8].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[8].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[8].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

                regions[9] = new RectangularRegion2D(rightBoundary4X, boundaryY, maxX, maxY, tol);
                regions[9].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[9].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            }
            else
            {
                throw new NotImplementedException();
            }

            // Create the model, crack & mesh with only 1 subdomain
            HolesBenchmark benchmark = CreateSingleSubdomainBenchmark();

            // Partition mesh into subdomains
            var regionsGeneral = new Dictionary<int, IRegion2D>();
            foreach (var pair in regions) regionsGeneral[pair.Key] = pair.Value;
            var partitioner = new GuidedPartioner2D<XNode, XContinuumElement2D>(benchmark.Mesh, regionsGeneral);
            Dictionary<int, List<XContinuumElement2D>> elementsOfSubdomains = partitioner.CreateSubdomains();

            // Replace the single subdomain that was already created with the ones from the mesh partition
            benchmark.Model.Subdomains.Clear();
            foreach (int subdomainID in elementsOfSubdomains.Keys)
            {
                benchmark.Model.Subdomains.Add(subdomainID, new XSubdomain(subdomainID));
                benchmark.Model.Subdomains[subdomainID].Elements.AddRange(elementsOfSubdomains[subdomainID]);
            }

            return benchmark;
        }

        private static HolesBenchmark CreateSingleSubdomainBenchmark()
        {
            double growthLength = 1.0; // mm. Must be sufficiently larger than the element size.

            var builder = new HolesBenchmark.Builder(meshPath, growthLength);
            builder.LeftLsmPlotDirectory = leftCrackPlotDirectory;
            builder.RightLsmPlotDirectory = rightCrackPlotDirectory;
            builder.SubdomainPlotDirectory = subdomainPlotDirectory;

            builder.HeavisideEnrichmentTolerance = 0.12;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            // If you modify the following two parameters significantly, then you will need to redefine which nodes are expected 
            // to be enriched.
            builder.TipEnrichmentRadius = 0.5;
            builder.BC = HolesBenchmark.BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY;

            builder.MaxIterations = 12;

            HolesBenchmark benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static ISolver DefineSolver(HolesBenchmark benchmark, SolverType solverType)
        {
            if (solverType == SolverType.Skyline)
            {
                var builder = new SkylineSolver.Builder();
                builder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.Feti1)
            {
                var factorizationTolerances = new Dictionary<int, double>();
                if (benchmark.Model.Subdomains.Count == 10)
                {
                    factorizationTolerances[0] = 1E-2;
                    factorizationTolerances[1] = 1E-2;
                    factorizationTolerances[2] = 1E-2;
                    factorizationTolerances[3] = 1E-2;
                    factorizationTolerances[4] = 1E-2;
                    factorizationTolerances[5] = 1E-2;
                    factorizationTolerances[6] = 1E-2;
                    factorizationTolerances[7] = 1E-2;
                    factorizationTolerances[8] = 1E-2;
                    factorizationTolerances[9] = 1E-2;
                }
                else
                {
                    throw new NotImplementedException();
                }
                var fetiMatrices = new DenseFeti1SubdomainMatrixManager.Factory();
                var builder = new Feti1Solver.Builder(fetiMatrices, factorizationTolerances);
                builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                builder.ProblemIsHomogeneous = true;
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.FetiDP)
            {
                Dictionary<int, HashSet<INode>> cornerNodes = null;

                if (benchmark.Model.Subdomains.Count == 10)
                {
                    cornerNodes = FindCornerNodesFromCrosspoints2D(benchmark.Model);
                }
                else
                {
                    throw new NotImplementedException();
                }

                // Must also specify corner nodes
                //var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes);
                var cornerNodeSelection = new CrackedFetiDPSubdomainCornerNodes(benchmark.Crack, cornerNodes);
                var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
                var builder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
                builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                builder.ProblemIsHomogeneous = true;
                return builder.BuildSolver(benchmark.Model);
            }
            else throw new ArgumentException("Invalid solver choice.");
        }

        private static void PlotSubdomains(HolesBenchmark benchmark, ISolver solver)
        {
            string subdomainPlotPath = subdomainPlotDirectory + "\\subdomains.vtk";
            string boundaryNodesPlotPath = subdomainPlotDirectory + "\\boundary_nodes.vtk";
            string cornerNodesPlotPath = subdomainPlotDirectory + "\\corner_nodes.vtk";

            if (solver is Feti1Solver feti1)
            {
                benchmark.Model.ConnectDataStructures();
                var writer = new MeshPartitionWriter();
                writer.WriteSubdomainElements(subdomainPlotPath, benchmark.Model);
                writer.WriteBoundaryNodes(boundaryNodesPlotPath, benchmark.Model);
            }
            else if (solver is FetiDPSolver fetiDP)
            {
                benchmark.Model.ConnectDataStructures();
                var writer = new MeshPartitionWriter();
                writer.WriteSubdomainElements(subdomainPlotPath, benchmark.Model);
                writer.WriteBoundaryNodes(boundaryNodesPlotPath, benchmark.Model);

                var allCornerNodes = new HashSet<INode>();
                foreach (IEnumerable<INode> cornerNodes in fetiDP.CornerNodesOfSubdomains.Values)
                {
                    allCornerNodes.UnionWith(cornerNodes);
                }
                writer.WriteSpecialNodes(cornerNodesPlotPath, "corner_nodes", allCornerNodes);
            }
            else throw new ArgumentException("Invalid solver");
        }

        private static void RunCrackPropagationAnalysis(HolesBenchmark benchmark, ISolver solver)
        {
            benchmark.Analyze(solver);

            // Write crack path
            Console.WriteLine("Left Crack path:");
            foreach (var point in benchmark.LeftCrack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine("\n");
            Console.WriteLine("Right Crack path:");
            foreach (var point in benchmark.RightCrack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine();

            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }
    }
}
