using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
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
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests;
using static ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019.Utilities;

namespace ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019
{
    public class Fillet
    {
        private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Fillet\Mesh\fillet.msh";
        private const string crackPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Fillet\Plots";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Fillet\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Fillet\solver_log.txt";

        public static void Run()
        {
            FilletBenchmark benchmarkSub1;
            FilletBenchmark benchmarkSub5;
            FilletBenchmark benchmarkSub7;

            // Skyline
            benchmarkSub1 = CreateSingleSubdomainBenchmark();
            ISolver skylineSolver = DefineSolver(benchmarkSub1, SolverType.Skyline);
            //Console.WriteLine("Uncracked analysis, 1 subdomain,  Skyline : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub1.Model, skylineSolver));
            //Console.WriteLine("Cracked analysis only 1 step, 1 subdomain,  Skyline : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub1.Model, benchmarkSub1.Crack, skylineSolver));
            RunCrackPropagationAnalysis(benchmarkSub1, skylineSolver);

            // FETI-1 5 subdomains
            //benchmarkSub5 = CreateMultiSubdomainBenchmark(5);
            //ISolver solverFeti1 = DefineSolver(benchmarkSub5, SolverType.Feti1);
            //PlotSubdomains(benchmarkSub5, solverFeti1);
            //Console.WriteLine("Uncracked analysis, 5 subdomains, FETI-1  : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub5.Model, solverFeti1));
            //Console.WriteLine("Cracked analysis only 1 step, 5 subdomains, FETI-1  : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub5.Model, benchmarkSub5.Crack, solverFeti1));
            //RunCrackPropagationAnalysis(benchmarkSub5, solverFeti1);

            // FETI-1 7 subdomains
            benchmarkSub7 = CreateMultiSubdomainBenchmark(7);
            ISolver solverFeti1 = DefineSolver(benchmarkSub7, SolverType.Feti1);
            //PlotSubdomains(benchmarkSub7, solverFeti1);
            //Console.WriteLine("Uncracked analysis, 7 subdomains, FETI-1  : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub7.Model, solverFeti1));
            //Console.WriteLine("Cracked analysis only 1 step, 5 subdomains, FETI-1  : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub7.Model, benchmarkSub7.Crack, solverFeti1));
            RunCrackPropagationAnalysis(benchmarkSub7, solverFeti1);

            // FETI-DP 5 subdomains
            //benchmarkSub5 = CreateMultiSubdomainBenchmark(5);
            //ISolver solverFetiDP = DefineSolver(benchmarkSub5, SolverType.FetiDP);
            //PlotSubdomains(benchmarkSub5, solverFetiDP);
            //Console.WriteLine("Uncracked analysis, 5 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub5.Model, solverFetiDP));
            //Console.WriteLine("Cracked analysis only 1 step, 5 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub5.Model, benchmarkSub5.Crack, solverFetiDP));
            //RunCrackPropagationAnalysis(benchmarkSub5, solverFetiDP);

            // FETI-DP 7 subdomains
            benchmarkSub7 = CreateMultiSubdomainBenchmark(7);
            ISolver solverFetiDP = DefineSolver(benchmarkSub7, SolverType.FetiDP);
            //PlotSubdomains(benchmarkSub7, solverFetiDP);
            //Console.WriteLine("Uncracked analysis, 7 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub7.Model, solverFetiDP));
            //Console.WriteLine("Cracked analysis only 1 step, 7 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub7.Model, benchmarkSub7.Crack, solverFetiDP));
            RunCrackPropagationAnalysis(benchmarkSub7, solverFetiDP);

            Console.Write("\nEnd");
        }

        private static FilletBenchmark CreateMultiSubdomainBenchmark(int numSubdomains)
        {
            // Define subdomain boundaries
            double tol = 1E-13;
            var regions = new Dictionary<int, RectangularRegion2D>();
            if (numSubdomains == 5)
            { // The crack will grow in only one subdomain
                // Uncracked: 3 at bottom and 1 at top
                regions[0] = new RectangularRegion2D(0.0, 0.0, 100.0, 75.0, tol);
                regions[0].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[1] = new RectangularRegion2D(275.0, 0.0, 375.0, 75.0, tol);
                regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[2] = new RectangularRegion2D(100.0, 0.0, 275.0, 75.0, tol);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[3] = new RectangularRegion2D(150.0, 110.0, 225.0, 150.0, tol);
                regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);

                // Cracked: 1 in the centre
                regions[4] = new RectangularRegion2D(100.0, 75.0, 275.0, 110.0, tol);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
            }
            else if (numSubdomains == 7)
            { // The uncracked subdomains are the same as before, but now the crack will grow in 3 subdomains.
                // Uncracked: 3 at bottom and 1 at top
                regions[0] = new RectangularRegion2D(0.0, 0.0, 100.0, 75.0, tol);
                regions[0].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[1] = new RectangularRegion2D(275.0, 0.0, 375.0, 75.0, tol);
                regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[2] = new RectangularRegion2D(100.0, 0.0, 275.0, 75.0, tol);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[3] = new RectangularRegion2D(150.0, 110.0, 225.0, 150.0, tol);
                regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);

                // Cracked: 3 in the centre
                double verticalBoundary0 = 172;
                double verticalBoundary1 = 205;
                regions[4] = new RectangularRegion2D(100.0, 75.0, verticalBoundary0, 110.0, tol);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[5] = new RectangularRegion2D(verticalBoundary0, 75.0, verticalBoundary1, 110.0, tol);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[6] = new RectangularRegion2D(verticalBoundary1, 75.0, 275.0, 110.0, tol);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            }
            else
            {
                throw new NotImplementedException();
            }

            // Create the model, crack & mesh with only 1 subdomain
            FilletBenchmark benchmark = CreateSingleSubdomainBenchmark();

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

        private static FilletBenchmark CreateSingleSubdomainBenchmark()
        {
            double growthLength = 5; // mm. Must be sufficiently larger than the element size.

            var builder = new FilletBenchmark.Builder(meshPath, growthLength);
            builder.LsmPlotDirectory = crackPlotDirectory;
            builder.SubdomainPlotDirectory = subdomainPlotDirectory;
            builder.HeavisideEnrichmentTolerance = 0.01;
            builder.RigidBCs = true;
            builder.MaxIterations = 13;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            FilletBenchmark benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static ISolver DefineSolver(FilletBenchmark benchmark, SolverType solverType)
        {
            if (solverType == SolverType.Skyline)
            {
                var builder = new SkylineSolver.Builder();
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.Feti1)
            {
                var factorizationTolerances = new Dictionary<int, double>();
                if (benchmark.Model.Subdomains.Count == 5)
                {
                    factorizationTolerances[0] = 1E-2;
                    factorizationTolerances[1] = 1E-2;
                    factorizationTolerances[2] = 1E-2;
                    factorizationTolerances[3] = 1E-2;
                    factorizationTolerances[4] = 1E1;
                }
                else if (benchmark.Model.Subdomains.Count == 7)
                {
                    factorizationTolerances[0] = 1E-2;
                    factorizationTolerances[1] = 1E-2;
                    factorizationTolerances[2] = 1E-2;
                    factorizationTolerances[3] = 1E-2;
                    factorizationTolerances[4] = 1E1;
                    factorizationTolerances[5] = 1E1;
                    factorizationTolerances[6] = 1E1;
                }
                else
                {
                    throw new NotImplementedException();
                }
                //var fetiMatrices = new DenseFeti1SubdomainMatrixManager.Factory();
                var fetiMatrices = new SkylineFeti1SubdomainMatrixManager.Factory();
                var builder = new Feti1Solver.Builder(fetiMatrices, factorizationTolerances);
                //builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                builder.PreconditionerFactory = new DirichletPreconditioner.Factory();
                builder.ProblemIsHomogeneous = true;
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.FetiDP)
            {
                benchmark.Model.ConnectDataStructures();
                Dictionary<int, HashSet<INode>> cornerNodes = null;
                
                if (benchmark.Model.Subdomains.Count == 5)
                {
                    //// This does not work for this mesh, as there are no crosspoints
                    ////cornerNodes = FindCornerNodesFromCrosspoints2D(benchmark.Model);

                    // The bottom and right subdomains do not need corner nodes, as they are fully fixed.
                    cornerNodes = new Dictionary<int, HashSet<INode>>();
                    cornerNodes[0] = new HashSet<INode>();
                    cornerNodes[1] = new HashSet<INode>();

                    // Find the 4 corners of the cracked subdomain
                    XSubdomain crackedSubdomain = benchmark.Model.Subdomains[4];
                    double tol = 1E-6;

                    INode[] leftSideNodes = crackedSubdomain.Nodes.Where(node => Math.Abs(node.X - 150.0) <= tol).ToArray();
                    INode cornerCrackedTopLeft = leftSideNodes[0];
                    foreach (INode node in leftSideNodes)
                    {
                        if (node.Y > cornerCrackedTopLeft.Y) cornerCrackedTopLeft = node;
                    }

                    INode[] rightSideNodes = crackedSubdomain.Nodes.Where(node => Math.Abs(node.X - 225.0) <= tol).ToArray();
                    INode cornerCrackedTopRight = rightSideNodes[0];
                    foreach (INode node in rightSideNodes)
                    {
                        if (node.Y > cornerCrackedTopRight.Y) cornerCrackedTopRight = node;
                    }

                    INode cornerCrackedBottomLeft = crackedSubdomain.Nodes[0];
                    foreach (INode node in crackedSubdomain.Nodes)
                    {
                        if (node.X < cornerCrackedBottomLeft.X) cornerCrackedBottomLeft = node;
                    }

                    INode cornerCrackedBottomRight = crackedSubdomain.Nodes[0];
                    foreach (INode node in crackedSubdomain.Nodes)
                    {
                        if (node.X > cornerCrackedBottomRight.X) cornerCrackedBottomRight = node;
                    }

                    // Fill the rest of the corner nodes
                    cornerNodes[2] = new HashSet<INode>(new INode[] { cornerCrackedBottomLeft, cornerCrackedBottomRight });
                    cornerNodes[3] = new HashSet<INode>(new INode[] { cornerCrackedTopLeft, cornerCrackedTopRight });
                    cornerNodes[4] = new HashSet<INode>(new INode[]
                    {
                        cornerCrackedBottomLeft, cornerCrackedBottomRight, cornerCrackedTopLeft, cornerCrackedTopRight
                    });
                }
                else if (benchmark.Model.Subdomains.Count == 7)
                {
                    cornerNodes = FindCornerNodesFromCrosspoints2D(benchmark.Model);
                }

                //var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes);
                var cornerNodeSelection = new CrackedFetiDPSubdomainCornerNodes(benchmark.Crack, cornerNodes);
                //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
                var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
                var builder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
                //builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                builder.PreconditionerFactory = new DirichletPreconditioner.Factory();
                builder.ProblemIsHomogeneous = true;
                return builder.BuildSolver(benchmark.Model);
            }
            else throw new ArgumentException("Invalid solver choice.");
        }

        private static void PlotSubdomains(FilletBenchmark benchmark, ISolver solver)
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

        private static void RunCrackPropagationAnalysis(FilletBenchmark benchmark, ISolver solver)
        {
            benchmark.Analyze(solver);

            // Write crack path
            Console.WriteLine("Crack path:");
            foreach (var point in benchmark.Crack.CrackPath)
            {
                Console.WriteLine($"{point.X} {point.Y}");
            }
            Console.WriteLine();

            solver.Logger.WriteToFile(solverLogPath, $"{solver.Name}_log", true);
            solver.Logger.WriteAggregatesToFile(solverLogPath, $"{solver.Name}_log", true);
        }
    }
}
