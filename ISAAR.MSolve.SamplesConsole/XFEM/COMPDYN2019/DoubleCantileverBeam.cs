using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests;
using static ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019.Utilities;

namespace ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019
{
    public class DoubleCantileverBeam
    {
        private const int numElementsY = 11;
        private const double tipEnrichementRadius = 0.0;
        private const string crackPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\DCB\Plots";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\DCB\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\DCB\solver_log.txt";

        public static void Run()
        {
            int numSubdomainsX = 1;
            int numSubdomainsY = 1;
            var solverType = SolverType.Skyline;

            DcbBenchmarkBelytschko benchmark = CreateBenchmark(numElementsY, numSubdomainsX, numSubdomainsY, tipEnrichementRadius);
            ISolver solver = DefineSolver(benchmark, solverType);
            RunCrackPropagationAnalysis(benchmark, solver);

            Console.Write("\nEnd");
        }

        private static DcbBenchmarkBelytschko CreateBenchmark(int numElementsY, int numSubdomainsX, int numSubdomainsY, double tipEnrichmentRadius)
        {
            var builder = new DcbBenchmarkBelytschko.Builder(numElementsY, numSubdomainsX, numSubdomainsY);
            builder.LsmPlotDirectory = crackPlotDirectory;
            builder.SubdomainPlotDirectory = subdomainPlotDirectory;
            builder.HeavisideEnrichmentTolerance = 0.01;
            builder.MaxIterations = 13;
            builder.TipEnrichmentRadius = tipEnrichmentRadius;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            DcbBenchmarkBelytschko benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            return benchmark;
        }

        private static ISolver DefineSolver(DcbBenchmarkBelytschko benchmark, SolverType solverType)
        {
            if (solverType == SolverType.Skyline)
            {
                var builder = new SkylineSolver.Builder();
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.Feti1)
            {
                double tol = 1E-7;
                var factorizationTolerances = new Dictionary<int, double>();
                foreach (int s in benchmark.Model.Subdomains.Keys) factorizationTolerances[s] = tol;
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
                Dictionary<int, HashSet<INode>> initialCorners = FindCornerNodesFromCrosspoints2D(benchmark.Model);
                var cornerNodeSelection = new CrackedFetiDPSubdomainCornerNodes(benchmark.Crack, initialCorners);
                var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
                var builder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
                //builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                builder.PreconditionerFactory = new DirichletPreconditioner.Factory();
                builder.ProblemIsHomogeneous = true;
                return builder.BuildSolver(benchmark.Model);
            }
            else throw new ArgumentException("Invalid solver choice.");
        }

        private static void PlotSubdomains(DcbBenchmarkBelytschko benchmark, ISolver solver)
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

        private static void RunCrackPropagationAnalysis(DcbBenchmarkBelytschko benchmark, ISolver solver)
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
