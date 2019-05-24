using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Tests;

namespace ISAAR.MSolve.SamplesConsole.XFEM.COMPDYN2019
{
    public class Fillet
    {
        private enum SolverType { Skyline, Feti1, FetiDP}

        public static void Run()
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Fillet\Mesh\fillet.msh";
            string subdomainPlotPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Fillet\Plots\subdomains.vtk";

            FilletBenchmark benchmarkSub1 = CreateSingleSubdomainBenchmark(meshPath);
            ISolver skylineSolver = DefineSolver(benchmarkSub1, SolverType.Skyline);
            double uNormSkyline = RunUncrackedAnalysis(benchmarkSub1, skylineSolver);

            FilletBenchmark benchmarkSub5 = CreateMultiSubdomainBenchmark(5, meshPath);
            PlotSubdomains(subdomainPlotPath, benchmarkSub5.Model);
            ISolver solverFeti = DefineSolver(benchmarkSub5, SolverType.Feti1);
            double uNormFeti = RunUncrackedAnalysis(benchmarkSub5, solverFeti);
        }

        private static FilletBenchmark CreateMultiSubdomainBenchmark(int numSubdomains, string meshPath)
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
                regions[4] = new RectangularRegion2D(100.0, 75.0, 158.333, 110.0, tol);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[5] = new RectangularRegion2D(158.333, 75.0, 216.667, 110.0, tol);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[6] = new RectangularRegion2D(216.667, 75.0, 275.0, 110.0, tol);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            }
            else
            {
                throw new NotImplementedException();
            }

            // Create the model, crack & mesh with only 1 subdomain
            FilletBenchmark benchmark = CreateSingleSubdomainBenchmark(meshPath);

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

        private static FilletBenchmark CreateSingleSubdomainBenchmark(string meshPath)
        {
            double growthLength = 5; // mm. Must be sufficiently larger than the element size.

            var builder = new FilletBenchmark.Builder(meshPath, growthLength);
            builder.HeavisideEnrichmentTolerance = 0.01;
            builder.RigidBCs = true;
            builder.NumSubdomains = 5;
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

                }
                else
                {
                    throw new NotImplementedException();
                }
                var builder = new Feti1Solver.Builder(factorizationTolerances);
                builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                builder.ProblemIsHomogeneous = true;
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.FetiDP)
            {
                // Must also specify corner nodes
                throw new NotImplementedException();
            }
            else throw new ArgumentException("Invalid solver choice.");
        }

        private static void RunCrackPropagationAnalysis(FilletBenchmark benchmark, ISolver solver)
        {
        }

        private static double RunUncrackedAnalysis(FilletBenchmark benchmark, ISolver solver)
        {
            var problem = new ProblemStructural(benchmark.Model, solver);
            var linearAnalyzer = new LinearAnalyzer(benchmark.Model, solver, problem);
            var staticAnalyzer = new StaticAnalyzer(benchmark.Model, solver, problem, linearAnalyzer);

            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Return norm2(displacements)
            if (benchmark.Model.Subdomains.Count == 1)
            {
                return solver.LinearSystems.First().Value.Solution.Norm2();
            }
            else if (solver is Feti1Solver fetiSolver)
            {
                var sudomainDisplacements = new Dictionary<int, IVectorView>();
                foreach (var ls in fetiSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
                return fetiSolver.GatherGlobalDisplacements(sudomainDisplacements).Norm2();
            }
            else throw new NotImplementedException("Invalid solver");
        }

        private static void PlotSubdomains(string plotPath, XModel model)
        {
            model.ConnectDataStructures();
            var writer = new VtkMeshPartitionWriter();
            var nodesPerSubdomain = new Dictionary<int, IReadOnlyList<XNode>>();
            var elementsPerSubdomain = new Dictionary<int, IReadOnlyList<IXFiniteElement>>();
            foreach (int subdomainID in model.Subdomains.Keys)
            {
                nodesPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Nodes;
                elementsPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Elements;
            }
            writer.WriteSubdomainElements(plotPath, nodesPerSubdomain, elementsPerSubdomain);
        }
    }
}
