using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Reordering;
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
        private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_8000.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_13738.msh";
        ////private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_22666.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_29052.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_60000.msh";
        //private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_100000.msh";
        ////private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_189740.msh";
        ////private const string meshPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Mesh\holes_357324.msh";
        private const string leftCrackPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Plots\Left";
        private const string rightCrackPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Plots\Right";
        private const string subdomainPlotDirectory = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\Plots\Subdomains";
        private const string solverLogPath = @"C:\Users\Serafeim\Desktop\COMPDYN2019\Holes\solver_log.txt";

        public static void Run()
        {
            HolesBenchmark benchmarkSub1;
            HolesBenchmark benchmarkSub10;
            HolesBenchmark benchmarkSub15;

            // Skyline
            //benchmarkSub1 = CreateSingleSubdomainBenchmark();
            //ISolver skylineSolver = DefineSolver(benchmarkSub1, SolverType.Skyline);
            //RunCrackPropagationAnalysis(benchmarkSub1, skylineSolver);
            //Console.WriteLine("Uncracked analysis, 1 subdomain, Skyline   : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub1.Model, skylineSolver));
            //Console.WriteLine("Cracked analysis only 1 step, 1 subdomain, Skyline   : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub1.Model, benchmarkSub1.Crack, skylineSolver));
            //Console.WriteLine("Skyline solver, 1 subdomain: ");

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
            RunCrackPropagationAnalysis(benchmarkSub10, solverFetiDP);
            //PlotSubdomains(benchmarkSub10, solverFetiDP);
            //Console.WriteLine("Uncracked analysis, 10 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunUncrackedAnalysis(benchmarkSub10.Model, solverFetiDP));
            //Console.WriteLine("Cracked analysis only 1 step, 10 subdomains, FETI-DP : norm2(globalU) = " +
            //    RunSingleCrackedStep(benchmarkSub10.Model, benchmarkSub10.Crack, solverFetiDP));
            //Console.WriteLine("FETI-DP, 10 subdomains: ");

            // FETI-DP 15 subdomains
            //benchmarkSub15 = CreateMultiSubdomainBenchmark(15);
            //ISolver solverFetiDP = DefineSolver(benchmarkSub15, SolverType.FetiDP);
            //RunCrackPropagationAnalysis(benchmarkSub15, solverFetiDP);

            Console.Write("\nEnd");
        }

        private static HolesBenchmark CreateMultiSubdomainBenchmark(int numSubdomains)
        {
            // Define subdomain boundaries
            double tol = 1E-13;
            var regions = new Dictionary<int, IRegion2D>();
            if (numSubdomains == 10)
            {
                double xMin = 0.0, yMin = 0.0, xMax = 20.0, yMax = 10.0;
                double x1 = 2.75, x2 = 3.5, x3 = 5.0, x4 = 8.0, x5 = 12.2, x6 = 15.5, x7 = 16.5, x8 = 16.75;
                double y1 = 1.85, y2 = 5.0, y3 = 8.85;

                // Left crack regions: 
                var region0 = new RectangularRegion2D(xMin, yMin, x2, y2, tol);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[0] = region0;

                var region1 = new RectangularRegion2D(x2, yMin, x4, y2, tol);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[1] = region1;

                var region2 = new RectangularRegion2D(x4, yMin, x5, y2, tol);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[2] = region2;

                var region3Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x5, yMin), new CartesianPoint(x6, yMin), new CartesianPoint(x6, y1),
                    new CartesianPoint(x8, y2), new CartesianPoint(x5, y2)
                };
                var region3Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region3Vertices[1], region3Vertices[2]),
                    new LineSegment2D(region3Vertices[2], region3Vertices[3]),
                    new LineSegment2D(region3Vertices[3], region3Vertices[4]),
                    new LineSegment2D(region3Vertices[4], region3Vertices[1])
                };
                regions[3] = new PolygonalRegion2D(region3Vertices, region3Boundaries);

                var region4Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x6, yMin), new CartesianPoint(xMax, yMin), new CartesianPoint(xMax, y2),
                    new CartesianPoint(x8, y2), new CartesianPoint(x6, y1)
                };
                var region4Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region4Vertices[2], region4Vertices[3]),
                    new LineSegment2D(region4Vertices[3], region4Vertices[4]),
                    new LineSegment2D(region4Vertices[4], region4Vertices[0])
                };
                regions[4] = new PolygonalRegion2D(region4Vertices, region4Boundaries);

                // Right crack regions:
                var region5Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(xMin, y2), new CartesianPoint(x1, y2), new CartesianPoint(x3, y3),
                    new CartesianPoint(x3, yMax), new CartesianPoint(xMin, yMax)
                };
                var region5Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region5Vertices[0], region5Vertices[1]),
                    new LineSegment2D(region5Vertices[1], region5Vertices[2]),
                    new LineSegment2D(region5Vertices[2], region5Vertices[3]),
                };
                regions[5] = new PolygonalRegion2D(region5Vertices, region5Boundaries);

                var region6Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x1, y2), new CartesianPoint(x4, y2), new CartesianPoint(x4, yMax),
                    new CartesianPoint(x3, yMax), new CartesianPoint(x3, y3)
                };
                var region6Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region6Vertices[0], region6Vertices[1]),
                    new LineSegment2D(region6Vertices[1], region6Vertices[2]),
                    new LineSegment2D(region6Vertices[3], region6Vertices[4]),
                    new LineSegment2D(region6Vertices[4], region6Vertices[0]),
                };
                regions[6] = new PolygonalRegion2D(region6Vertices, region6Boundaries);

                var region7 = new RectangularRegion2D(x4, y2, x5, yMax, tol);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[7] = region7;

                var region8 = new RectangularRegion2D(x5, y2, x7, yMax, tol);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[8] = region8;

                var region9 = new RectangularRegion2D(x7, y2, xMax, yMax, tol);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[9] = region9;
            }
            else if (numSubdomains == 15)
            {
                double xMin = 0.0, yMin = 0.0, xMax = 20.0, yMax = 10.0;
                double x1 = 2.75, x2 = 3.5, x3 = 5.0, x4 = 8.0, x5 = 12.2, x6 = 15.5, x7 = 16.5, x8 = 16.75;
                double y1 = 1.85, y2 = 3.25, y3 = 6.5, y4 = 8.85;

                // Bottom row: 
                var region0 = new RectangularRegion2D(xMin, yMin, x2, y2, tol);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region0.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[0] = region0;

                var region1 = new RectangularRegion2D(x2, yMin, x4, y2, tol);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region1.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[1] = region1;

                var region2 = new RectangularRegion2D(x4, yMin, x5, y2, tol);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region2.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[2] = region2;

                var region3Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x5, yMin), new CartesianPoint(x6, yMin), new CartesianPoint(x6, y1),
                    new CartesianPoint(x8, y2), new CartesianPoint(x5, y2)
                };
                var region3Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region3Vertices[1], region3Vertices[2]),
                    new LineSegment2D(region3Vertices[2], region3Vertices[3]),
                    new LineSegment2D(region3Vertices[3], region3Vertices[4]),
                    new LineSegment2D(region3Vertices[4], region3Vertices[1])
                };
                regions[3] = new PolygonalRegion2D(region3Vertices, region3Boundaries);

                var region4Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x6, yMin), new CartesianPoint(xMax, yMin), new CartesianPoint(xMax, y2),
                    new CartesianPoint(x8, y2), new CartesianPoint(x6, y1)
                };
                var region4Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region4Vertices[2], region4Vertices[3]),
                    new LineSegment2D(region4Vertices[3], region4Vertices[4]),
                    new LineSegment2D(region4Vertices[4], region4Vertices[0])
                };
                regions[4] = new PolygonalRegion2D(region4Vertices, region4Boundaries);

                // Middle row:
                var region5 = new RectangularRegion2D(xMin, y2, x2, y3, tol);
                region5.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region5.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region5.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[5] = region5;

                var region6 = new RectangularRegion2D(x2, y2, x4, y3, tol);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region6.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[6] = region6;

                var region7 = new RectangularRegion2D(x4, y2, x5, y3, tol);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region7.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[7] = region7;

                var region8 = new RectangularRegion2D(x5, y2, x7, y3, tol);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region8.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[8] = region8;

                var region9 = new RectangularRegion2D(x7, y2, xMax, y3, tol);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region9.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[9] = region9;

                // Top row:
                var region10Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(xMin, y3), new CartesianPoint(x1, y3), new CartesianPoint(x3, y4),
                    new CartesianPoint(x3, yMax), new CartesianPoint(xMin, yMax)
                };
                var region10Boundaries = new LineSegment2D[3]
                {
                    new LineSegment2D(region10Vertices[0], region10Vertices[1]),
                    new LineSegment2D(region10Vertices[1], region10Vertices[2]),
                    new LineSegment2D(region10Vertices[2], region10Vertices[3]),
                };
                regions[10] = new PolygonalRegion2D(region10Vertices, region10Boundaries);

                var region11Vertices = new CartesianPoint[]
                {
                    new CartesianPoint(x1, y3), new CartesianPoint(x4, y3), new CartesianPoint(x4, yMax),
                    new CartesianPoint(x3, yMax), new CartesianPoint(x3, y4)
                };
                var region11Boundaries = new LineSegment2D[4]
                {
                    new LineSegment2D(region11Vertices[0], region11Vertices[1]),
                    new LineSegment2D(region11Vertices[1], region11Vertices[2]),
                    new LineSegment2D(region11Vertices[3], region11Vertices[4]),
                    new LineSegment2D(region11Vertices[4], region11Vertices[0]),
                };
                regions[11] = new PolygonalRegion2D(region11Vertices, region11Boundaries);

                var region12 = new RectangularRegion2D(x4, y3, x5, yMax, tol);
                region12.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region12.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region12.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[12] = region12;

                var region13 = new RectangularRegion2D(x5, y3, x7, yMax, tol);
                region13.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region13.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                region13.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
                regions[13] = region13;

                var region14 = new RectangularRegion2D(x7, y3, xMax, yMax, tol);
                region14.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
                region14.AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
                regions[14] = region14;
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
                benchmark.Partitioner = new TipAdaptivePartitioner(benchmark.Crack);
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
                //var fetiMatrices = new DenseFeti1SubdomainMatrixManager.Factory();
                var fetiMatrices = new SkylineFeti1SubdomainMatrixManager.Factory(new OrderingAmdSuiteSparse()); 
                var builder = new Feti1Solver.Builder(fetiMatrices, factorizationTolerances);
                //builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                builder.PreconditionerFactory = new DirichletPreconditioner.Factory();
                builder.ProblemIsHomogeneous = true;
                return builder.BuildSolver(benchmark.Model);
            }
            else if (solverType == SolverType.FetiDP)
            {
                //benchmark.Partitioner = new TipAdaptivePartitioner(benchmark.Crack);
                Dictionary<int, HashSet<INode>> cornerNodes = null;

                if (benchmark.Model.Subdomains.Count == 10 || benchmark.Model.Subdomains.Count == 15)
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
                //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
                //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
                var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory(new OrderingAmdSuiteSparse());
                var builder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
                //builder.PreconditionerFactory = new LumpedPreconditioner.Factory();
                //builder.PreconditionerFactory = new DiagonalDirichletPreconditioner.Factory();
                builder.PreconditionerFactory = new DirichletPreconditioner.Factory();
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
