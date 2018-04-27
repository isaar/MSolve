using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class DCBSolvers
    {
        private delegate ISolver CreateSolver(Model2D model);

        public static void Main()
        {
            int numRepetitions = 10;

            // Define the benchmark problem
            double growthLength = 0.3;
            double elementSize = 0.08;
            var builder = new DCB.Builder(elementSize, growthLength);
            builder.UniformMesh = false;
            builder.UseLSM = true;
            //TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip elements are 
            //      found. It happens at iteration 10.
            builder.MaxIterations = 10;
            builder.KnownPropagation = null; // TODO: enter the fixed propagator here, perhaps by solving the benchmark once.

            // Define the solvers
            string[] solvers = { "Skyline", "Jacobi Preconditioned CG", "SuiteSparse Cholesky" };
            CreateSolver[] callbacks = { CreateSkylineSolver, CreatePCGSolver, CreateCholeskySuiteSparseSolver };
            var totalTimes = new long[solvers.Length];

            // Run the analyses by alternating them to reduce bias.
            for (int t = 0; t < numRepetitions; ++t) 
            {
                Console.WriteLine($"Repetition: {t}");
                for (int s = 0; s < solvers.Length; ++s)
                {
                    DCB benchmark = builder.BuildBenchmark();
                    benchmark.InitializeModel();
                    var solver = callbacks[s](benchmark.Model);
                    benchmark.Analyze(solver);
                    long totalTime = solver.Logger.CalcTotalTime();
                    Console.WriteLine($"Solver {solvers[s]}: total time = {totalTime} ms.");
                    totalTimes[s] += totalTime;
                }
                Console.WriteLine();
            }

            // Print statistics
            for (int s = 0; s < solvers.Length; ++s)
            {
                Console.WriteLine($"Solver {solvers[s]}: Total time = {totalTimes[s] / numRepetitions} ms");
            }
        }

        private static ISolver CreateSkylineSolver(Model2D model)
        {
            return new SkylineSolver(model);
        }

        private static ISolver CreatePCGSolver(Model2D model)
        {
            return new PCGSolver(model, 1, 1e-8);
        }

        private static ISolver CreateCholeskySuiteSparseSolver(Model2D model)
        {
            return new CholeskySuiteSparseSolver(model);
        }
    }
}
