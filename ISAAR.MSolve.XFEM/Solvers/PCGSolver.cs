using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class PCGSolver: SolverBase
    {
        private readonly double maxIterationsOverOrder; //roughly
        private readonly double tolerance;

        public PCGSolver(Model2D model, double maxIterationsOverOrder, double tolerance) : base(model)
        {
            this.maxIterationsOverOrder = maxIterationsOverOrder;
            this.tolerance = tolerance;
        }

        public override void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Interleaced and separate dof enumerators seem to have similar performance.
            DofOrderer = InterleavedDofOrderer.Create(model);
            //DofOrderer = DofOrdererSeparate.Create(model);

            var assembler = new GlobalCSRAssembler();
            (DOKRowMajor Kuu, DOKRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);


            int maxIterations = (int)Math.Ceiling(Kuu.NumColumns * maxIterationsOverOrder);
            var pcg = new PreconditionedConjugateGradient(maxIterations, tolerance);

            // Preconditioner could be abstracted, but I think it depends on the solver.
            (double[] diagonal, int firstZeroIdx) = Kuu.GetDiagonalAsArray();
            var preconditioner = new JacobiPreconditioner(diagonal);
            //var preconditioner = new IdentityPreconditioner(true);

            (Vector x, IterativeStatistics statistics) = pcg.Solve(Kuu.BuildCSRMatrix(true), rhs, preconditioner);
            Solution = x;

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }
    }
}
