using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class MinresSolver : SolverBase
    {
        private readonly double maxIterationsOverOrder; //roughly
        private readonly double tolerance;

        public MinresSolver(Model2D model, double maxIterationsOverOrder, double tolerance) : base(model)
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
            //if (!Kuu.IsSymmetric(1e-10)) throw new AsymmetricMatrixException(
            //    "Stiffness matrix corresponding to free-free dofs is not symmetric");

            // Without preconditioning
            int maxIterations = (int)Math.Ceiling(Kuu.NumColumns * maxIterationsOverOrder);
            var minres = new MinimumResidual(maxIterations, tolerance, false, false);
            (Vector x, MinresStatistics statistics) = minres.Solve(Kuu.BuildCSRMatrix(true), rhs);

            // With preconditioning
            //var minres = new PreconditionedMinimumResidual(maxIterations, tolerance, false, false);
            //(double[] diagonal, int firstZeroIdx) = Kuu.GetDiagonalAsArray(); // Preconditioner could be abstracted, but I think it depends on the solver.
            //var preconditioner = new JacobiPreconditioner(diagonal);
            //var preconditioner = new IdentityPreconditioner(true);
            //(Vector x, MinresStatistics statistics) = minres.Solve(Kuu.BuildCSRMatrix(true), rhs, preconditioner);

            #region Debugging
            //var writer = new MatlabWriter();
            //string directory = @"C:\Users\Serafeim\Desktop\GRACM\Minres_debugging\";
            //writer.WriteSparseMatrix(Kuu, directory + "matrix.txt");
            //writer.WriteFullVector(rhs, directory + "rhs.txt");
            //writer.WriteFullVector(x, directory + "solution.txt");
            #endregion

            Console.WriteLine(statistics);
            Solution = x;

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }
    }
}
