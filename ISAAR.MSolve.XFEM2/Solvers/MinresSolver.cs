using System;
using ISAAR.MSolve.LinearAlgebra.Iterative.MinimumResidual;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class MinresSolver : SolverBase
    {
        private readonly double maxIterationsOverOrder; //roughly
        private readonly int reorthoSize;
        private readonly double tolerance;

        public MinresSolver(Model2D model, double maxIterationsOverOrder, double tolerance, int reorthoSize) : base(model)
        {
            this.maxIterationsOverOrder = maxIterationsOverOrder;
            this.tolerance = tolerance;
            this.reorthoSize = reorthoSize;
        }

        public override void Solve()
        {
            // Interleaced and separate dof enumerators seem to have similar performance.
            DofOrderer = InterleavedDofOrderer.Create(model);
            //DofOrderer = DofOrdererSeparate.Create(model);

            var assembler = new GlobalCSRAssembler();
            (DokRowMajor Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);
            //if (!Kuu.IsSymmetric(1e-10)) throw new AsymmetricMatrixException(
            //    "Stiffness matrix corresponding to free-free dofs is not symmetric");

            // Without preconditioning
            int maxIterations = (int)Math.Ceiling(Kuu.NumColumns * maxIterationsOverOrder);
            var minres = new MinRes(maxIterations, tolerance, reorthoSize, false, false);
            (IVector x, MinresStatistics statistics) = minres.Solve(Kuu.BuildCsrMatrix(true), rhs);

            // With preconditioning
            //var minres = new PreconditionedMinimumResidual(maxIterations, tolerance, false, false);
            //(double[] diagonal, int firstZeroIdx) = Kuu.GetDiagonalAsArray(); // Preconditioner could be abstracted, but I think it depends on the solver.
            //var preconditioner = new JacobiPreconditioner(diagonal);
            //var preconditioner = new IdentityPreconditioner(true);
            //(Vector x, MinresStatistics statistics) = minres.Solve(Kuu.BuildCSRMatrix(true), rhs, preconditioner);

            #region Debugging
            double relres = rhs.Subtract(Kuu.BuildCsrMatrix(true).Multiply(x)).Norm2() / rhs.Norm2();
            //var writer = new MatlabWriter();
            //string directory = @"C:\Users\Serafeim\Desktop\GRACM\Minres_debugging\";
            //writer.WriteSparseMatrix(Kuu, directory + "matrix.txt");
            //writer.WriteFullVector(rhs, directory + "rhs.txt");
            //writer.WriteFullVector(x, directory + "solution.txt");
            #endregion

            Console.WriteLine(statistics);
            Solution = x;
        }
    }
}
