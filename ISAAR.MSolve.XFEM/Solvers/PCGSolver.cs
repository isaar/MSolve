using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Solvers.Algorithms;
using ISAAR.MSolve.XFEM.Solvers.Preconditioning;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class PCGSolver: SolverBase
    {
        private readonly double maxIterationsOverOrder; //roughly
        private readonly double tolerance;

        public PCGSolver(Model2D model, double maxIterationsOverOrder, double tolerance): base(model)
        {
            this.maxIterationsOverOrder = maxIterationsOverOrder;
            this.tolerance = tolerance;
        }

        public override void Solve()
        {
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            var assembler = new GlobalCSRAssembler();
            (DOKRowMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);


            int maxIterations = (int)Math.Ceiling(Kuu.NumColumns * maxIterationsOverOrder);
            var pcg = new PCGAlgorithm(maxIterations, tolerance);

            // Preconditioner could be abstracted, but I think it depends on the solver.
            (double[] diagonal, int firstZeroIdx) = Kuu.GetDiagonalAsArray();
            var preconditioner = new JacobiPreconditioner(diagonal);
            //var preconditioner = new IdentityPreconditioner(true);

            (Vector x, IterativeStatistics statistics) = pcg.Solve(Kuu.BuildCSRMatrix(true), rhs, preconditioner);
            Solution = x;
        }
    }
}
