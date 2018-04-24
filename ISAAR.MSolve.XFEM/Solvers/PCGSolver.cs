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
    class PCGSolver: ISolver
    {
        private readonly Model2D model;
        private readonly double maxIterationsOverOrder; //roughly
        private readonly double tolerance;

        public PCGSolver(Model2D model, double maxIterationsOverOrder, double tolerance)
        {
            this.model = model;
            this.maxIterationsOverOrder = maxIterationsOverOrder;
            this.tolerance = tolerance;
        }

        public IDOFEnumerator DOFEnumerator { get; private set; }

        public Vector Solution { get; private set; }

        public void Initialize() { }

        public void Solve()
        {
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            (DOKRowMajor matrix, Vector rhs) = ReduceToSimpleLinearSystem();

            int maxIterations = (int)Math.Ceiling(matrix.NumColumns * maxIterationsOverOrder);
            var pcg = new PCGAlgorithm(maxIterations, tolerance);

            // Preconditioner could be abstracted, but I think it dpends on the solver.
            (double[] diagonal, int firstZeroIdx) = matrix.GetDiagonalAsArray();
            var preconditioner = new JacobiPreconditioner(diagonal);
            //var preconditioner = new IdentityPreconditioner(true);

            (Vector x, IterativeStatistics statistics) = pcg.Solve(matrix.BuildCSRMatrix(true), rhs, preconditioner);
            Solution = x;
        }

        /// <summary>
        /// The extended linear system is:
        /// [Kcc Kcu; Kuc Kuu] * [uc; uu] = [Fc; Fu]
        /// where c are the standard constrained dofs, f are the standard free dofs, e are the enriched dofs and 
        /// u = Union(f,c) are both the dofs with unknown left hand side vectors: uu = [uf; ue].
        /// To solve the system (for the unknowns ul):
        /// i) Kuu * uu = Fu - Kuc * uc = Feff
        /// ii) uu = Kuu \ Feff 
        /// </summary>
        /// <returns></returns>
        private (DOKRowMajor matrix, Vector rhs) ReduceToSimpleLinearSystem()
        {
            var assembler = new GlobalCSRAssembler();
            (DOKRowMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);

            // TODO: Perhaps a dedicated class should be responsible for these vectors
            Vector Fu = model.CalculateFreeForces(DOFEnumerator);
            Vector uc = model.CalculateConstrainedDisplacements(DOFEnumerator);
            Vector Feff = Fu - Kuc * uc;
            return (Kuu, Feff);
        }
    }
}
