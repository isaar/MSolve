using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class PCGSolver: ISolver
    {
        private readonly Model2D model;
        private readonly double maxIterationsOverOrder; //roughly
        private readonly double tolerance;
        private int iteration;

        public PCGSolver(Model2D model, double maxIterationsOverOrder, double tolerance)
        {
            this.model = model;
            this.maxIterationsOverOrder = maxIterationsOverOrder;
            this.tolerance = tolerance;
            Logger = new SolverLogger("PcgSolver");

        }
        public IDofOrderer DofOrderer { get; protected set; }

        public SolverLogger Logger { get; }

        public Vector Solution { get; protected set; }

        public void Initialize()
        {
            iteration = 0;
        }

        public void Solve()
        {
            ++iteration;
            var watch = new Stopwatch();
            watch.Start();

            // Interleaved and separate dof enumerators seem to have similar performance.
            DofOrderer = InterleavedDofOrderer.Create(model);
            //DofOrderer = DofOrdererSeparate.Create(model);

            var assembler = new GlobalCSRAssembler();
            (DokRowMajor Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            //if (!Kuu.IsSymmetric(1e-10)) throw new AsymmetricMatrixException(
            //    "Stiffness matrix corresponding to free-free dofs is not symmetric");
            Vector rhs = CalcEffectiveRhs(Kuc);

            watch.Stop();
            Logger.LogDuration(iteration, "Linear system assembly", watch.ElapsedMilliseconds);

            // Preconditioner
            watch.Restart();
            (double[] diagonal, int firstZeroIdx) = Kuu.GetDiagonalAsArray(); // Preconditioner could be abstracted, but I think it depends on the solver.
            var preconditioner = new JacobiPreconditioner(diagonal);
            //var preconditioner = new IdentityPreconditioner(true);
            watch.Stop();
            Logger.LogDuration(iteration, "Building the preconditioner", watch.ElapsedMilliseconds);

            // PCG
            watch.Restart();
            int maxIterations = (int)Math.Ceiling(Kuu.NumColumns * maxIterationsOverOrder);
            var pcg = new PreconditionedConjugateGradient(maxIterations, tolerance);
            (Vector x, CGStatistics statistics) = pcg.Solve(Kuu.BuildCsrMatrix(true), rhs, preconditioner);
            //var cg = new ConjugateGradient(maxIterations, tolerance);
            //(Vector x, IterativeStatistics statistics) = cg.Solve(Kuu.BuildCSRMatrix(true), rhs);
            watch.Stop();
            Logger.LogDuration(iteration, "PCG", watch.ElapsedMilliseconds);

            Console.WriteLine(statistics);
            Solution = x;

            Logger.LogDofs(iteration, DofOrderer.NumStandardDofs + DofOrderer.NumEnrichedDofs);
            #region Debugging
            //CheckPCG(model, DofOrderer, Kuu, Solution);
            #endregion
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
        private Vector CalcEffectiveRhs(DokRowMajor globalUnconstrainedConstrained)
        {
            Vector Fu = model.CalculateFreeForces(DofOrderer);
            Vector uc = model.CalculateConstrainedDisplacements(DofOrderer);
            Vector Feff = Fu - globalUnconstrainedConstrained.MultiplyRight(uc);
            return Feff;
        }

        private static void CheckPCG(Model2D model, IDofOrderer dofOrderer, DokRowMajor Kuu, Vector solution)
        {
            var assembler = new GlobalDOKAssembler();
            (DokSymmetric KuuChol, DokRowMajor KucChol) = assembler.BuildGlobalMatrix(model, dofOrderer);
            if (!KuuChol.Equals(Kuu.BuildCsrMatrix(true), 1e-10)) throw new Exception("Incorrect stiffness matrix assembly");
        }
    }
}
