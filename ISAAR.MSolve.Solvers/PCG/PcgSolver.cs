using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using LegacyVector = ISAAR.MSolve.Numerical.LinearAlgebra.Vector;

//TODO: perhaps the user should choose the PCG settings himself and pass it. In this case, this should be named IterativeSolver.
//TODO: the maxIterations of PCG should be able to use the order of the matrix as a default value.
//TODO: IIndexable2D is not a good choice if all solvers must cast it to the matrix types the operate on.
namespace ISAAR.MSolve.Solvers.PCG
{
    /// <summary>
    /// Iterative solver for models with only 1 subdomain. Uses the Proconditioned Conjugate Gradient algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgSolver: ISolver
    {
        private const string name = "PcgSolver"; // for error messages
        private readonly PreconditionedConjugateGradient pcgAlgorithm;
        private readonly IPreconditionerBuilder preconditionerBuilder;
        private readonly LinearSystem_v2<IMatrix, LegacyVector> linearSystem;
        private IPreconditioner preconditioner;

        public PcgSolver(IReadOnlyList<LinearSystem_v2<IMatrix, LegacyVector>> linearSystems, 
            int maxIterations, double residualTolerance, IPreconditionerBuilder preconditionerBuilder)
        {
            if (linearSystems.Count != 1) throw new InvalidSolverException(name + " can be used if there is only 1 subdomain.");
            pcgAlgorithm = new PreconditionedConjugateGradient(maxIterations, residualTolerance);
            this.preconditionerBuilder = preconditionerBuilder;
        }

        public void Initialize()
        {
        }

        /// <summary>
        /// Solves the linear system with PCG method. If the matrix has been modified, a new preconditioner will be computed.
        /// </summary>
        public void Solve()
        {
            if (linearSystem.IsMatrixModified)
            {
                preconditioner = preconditionerBuilder.BuildPreconditioner(linearSystem.Matrix);
                linearSystem.IsMatrixModified = false;
            }
            (Vector solution, CGStatistics stats) = 
                pcgAlgorithm.Solve(linearSystem.Matrix, Vector.CreateFromLegacyVector(linearSystem.RhsVector), preconditioner);
            linearSystem.Solution = solution.ToLegacyVector();
        }
    }
}
