using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Use a dedicated DiagonalMatrix class, instead of passing in double[] or Vector. It will also implement the inverse and 
//      multiplication routines.
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning
{
    /// <summary>
    /// Implements the Jacobi or diagonal preconditioner for a square matrix. If A is the original matrix, the Jacobi  
    /// preconditioner is a matrix M, such that it oncly contains the diagonal of A and inverse(M) is also diagonal with 
    /// entries: 1/A[0,0], 1/A[1,1], ... The Jacobi preconditioner is cheapest to build and apply, but doesn't improve 
    /// convergence as much as other preconditioners.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class JacobiPreconditioner: IPreconditioner
    {
        public const double Tolerance = 1e-10;
        private readonly double[] inverseDiagonal;

        /// <summary>
        /// Initializes a new instance of <see cref="JacobiPreconditioner"/> for the linear system's matrix whose main diagonal
        /// is provided in <paramref name="diagonal"/>.
        /// </summary>
        /// <param name="diagonal">The main diagonal of the original matrix of the linear system. Constraints: 
        ///     all its entries must be non-zero.</param>
        /// <param name="tolerance">The value under which a diagonal entry will be considered as zero.</param>
        /// <exception cref="SingularMatrixException">If there is a zero diagonal entry.</exception>
        public JacobiPreconditioner(double[] diagonal, double tolerance = Tolerance)
        {
            Order = diagonal.Length;
            inverseDiagonal = new double[Order];
            for (int i = 0; i < Order; ++i)
            {
                double val = diagonal[i];
                if (Math.Abs(val) <= tolerance) throw new SingularMatrixException($"Zero diagonal entry at index {i}");
                inverseDiagonal[i] = 1.0 / val;
            }
        }

        /// <summary>
        /// The number of rows/columns of the preconditioner and the original matrix
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// See <see cref="IPreconditioner.SolveLinearSystem(Vector)"/>
        /// </summary>
        /// <param name="rhs"></param>
        public Vector SolveLinearSystem(Vector rhs)
        {
            Preconditions.CheckSystemSolutionDimensions(Order, Order, rhs.Length);
            double[] solution = new double[Order];
            for (int i = 0; i < Order; ++i) solution[i] = inverseDiagonal[i] * rhs[i];
            return Vector.CreateFromArray(solution);
        }
    }
}
