using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Use a dedicated DiagonalMatrix class, instead of passing in double[] or Vector. It will implement the inverse and 
//      multiplication routines. It will also handle distributed matrices. E.g. IDiagonal IMatrixView.GetDiagonal() which will 
//      then have an IDiagonalMatrix.Inverse(). The problem is how we will go from CSR to DiagonalMatrix. Perhaps it would be 
//      better to use the DOK instead.
//TODO: Alternative: instead of demanding the caller to extract the diagonal, this class should read the matrix and only access 
//      its diagonal. I think this alternative is less flexible and more difficult to implement.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning
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
        public const double DefaultTolerance = 1e-10;
        private readonly double[] inverseDiagonal;

        /// <summary>
        /// Initializes a new instance of <see cref="JacobiPreconditioner"/> for the linear system's matrix whose main diagonal
        /// is provided in <paramref name="diagonal"/>.
        /// </summary>
        /// <param name="diagonal">
        /// The main diagonal of the original matrix of the linear system. Constraints: all its entries must be non-zero.
        /// </param>
        /// <param name="tolerance">The value under which a diagonal entry will be considered as zero.</param>
        /// <exception cref="SingularMatrixException">If there is a zero diagonal entry.</exception>
        public JacobiPreconditioner(double[] diagonal, double tolerance = DefaultTolerance)
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
        public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector)
        {
            Preconditions.CheckSystemSolutionDimensions(Order, rhsVector.Length);
            for (int i = 0; i < Order; ++i) lhsVector.Set(i, inverseDiagonal[i] * rhsVector[i]);
        }

        /// <summary>
        /// Creates instances of <see cref="JacobiPreconditioner"/>.
        /// </summary>
        public class Factory: IPreconditionerFactory
        {
            private readonly double tolerance;

            /// <summary>
            /// Initializes a new instance of <see cref="JacobiPreconditioner.Factory"/> with the specified settings.
            /// </summary>
            /// <param name="tolerance">The value under which a diagonal entry will be considered as zero.</param>
            public Factory(double tolerance = JacobiPreconditioner.DefaultTolerance) => this.tolerance = tolerance;

            /// <summary>
            /// See <see cref="IPreconditionerFactory.CreatePreconditionerFor(IMatrixView)"/>.
            /// </summary>
            public IPreconditioner CreatePreconditionerFor(IMatrixView matrix) 
                => new JacobiPreconditioner(matrix.GetDiagonalAsArray(), tolerance);
        }
    }
}
