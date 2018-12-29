using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;

// TODO: check if the last minor is non-negative, during factorization. Is it possible that it isn't. Does it affect system 
// solution or inversion?
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    /// <summary>
    /// Cholesky factorization of a symmetric positive definite matrix, stored in packed column major format. Only the upper
    /// triangle part of the matrix is stored and factorized. Uses LAPACK.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CholeskyPacked: ITriangulation
    {
        private readonly double[] data;

        private CholeskyPacked(int order, double[] upperData)
        {
            this.data = upperData;
            this.Order = order;
        }

        /// <summary>
        /// If true, the internal data of this object are overwritten and used by another object. No property or method of
        /// this object must be called as it would throw exceptions or lead to data corruption. If false, this object can be 
        /// used normally.
        /// </summary>
        public bool IsOverwritten { get; private set; }

        /// <summary>
        /// The number of rows/columns of the original square matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Calculates the cholesky factorization of a symmetric positive definite matrix, such that A = transpose(U) * U. 
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix.</param>
        /// <param name="matrix">The entries of the original symmetric matrix in packed column major layout.</param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not symmetric positive definite.</exception>
        public static CholeskyPacked Factorize(int order, double[] matrix)
        {
            // Call LAPACK
            int indefiniteMinorIdx = LapackLinearEquations.Dpptrf(StoredTriangle.Upper, order, matrix, 0);

            // Check LAPACK execution
            if (indefiniteMinorIdx < 0) return new CholeskyPacked(order, matrix);
            else
            {
                string msg = $"The leading minor of order {indefiniteMinorIdx} (and therefore the matrix itself) is not"
                    + " positive-definite, and the factorization could not be completed.";
                throw new IndefiniteMatrixException(msg);
            }
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        /// <remarks>
        /// A = U^T*U => det(A) = det(U^T)* det(U) => det(A) = (det(U))^2, where det(U) = U[0,0] * U[1,1] * ... * U[n,n]
        /// </remarks>
        public double CalcDeterminant()
        {
            CheckOverwritten();
            double det = 1.0;
            for (int i = 0; i < Order; ++i)
            {
                det *= data[i + (i * (i + 1)) / 2];
            }
            return det * det;
        }

        /// <summary>
        /// Explicitly creates the upper triangular matrix U that resulted from the Cholesky factorization: A = transpose(U) * U,
        /// where A and U are n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public TriangularUpper GetFactorU()
        {
            CheckOverwritten();
            return TriangularUpper.CreateFromArray(Order, data, true);
        }

        /// <summary>
        /// Calculates the inverse of the original matrix and returns it in a new <see cref="SymmetricMatrix"/> instance. 
        /// WARNING: If <paramref name="inPlace"/> is set to true, this object must not be used again, otherwise a 
        /// <see cref="InvalidOperationException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">False, to copy the internal factorization data before inversion. True, to overwrite it with
        ///     the inverse matrix, thus saving memory and time. However, that will make this object unusable, so you MUST NOT 
        ///     call any other members afterwards.</param>
        public SymmetricMatrix Invert(bool inPlace)
        {
            CheckOverwritten();

            // Call LAPACK
            double[] inverse; // if A is posdef, so is inv(A)
            if (inPlace)
            {
                inverse = data;
                IsOverwritten = true;
            }
            else
            {
                inverse = new double[data.Length];
                Array.Copy(data, inverse, data.Length);
            }
            int indefiniteMinorIdx = LapackLinearEquations.Dpptri(StoredTriangle.Upper, Order, inverse, 0);

            // Check LAPACK execution
            if (indefiniteMinorIdx < 0)
            {
                return SymmetricMatrix.CreateFromArray(inverse, Order, DefiniteProperty.PositiveDefinite);
            }
            else  // this should not have happened
            {
                throw new IndefiniteMatrixException($"The entry ({indefiniteMinorIdx}, {indefiniteMinorIdx}) of the factor U"
                    + " is 0 and the inverse could not be computed.");
            }
        }

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
        /// </summary>
        /// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid arguments.</exception>
        public void SolveLinearSystem(Vector rhs, Vector solution)
        {
            CheckOverwritten();
            Preconditions.CheckSystemSolutionDimensions(Order, rhs.Length);
            Preconditions.CheckMultiplicationDimensions(Order, solution.Length);

            // Call LAPACK
            solution.CopyFrom(rhs);
            int numRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int leadingDimB = Order; // column major ordering: leading dimension of b is n 
            LapackLinearEquations.Dpptrs(StoredTriangle.Upper, Order, numRhs, data, 0, solution.RawData, 0, leadingDimB);
        }

        private void CheckOverwritten()
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");
        }
    }
}
