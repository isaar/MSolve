using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;

//TODO: When returning L & U, use Triangular matrices. Also return P. Also L, U should be TriangularLower and TriangularUpper
//TODO: Is the determinant affected by the permutation P? I think P changes the sign, depending on how many row exhanges there 
//      are. I did not take it into account in the implementation or the documentation.
//TODO: Computing the determinant may overflow the double variable. It should be done in a C .dll and the overflow should be 
//      handled there by having a naive (fast) and a safe version.
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// The LU factorization of a matrix A with partial pivoting (row exchanges) consists of a lower triangular matrix L 
    /// (with 1 in its diagonal entries), an upper triangular matrix U and a permutation matrix P, such that A = P*L*U. This 
    /// class stores L,U,P in an efficient manner and provides common methods to use them. A must be square. Uses LAPACK.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class LUFactorization: ITriangulation
    {
        /// <summary>
        /// The default value under which a diagonal entry (pivot) is considered to be 0 during Cholesky factorization.
        /// </summary>
        private const double PivotTolerance = 1e-13;  //TODO: Perhaps a smaller tolerance is appropriate, since the "almost zero" will propagate during back & forward substitution.

        private readonly double[] lowerUpper;

        /// <summary>
        /// Its entries are in 1-based indexing. This is not the actual permutation vector. Instead it records the row exchanges. 
        /// E.g. { 4, 4, 4, 4} means that rows 0, 1, 2 have been switched with the last row (the last row doesn't change).
        /// </summary>
        private readonly int[] rowExchanges;
        private readonly int firstZeroPivot;
        private readonly double pivotTolerance;

        private LUFactorization(int order, double[] lowerUpper, int[] permutation, int firstZeroPivot, bool isSingular, 
            double pivotTolerance)
        {
            this.Order = order;
            this.lowerUpper = lowerUpper;
            this.rowExchanges = permutation;
            this.firstZeroPivot = firstZeroPivot;
            this.pivotTolerance = pivotTolerance;
            this.IsSingular = isSingular;
            this.IsOverwritten = false;
        }

        /// <summary>
        /// If true, the internal data of this object are overwritten and used by another object. No property or method of
        /// this object must be called as it would throw exceptions or lead to data corruption. If false, this object can be 
        /// used normally.
        /// </summary>
        public bool IsOverwritten { get; private set; }

        /// <summary>
        /// If true, the original matrix before the factorization is not invertible. In this case, 
        /// <see cref="SolveLinearSystem(Vector)"/> and <see cref="Invert(bool)"/> will throw an exception. If false, the 
        /// original matrix is invertible and those methods are safe to call.
        /// </summary>
        public bool IsSingular { get; }

        /// <summary>
        /// The number of rows/columns of the original square matrix. 
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Calculates the LUP factorization of a square matrix, such that A = P * L * U. Requires an extra O(n) available 
        /// memory, where n is the <paramref name="order"/>.
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix.</param>
        /// <param name="matrix">The internal buffer stroring the matrix entries in column major order. It will 
        ///     be overwritten.</param>
        /// <param name="pivotTolerance">If a diagonal entry (called pivot) is &lt;= <paramref name="pivotTolerance"/> it will be  
        ///     considered as zero and a permutation will be used to find a non-zero pivot (the process is called pivoting).
        ///     </param>
        public static LUFactorization Factorize(int order, double[] matrix,
            double pivotTolerance = LUFactorization.PivotTolerance)
        {
            int[] rowExchanges = new int[order];
            int firstZeroPivot = LapackLinearEquations.Dgetrf(order, order, matrix, 0, order, rowExchanges, 0, pivotTolerance);
            return new LUFactorization(order, matrix, rowExchanges, firstZeroPivot, firstZeroPivot > 0, pivotTolerance);
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        /// <remarks>
        /// det(A) = det(L*U) = det(L)*det(U). Since all these are triangular matrices their determinants is the product of their 
        /// diagonal entries: det(L) = 1*1*...*1 = 1. Thus det(A) = det(U) = U1*U2*...*Un
        /// The sign is positive for an even number of row exchanges and negative for an odd number. The row exchanges are 
        /// recorded in <see cref="rowExchanges"/>. For more details see 
        /// https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/309460.
        /// </remarks>
        public double CalcDeterminant() //TODO: this should be implemented in a dll
        {
            CheckOverwritten();
            if (IsSingular) return 0.0;
            else
            {
                double det = 1.0;
                for (int i = 0; i < Order; ++i)
                {
                    if (rowExchanges[i] == i + 1) det *= lowerUpper[i * Order + i];
                    else det *= - lowerUpper[i * Order + i];
                }
                return det;
            }
        }

        /// <summary>
        /// Explicitly creates the lower triangular matrix L that resulted from the LU factorization: A = P * L * U,
        /// where A, L, U and P are n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetFactorL()
        {
            CheckOverwritten();
            double[] l = Conversions.FullColMajorToFullLowerColMajor(lowerUpper, true);
            return Matrix.CreateFromArray(l, Order, Order, false);
        }

        /// <summary>
        /// Explicitly creates the upper triangular matrix U that resulted from the LU factorization: A = P * L * U,
        /// where A, L, U and P are n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetFactorU()
        {
            CheckOverwritten();
            double[] u = Conversions.FullColMajorToFullUpperColMajor(lowerUpper, false);
            return Matrix.CreateFromArray(u, Order, Order, false);
        }

        /// <summary>
        /// Calculates the inverse of the original square matrix and returns it in a new <see cref="Matrix"/> instance. This
        /// only works if the original matrix is not singular, which can be checked through <see cref="IsSingular"/>.
        /// WARNING: If <paramref name="inPlace"/> is set to true, this object must not be used again, otherwise a 
        /// <see cref="InvalidOperationException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">False, to copy the internal factorization data before inversion. True, to overwrite it with
        ///     the inverse matrix, thus saving memory and time. However, that will make this object unusable, so you MUST NOT 
        ///     call any other members afterwards.</param>
        /// <exception cref="SingularMatrixException">Thrown if the original matrix is not invertible.</exception>
        public Matrix Invert(bool inPlace)
        {
            // Check if the matrix is suitable for inversion
            CheckOverwritten();
            if (IsSingular) throw new SingularMatrixException("The factorization has been completed, but U is singular."
                + $" The first zero pivot is U[{firstZeroPivot}, {firstZeroPivot}] = 0.");

            // Copy if the matrix if the inversion will be in place.
            double[] inverse;
            if (inPlace)
            {
                inverse = lowerUpper;
                IsOverwritten = true;
            }
            else
            {
                inverse = new double[lowerUpper.Length];
                Array.Copy(lowerUpper, inverse, lowerUpper.Length);
            }

            // Call LAPACK
            int firstZeroDiagonal = LapackLinearEquations.Dgetri(Order, inverse, 0, Order, rowExchanges, 0, pivotTolerance);
            if (firstZeroDiagonal > 0) // This should not have happened though
            {
                throw new SingularMatrixException($"The ({firstZeroDiagonal}, {firstZeroDiagonal}) element of factor U is zero,"
                    + " U is singular and the inversion could not be completed.");
            }

            return Matrix.CreateFromArray(inverse, Order, Order, false);
        }

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
        /// </summary>
        /// <remarks>
        /// This method is not garanteed to succeed. A singular matrix can be factorized as A=P*L*U, but not all linear systems
        /// with a singular matrix can be solved.
        /// </remarks>
        /// <exception cref="SingularMatrixException">Thrown if the original matrix is not invertible.</exception>
        /// <exception cref="LapackException">Thrown if the call to LAPACK fails due to invalid arguments.</exception>
        public void SolveLinearSystem(Vector rhs, Vector solution)
        {
            CheckOverwritten();
            Preconditions.CheckSystemSolutionDimensions(Order, rhs.Length);
            Preconditions.CheckMultiplicationDimensions(Order, solution.Length);

            // Check if the matrix is singular first
            if (IsSingular)
            {
                string msg = "The factorization has been completed, but U is singular."
                    + $" The first zero pivot is U[{firstZeroPivot}, {firstZeroPivot}] = 0.";
                throw new SingularMatrixException(msg);
            }

            // Back & forward substitution using LAPACK
            int n = Order;
            solution.CopyFrom(rhs); //double[] solution = rhs.CopyToArray();
            int numRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int leadingDimB = n; // column major ordering: leading dimension of b is n 
            LapackLinearEquations.Dgetrs(TransposeMatrix.NoTranspose, n, numRhs, lowerUpper, 0, n, rowExchanges, 0, 
                solution.RawData, 0, leadingDimB);
        }

        private void CheckOverwritten()
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");
        }
    }
}
