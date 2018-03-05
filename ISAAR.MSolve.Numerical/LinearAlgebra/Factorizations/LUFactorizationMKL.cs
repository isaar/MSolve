using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    /// <summary>
    /// The LU factorization of a matrix A consists of a lower triangular matrix L (with 1 in its diagonal entries), an upper 
    /// triangular matrix U and a permutation matrix P, such that A = P*L*U. This class stores L,U,P in an efficient manner and
    /// provides common methods to use them. 
    /// </summary>
    public class LUFactorizationMKL: IFactorizationMKL
    {
        // Perhaps a smaller tolerance is appropriate, since the "almost zero" will propagate during back & forward substitution.
        private const double PivotTolerance = 1e-13;

        private readonly double[] lowerUpper;
        private readonly int[] permutation;
        private readonly int firstZeroPivot;

        private LUFactorizationMKL(int order, double[] lowerUpper, int[] permutation, int firstZeroPivot, bool isSingular)
        {
            this.Order = order;
            this.lowerUpper = lowerUpper;
            this.permutation = permutation;
            this.firstZeroPivot = firstZeroPivot;
            this.IsSingular = isSingular;
        }

        public bool IsSingular { get; }

        /// <summary>
        /// The number of rows or columns of the matrix. 
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Calculates the LUP factorization of a square matrix, such that A = P * L * U. Requires an extra O(n^2 + n) 
        /// available memory.
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix.</param>
        /// <param name="originalColMajorMatrix">The internal buffer stroring the matrix entries in column major order. It will 
        ///     be copied and not altered, thus you should not copy it yourself.</param>
        /// <param name="pivotTolerance"></param>
        /// <returns></returns>
        public static LUFactorizationMKL CalcFactorization(int order, double[] originalColMajorMatrix, 
            double pivotTolerance = LUFactorizationMKL.PivotTolerance)
        {
            // Copy matrix. This may exceed available memory and needs an extra O(n^2) accesses. 
            // To avoid these, use the ~InPlace version.
            double[] lowerUpper = new double[originalColMajorMatrix.Length];
            Array.Copy(originalColMajorMatrix, lowerUpper, originalColMajorMatrix.Length);

            // Call MKL
            int[] permutation = new int[order];
            int info = MKLUtilities.DefaultInfo;
            Lapack.Dgetrf(ref order, ref order, ref lowerUpper[0], ref order, ref permutation[0], ref info);

            // Check MKL execution
            if (info == MKLUtilities.DefaultInfo)
            {
                // first check the default info value, since it lies in the other intervals.
                // info == dafeult => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                throw new MKLException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (info < 0)
            {
                string msg = string.Format("The {0}th parameter has an illegal value.", -info)
                    + " Please contact the developer responsible for the linear algebra project.";
                throw new MKLException(msg);
            }

            int firstZeroPivot = int.MinValue;
            if (info > 0) firstZeroPivot = info - 1;
            else if (Math.Abs(lowerUpper[order * order - 1]) <= pivotTolerance)
            {
                // False Negative: info = 0, but LAPACK doesn't check the last diagonal entry!
                firstZeroPivot = order - 1;
            }
            return new LUFactorizationMKL(order, lowerUpper, permutation, firstZeroPivot, (firstZeroPivot >= 0));
        }

        /// <summary>
        /// Calculates the determinant of the original matrix. det(A) = det(L*U) = det(L)*det(U). 
        /// Since all these are triangular matrices their determinants is the product of their diagonal entries:
        /// det(L) = 1*1*...*1 = 1. Thus det(A) = det(U) = U1*U2*...*Un
        /// </summary>
        /// <returns>The determinant of the original matrix.</returns>
        public double CalcDeterminant()
        {
            if (IsSingular) return 0.0;
            else
            {
                double det = 1.0;
                for (int i = 0; i < Order; ++i)
                {
                    det *= lowerUpper[i * Order + i];
                }
                return det;
            }
        }

        // Explicitly composes and returns the L and U matrices. TODO: 1) Use matrix classes instead of arrays, 2) Also return P
        public (MatrixMKL lower, MatrixMKL upper) Expand()
        {
            double[] l = Conversions.FullColMajorToFullLowerColMajor(lowerUpper, true);
            double[] u = Conversions.FullColMajorToFullUpperColMajor(lowerUpper, false);
            return (MatrixMKL.CreateFromArray(l, Order, Order, false), MatrixMKL.CreateFromArray(u, Order, Order, false));
        }

        /// <summary>
        /// Inverts the original matrix. The inverse is stored in a different buffer than the factorized matrix.
        /// </summary>
        /// <returns></returns>
        public MatrixMKL Invert()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Inverts the original matrix without creating a new buffer. The inverse is replaces the factorized matrix. Therefore 
        /// this object should not be used again.
        /// </summary>
        /// <returns></returns>
        public MatrixMKL InvertInPlace()
        {
            // would be cool if I could set the data structures to 0
            throw new NotImplementedException();
        }

        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            // Check if the matrix is singular first
            if (IsSingular)
            {
                string msg = "The factorization has been completed, but U is singular."
                    + $" The first zero pivot is U[{firstZeroPivot}, {firstZeroPivot}] = 0.";
                throw new SingularMatrixException(msg);
            }

            // Back & forward substitution using MKL
            int n = Order;
            double[] b = rhs.CopyToArray();
            int info = MKLUtilities.DefaultInfo;
            int nRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int ldb = n; // column major ordering: leading dimension of b is n 
            Lapack.Dgetrs("N", ref n, ref nRhs, ref lowerUpper[0], ref n, ref permutation[0], ref b[0], ref ldb, ref info);

            // Check MKL execution
            if ((info == MKLUtilities.DefaultInfo) || (info > 0)) 
            {
                // first check the default info value, since it lies in the other intervals.
                // info == dafeult => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                throw new MKLException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (info < 0)
            {
                string msg = $"The {-info}th parameter has an illegal value." 
                    + " Please contact the developer responsible for the linear algebra project.";
                throw new MKLException(msg);
            }

            return VectorMKL.CreateFromArray(b, false);
        } 
    }
}
