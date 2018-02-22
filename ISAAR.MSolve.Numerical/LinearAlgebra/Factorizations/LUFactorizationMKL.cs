using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    /// <summary>
    /// The LU factorization of a matrix A consists of a lower triangular matrix L (with 1 in its diagonal entries), an upper 
    /// triangular matrix U and a permutation matrix P, such that A = P*L*U. This class stores L,U,P in an efficient manner and
    /// provides common methods to use them. 
    /// </summary>
    public class LUFactorizationMKL: IFactorizationMKL
    {
        private readonly double[] lu;
        private readonly int[] p;
        private readonly int firstZeroPivot;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="order"></param>
        /// <param name="lowerUpper"></param>
        /// <param name="permutation"></param>
        /// <param name="firstZeroPivot">Will be negative if the matrix is invertible.</param>
        internal LUFactorizationMKL(int order, double[] lowerUpper, int[] permutation, int firstZeroPivot)
        {
            this.Order = order;
            this.lu = lowerUpper;
            this.p = permutation;
            this.firstZeroPivot = firstZeroPivot;
            if (firstZeroPivot < 0) this.IsSingular = false;
            else this.IsSingular = true;
        }

        public bool IsSingular { get; }

        /// <summary>
        /// The number of rows or columns of the matrix. 
        /// </summary>
        public int Order { get; }

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
                    det *= lu[i * Order + i];
                }
                return det;
            }
        }

        // Explicitly composes and returns the L and U matrices. TODO: 1) Use matrix classes instead of arrays, 2) Also return P
        public (SquareMatrixMKL lower, SquareMatrixMKL upper) Expand()
        {
            double[] l = Conversions.FullColMajorToFullLowerColMajor(lu, true);
            double[] u = Conversions.FullColMajorToFullUpperColMajor(lu, false);
            return (SquareMatrixMKL.CreateFromArray(l, Order, false), SquareMatrixMKL.CreateFromArray(u, Order, false));
        }

        public DenseVector SolveLinearSystem(DenseVector rhs)
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
            Lapack.Dgetrs("N", ref n, ref nRhs, ref lu[0], ref n, ref p[0], ref b[0], ref ldb, ref info);

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

            return DenseVector.CreateFromArray(b, false);
        }
    }
}
