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
using ISAAR.MSolve.Numerical.MKL;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    /// <summary>
    /// Cholesky factorized of a matrix stored in full column major format. Only the upper part of the matrix is factorized. The 
    /// lower part is still stored in the array but it is ignored.
    /// </summary>
    public class CholeskyFull
    {
        private readonly double[] data;

        private CholeskyFull(int order, double[] data)
        {
            this.Order = order;
            this.data = data;
        }

        /// <summary>
        /// The number of rows = number of columns of the original matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Calculates the cholesky factorization of a square matrix, such that A = transpose(U) * U. If the matrix is not 
        /// positive definite an <see cref="IndefiniteMatrixException"/> will be thrown.
        /// </summary>
        /// <param name="order"></param>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static CholeskyFull Factorize(int order, double[] matrix)
        {
            // Call MKL
            int info = MKLUtilities.DefaultInfo;
            Lapack.Dpotrf("U", ref order, ref matrix[0], ref order, ref info);

            // Check MKL execution
            if (info == MKLUtilities.DefaultInfo)
            {
                // first check the default info value, since it lies in the other intervals.
                // info == default => the MKL call did not succeed. 
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
            else if (info > 0)
            {
                string msg = $"The leading minor of order {info-1} (and therefore the matrix itself) is not"
                + " positive-definite, and the factorization could not be completed.";
                throw new IndefiniteMatrixException(msg);
            }

            return new CholeskyFull(order, matrix);
        }

        /// <summary>
        /// Calculates the determinant of the original matrix. A = U^T*U => det(A) = det(U^T)* det(U) => det(A) = (det(U))^2,
        /// where det(U) = U[0,0] * U[1,1] * ... * U[n,n]
        /// </summary>
        /// <returns></returns>
        public double CalcDeterminant()
        {
            double det = 1.0;
            for (int i = 0; i < Order; ++i)
            {
                det *= data[i * Order + i];
            }
            return det * det;
        }

        public Matrix GetFactorU()
        {
            double[] u = Conversions.FullColMajorToFullUpperColMajor(data, false);
            return Matrix.CreateFromArray(u, Order, Order, false);
        }

        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            Preconditions.CheckSystemSolutionDimensions(this.Order, this.Order, rhs.Length);

            // Back & forward substitution using MKL
            int n = Order;
            double[] b = rhs.CopyToArray();
            int info = MKLUtilities.DefaultInfo;
            int nRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int ldb = n; // column major ordering: leading dimension of b is n 
            Lapack.Dpotrs("U", ref n, ref nRhs, ref data[0], ref n, ref b[0], ref ldb, ref info);

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
