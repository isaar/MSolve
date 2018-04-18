using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.ArrayManipulations;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.MKL;


// TODO: check if the last minor is non-negative, during factorization. Is it possible that it isn't. Does it affect system 
// solution or inversion?
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
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

        public bool IsOverwritten { get; private set; }

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
            if (info == 0) return new CholeskyFull(order, matrix);
            else if (info > 0)
            {
                string msg = $"The leading minor of order {info-1} (and therefore the matrix itself) is not"
                + " positive-definite, and the factorization could not be completed.";
                throw new IndefiniteMatrixException(msg);
            }
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0
        }

        /// <summary>
        /// Calculates the determinant of the original matrix. A = U^T*U => det(A) = det(U^T)* det(U) => det(A) = (det(U))^2,
        /// where det(U) = U[0,0] * U[1,1] * ... * U[n,n]
        /// </summary>
        /// <returns></returns>
        public double CalcDeterminant()
        {
            CheckOverwritten();
            double det = 1.0;
            for (int i = 0; i < Order; ++i)
            {
                det *= data[i * Order + i];
            }
            return det * det;
        }

        public Matrix GetFactorU()
        {
            CheckOverwritten();
            double[] u = Conversions.FullColMajorToFullUpperColMajor(data, false);
            return Matrix.CreateFromArray(u, Order, Order, false);
        }

        /// <summary>
        /// Inverts the original square matrix. The matrix must be positive definite, otherwise an
        /// <see cref="InvalidOperationException"/> will be thrown. If <paramref name="inPlace"/> is set to true, this object must 
        /// not be used again, otherwise a <see cref="InvalidOperationException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">False, to copy the internal factorization data before inversion. True, to overwrite it with
        ///     the inverse matrix, thus saving memory and time. However, that will make this object unusable, so you MUST NOT 
        ///     call any other members afterwards.</param>
        /// <returns></returns>
        public Matrix Invert(bool inPlace)
        {
            CheckOverwritten();

            // Call MKL
            int info = MKLUtilities.DefaultInfo;
            double[] inverse;
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
            info = LAPACKE.Dpotri(LAPACKE.LAPACK_COL_MAJOR, LAPACKE.LAPACK_UPPER, Order, inverse, Order);
            Conversions.CopyUpperToLowerColMajor(inverse, Order); // So far the lower triangle was the same as the original matrix

            // Check MKL execution
            if (info == 0) return Matrix.CreateFromArray(inverse, Order, Order, false);
            else if (info > 0) // this should not have happened
            {
                throw new IndefiniteMatrixException($"The leading minor of order {info - 1} (and therefore the matrix itself)"
                + "is not positive-definite, and the factorization could not be completed.");
            }
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0
        }

        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            CheckOverwritten();
            Preconditions.CheckSystemSolutionDimensions(this.Order, this.Order, rhs.Length);

            // Back & forward substitution using MKL
            int n = Order;
            double[] b = rhs.CopyToArray();
            int info = MKLUtilities.DefaultInfo;
            int nRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int ldb = n; // column major ordering: leading dimension of b is n 
            Lapack.Dpotrs("U", ref n, ref nRhs, ref data[0], ref n, ref b[0], ref ldb, ref info);

            // Check MKL execution
            if (info == 0) return VectorMKL.CreateFromArray(b, false);
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }

        private void CheckOverwritten()
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");
        }
    }
}
