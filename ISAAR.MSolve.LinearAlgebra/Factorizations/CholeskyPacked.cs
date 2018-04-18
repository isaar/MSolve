using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.MKL;

// TODO: check if the last minor is non-negative, during factorization. Is it possible that it isn't. Does it affect system 
// solution or inversion?
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    public class CholeskyPacked: IFactorization
    {
        private readonly double[] data;

        private CholeskyPacked(int order, double[] upperData)
        {
            this.data = upperData;
            this.Order = order;
        }

        public bool IsOverwritten { get; private set; }

        /// <summary>
        /// The number of rows or columns of the matrix. 
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Calculates the cholesky factorization of a square matrix, such that A = transpose(U) * U. If the matrix is not 
        /// positive definite an <see cref="IndefiniteMatrixException"/> will be thrown.
        /// </summary>
        /// <param name="order"></param>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static CholeskyPacked Factorize(int order, double[] matrix)
        {
            // Call MKL
            int info = MKLUtilities.DefaultInfo;
            Lapack.Dpptrf("U", ref order, ref matrix[0], ref info);

            // Check MKL execution
            if (info == 0) return new CholeskyPacked(order, matrix);
            else if (info > 0)
            {
                string msg = "The leading minor of order " + (info - 1) + " (and therefore the matrix itself) is not"
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
                det *= data[i + (i * (i + 1)) / 2];
            }
            return det * det;
        }

        public TriangularUpper GetFactorU()
        {
            CheckOverwritten();
            return TriangularUpper.CreateFromArray(data, true);
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
        public SymmetricMatrix Invert(bool inPlace)
        {
            CheckOverwritten();

            // Call MKL
            int info = MKLUtilities.DefaultInfo;
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
            info = LAPACKE.Dpptri(LAPACKE.LAPACK_COL_MAJOR, LAPACKE.LAPACK_UPPER, Order, inverse);
            
            // Check MKL execution
            if (info == 0) return SymmetricMatrix.CreateFromArray(inverse, Order, DefiniteProperty.PositiveDefinite);
            else if (info > 0) // this should not have happened
            {
                throw new IndefiniteMatrixException($"The leading minor of order {info - 1} (and therefore the matrix itself)"
                + "is not positive-definite, and the factorization could not be completed.");
            }
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0
        }

        public Vector SolveLinearSystem(Vector rhs)
        {
            CheckOverwritten();

            // Call MKL
            int n = Order;
            double[] b = rhs.CopyToArray();
            int info = MKLUtilities.DefaultInfo;
            int nRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int ldb = n; // column major ordering: leading dimension of b is n 
            Lapack.Dpptrs("U", ref n, ref nRhs, ref data[0], ref b[0], ref ldb, ref info);

            // Check MKL execution
            if (info == 0) return Vector.CreateFromArray(b, false);
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }

        private void CheckOverwritten()
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");
        }
    }
}
