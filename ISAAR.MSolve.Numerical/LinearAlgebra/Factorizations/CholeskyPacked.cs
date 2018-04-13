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

// TODO: check if the last minor is non-negative, during factorization. Is it possible that it isn't. Does it affect system 
// solution or inversion?
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public class CholeskyPacked: IFactorization
    {
        private readonly double[] data;

        internal CholeskyPacked(double[] upperData, int order)
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
        /// Calculates the determinant of the original matrix. A = U^T*U => det(A) = det(U^T)* det(U) => det(A) = (det(U))^2,
        /// where det(U) = U[0,0] * U[1,1] * ... * U[n,n]
        /// </summary>
        /// <returns></returns>
        public double CalcDeterminant()
        {
            double det = 1.0;
            for (int i = 0; i < Order; ++i)
            {
                det *= data[i + (i * (i + 1)) / 2];
            }
            return det * det;
        }

        public TriangularUpper GetUpperTriangle()
        {
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
            // Check if the matrix is suitable for inversion
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");

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
            else if ((info == MKLUtilities.DefaultInfo) || (info > 0))
            {
                // first check the default info value, since it lies in the other intervals.
                // info == dafeult => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                throw new MKLException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (info < 0)
            {
                throw new MKLException($"The {-info}th parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else // (info > 0) this should not happen
            {
                throw new IndefiniteMatrixException($"The leading minor of order {info - 1} (and therefore the matrix itself)"
                + "is not positive-definite, and the factorization could not be completed.");
            }
        }

        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            // Call MKL
            int n = Order;
            double[] b = rhs.CopyToArray();
            int info = MKLUtilities.DefaultInfo;
            int nRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int ldb = n; // column major ordering: leading dimension of b is n 
            Lapack.Dpptrs("U", ref n, ref nRhs, ref data[0], ref b[0], ref ldb, ref info);

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
