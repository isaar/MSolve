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
    public class CholeskyFactorization: IFactorization
    {
        private readonly double[] data;

        internal CholeskyFactorization(double[] upperData, int order)
        {
            this.data = upperData;
            this.Order = order;
        }

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
            return det;
        }

        public TriangularUpper GetUpperTriangle()
        {
            return TriangularUpper.CreateFromArray(data, true);
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
