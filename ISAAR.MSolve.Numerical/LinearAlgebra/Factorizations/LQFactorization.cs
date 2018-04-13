using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.MKL;

// TODO: Handle the case when n < m. That would be a QL factorization.
// TODO: I think least squares only work if the columns of A are independent. This is not taken care of so far.
// TODO: Supposedly MKL's LQ is much slower than QR of transpose(A) (see 
// https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/505731). Since they both solve the same problem, it 
// might be worth looking into it.
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    /// <summary>
    /// LQ factorization of a full matrix A. A = L*Q. Not that this is equivalent to taking the QR factorization of A^T:
    /// A^T = Qo * Ro => A = Ro^T * Qo^T => L = Ro^T and Q = Qo^T. For more see:
    /// https://software.intel.com/en-us/mkl-developer-reference-c-orthogonal-factorizations-lapack-computational-routines#E832D468-0891-40EC-9468- ,
    /// </summary>
    public class LQFactorization
    {
        private readonly double[] reflectorsAndL;
        private readonly double[] tau;

        private LQFactorization(int numRows, int numCols, double[] reflectorsAndL, double[] tau)
        {
            this.NumRows = numRows;
            this.NumColumns = numCols;
            this.reflectorsAndL = reflectorsAndL;
            this.tau = tau;
        }

        /// <summary>
        /// The number of columns of the original matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the original matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// Calculates the LQ factorization of a matrix, such that A = L*Q. Requires an extra 
        /// min(<paramref name="numRows"/>, <paramref name="numCols"/>) available memory.
        /// </summary>
        /// <param name="numRows">The number of rows of the matrix.</param>
        /// <param name="numCols">The number of columns of the matrix.</param>
        /// <param name="matrix">The internal buffer storing the matrix entries in column major order. It will 
        ///     be overwritten with factorization data.</param>
        /// <returns></returns>
        public static LQFactorization Factorize(int numRows, int numCols, double[] matrix)
        {
            if (numRows > numCols)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }

            // Call MKL
            int m = numRows;
            int n = numCols;
            int ldA = m;
            int minDim = Math.Min(m, n); //TODO: this is known to be numRows (though it may change in the future)
            double[] tau = new double[minDim];
            int info = MKLUtilities.DefaultInfo;
            info = LAPACKE.Dgelqf(LAPACKE.LAPACK_COL_MAJOR, m, n, matrix, ldA, tau);

            // Check MKL execution
            if (info == 0) return new LQFactorization(numRows, numCols, matrix, tau);
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }

        public Matrix GetFactorL()
        {
            if (NumRows > NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }
            double[] r = Conversions.FullColMajorToFullLowerColMajorRect(NumRows, NumColumns, reflectorsAndL);
            return Matrix.CreateFromArray(r, NumRows, NumColumns, false);
        }

        public Matrix GetFactorQ()
        {
            if (NumRows > NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }
            int n = NumColumns;
            int m = n;
            int k = NumRows;
            int ldA = m;

            // We need a larger buffer for Q (n-by-n) > reflectorsAndL (p-by-n)
            double[] q = ArrayManipulations.IncreaseRowsColMajor(NumRows, NumColumns, NumColumns, reflectorsAndL);
            int info = MKLUtilities.DefaultInfo;
            info = LAPACKE.Dorglq(LAPACKE.LAPACK_COL_MAJOR, m, n, k, q, ldA, tau);

            // Check MKL execution
            if (info == 0) return Matrix.CreateFromArray(q, NumColumns, NumColumns, false);
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }
        
    }
}
