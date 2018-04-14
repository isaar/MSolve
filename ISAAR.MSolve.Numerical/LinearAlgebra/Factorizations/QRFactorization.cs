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

// TODO: Perhaps I should use the QR with pivoting
// TODO: Add the option to specify if the diagonal entries of A are non diagonal.
// TODO: Handle the case when m < n. That would be a RQ factorization.
// TODO: I think least squares only work if the columns of A are independent. This is not taken care of so far.
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    /// <summary>
    /// QR factorization of a full matrix A. A = Q*R. For more see:
    /// https://software.intel.com/en-us/mkl-developer-reference-c-orthogonal-factorizations-lapack-computational-routines#E832D468-0891-40EC-9468- ,
    /// http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga3766ea903391b5cf9008132f7440ec7b.html#ga3766ea903391b5cf9008132f7440ec7b.
    /// </summary>
    public class QRFactorization
    {
        private readonly double[] reflectorsAndR;
        private readonly double[] tau;

        private QRFactorization(int numRows, int numCols, double[] reflectorsAndR, double[] tau)
        {
            this.NumRows = numRows;
            this.NumColumns = numCols;
            this.reflectorsAndR = reflectorsAndR;
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
        /// Calculates the QR factorization of a matrix, such that A = Q*R. Requires an extra 
        /// min(<paramref name="numRows"/>, <paramref name="numCols"/>) available memory.
        /// </summary>
        /// <param name="numRows">The number of rows of the matrix.</param>
        /// <param name="numCols">The number of columns of the matrix.</param>
        /// <param name="matrix">The internal buffer storing the matrix entries in column major order. It will 
        ///     be overwritten with factorization data.</param>
        /// <returns></returns>
        public static QRFactorization Factorize(int numRows, int numCols, double[] matrix)
        {
            if (numRows < numCols)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }

            // Call MKL
            int m = numRows;
            int n = numCols;
            int ldA = m;
            int minDim = Math.Min(m, n); //TODO: this is known to be numCols (though it may change in the future)
            double[] tau = new double[minDim];
            int info = MKLUtilities.DefaultInfo;
            info = LAPACKE.Dgeqrf(LAPACKE.LAPACK_COL_MAJOR, m, n, matrix, ldA, tau);

            // Check MKL execution
            if (info == 0) return new QRFactorization(numRows, numCols, matrix, tau);
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }

        public Matrix GetFactorQ()
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }
            int m = NumRows;
            int n = m;
            int k = NumColumns;
            int ldA = m;
            // We need a larger buffer for Q (m-by-m) > reflectorsAndR (m-by-p)
            double[] q = ArrayManipulations.ResizeColsColMajor(NumRows, NumColumns, NumRows, reflectorsAndR);
            int info = MKLUtilities.DefaultInfo;
            info = LAPACKE.Dorgqr(LAPACKE.LAPACK_COL_MAJOR, m, n, k, q, ldA, tau);

            // Check MKL execution
            if (info == 0) return Matrix.CreateFromArray(q, NumRows, NumRows, false);
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }

        public Matrix GetFactorR()
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }
            double[] r = Conversions.FullColMajorToFullUpperColMajorRect(NumRows, NumColumns, reflectorsAndR);
            return Matrix.CreateFromArray(r , NumRows, NumColumns, false);
        }

        /// <summary>
        /// Warning: the columns of the original matrix must be independent for this to work.
        /// </summary>
        /// <param name="rhs">A vector that may lie outside the column space of the original matrix. Its length must be equal to 
        ///     <see cref="NumRows"/></param>
        /// <returns></returns>
        public VectorMKL SolveLeastSquares(VectorMKL rhs) //TODO: perhaps I should use the driver routines of LAPACKE
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }
            Preconditions.CheckSystemSolutionDimensions(NumRows, NumColumns, rhs.Length);

            // Least squares: x = inv(A^T * A) * A^T * b = inv(R) * Q^T * b, where b is the right hand side vector. 
            // Step 1: c = Q^T * b. Q is m-by-m, b is m-by-1 => c is m-by-1
            double[] c = rhs.CopyToArray();
            int m = rhs.Length;
            int nRhs = 1; // rhs = m-by-1
            int k = tau.Length;
            int ldA = m;
            int ldC = m;
            int infoMult = MKLUtilities.DefaultInfo;
            infoMult = LAPACKE.Dormqr(LAPACKE.LAPACK_COL_MAJOR, LAPACKE.LAPACK_SIDE_LEFT, LAPACKE.LAPACK_TRANSPOSE,
                m, nRhs, k, reflectorsAndR, ldA, tau, c, ldC);

            // Check MKL execution
            if (infoMult != 0) throw MKLUtilities.ProcessNegativeInfo(infoMult); // info < 0. This function does not return info > 0
            
            // Step 2: R * x = c, with R being m-by-n and upper trapezoidal (because m >= n).
            // Decomposing R: [R1; 0] * x = [c1 ; c2 ] => R1 * x = c1 => R1 * x = c1, with R1 being n-by-n, upper triangular
            // and stored in the factorized col major data. c1 and x are both n-by-1. The information stored in c2 is lost due to
            // the least squares approximation.
            // TODO: I do not really need to discard the extra m-n terms of c2, but I think it is unsafe to carry them around and
            // risk some method of Vector using the length of the internal buffer, instead of its Length propert.
            if (NumRows > NumColumns) Array.Resize<double>(ref c, NumColumns); 
            int n = NumColumns; // Order of matrix R1
            int ldR = NumRows; // R1 is stored in the upper trapezoid of a NumRows * NumColumns col major array.
            int incC = 1; // step in rhs array, which is the same c used for Q^T * b
            CBlas.Dtrsv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit,
                n, ref reflectorsAndR[0], ldR, ref c[0], incC);
            // TODO: Check output of BLAS somehow. E.g. Zero diagonal entries will result in NaN in the result vector.
            return VectorMKL.CreateFromArray(c, false);
        }
    }
}
