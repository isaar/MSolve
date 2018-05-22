using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.ArrayManipulations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.MKL;
using ISAAR.MSolve.LinearAlgebra.Commons;

// TODO: Handle the case when n < m. That would be a QL factorization.
// TODO: I think least squares only work if the columns of A are independent. This is not taken care of so far.
// TODO: Supposedly MKL's LQ is much slower than QR of transpose(A) (see 
// https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/505731). Since they both solve the same problem, it 
// might be worth looking into it.
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
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
            double[] r = Conversions.RectColMajorToRectLowerColMajor(NumRows, NumColumns, reflectorsAndL);
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
            double[] q = ArrayColMajor.IncreaseRows(NumRows, NumColumns, NumColumns, reflectorsAndL);
            int info = MKLUtilities.DefaultInfo;
            info = LAPACKE.Dorglq(LAPACKE.LAPACK_COL_MAJOR, m, n, k, q, ldA, tau);

            // Check MKL execution
            if (info == 0) return Matrix.CreateFromArray(q, NumColumns, NumColumns, false);
            else throw MKLUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }

        /// <summary>
        /// From the infinite solutions of A*x=b, find the one with the minimum norm2(x). Warning: the rows of the original 
        /// matrix must be independent for this to work: If A (m-by-n, with m $lt;= n) has full row rank, r = m => the column 
        /// space of A is an m dimensional subspace of R^m (every vector A*z is m-by-1) => the column space of A is the whole
        /// R^m and A*x=b, always has at least one solution. The nullspace of A is an (n-m) dimensional subspace of R^n. Thus 
        /// there are infinite solutions to A*x=b.
        /// </summary>
        /// <param name="rhs">The right hand side vector. Its length must be that may lie outside the column space of the original matrix. Its length must be equal to 
        ///     <see cref="NumRows"/></param>
        /// <returns></returns>
        public Vector SolveMinNorm(Vector rhs) //TODO: perhaps I should use the driver routines of LAPACKE
        {
            if (NumRows > NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }
            Preconditions.CheckSystemSolutionDimensions(NumRows, NumColumns, rhs.Length);

            // Min norm: A * x = b => L * Q * x = b, where b is the right hand side vector. 
            // Step 1: L * c = b. L is m-by-n, b is m-by-1 => c is n-by-1. Reminder n>=m.
            // Decomposing L: [L1 0] * [c1; c2] = b => L1*c1 = b, where: 
            // L1 is m-by-m, lower triangular and stored in the factorized col major dat
            // c1 is m-by-1 and can be found by forward substitution
            // c2 is (n-m)-by-1 and can be any vector whatsoever.
            int ldA = NumRows; // Also k = min(NumRows, NumColumns) = NumRows = ldA
            int incC = 1; // step in rhs array
            double[] c = new double[NumColumns];
            rhs.CopyToArray(0, c, 0, NumRows);
            CBlas.Dtrsv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasLower, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit,
                NumRows, ref reflectorsAndL[0], ldA, ref c[0], incC);
            // TODO: Check output of BLAS somehow. E.g. Zero diagonal entries will result in NaN in the result vector.

            // Step 2: x = Q^T * c. Q is n-by-n, c is n-by-1, x is n-by-1
            int m = NumColumns; // c has been grown past the limits of rhs.
            int nRhs = 1; // rhs = m-by-1
            int k = tau.Length;
            int ldC = m;
            int infoMult = MKLUtilities.DefaultInfo;
            infoMult = LAPACKE.Dormlq(LAPACKE.LAPACK_COL_MAJOR, LAPACKE.LAPACK_SIDE_LEFT, LAPACKE.LAPACK_TRANSPOSE,
                m, nRhs, k, reflectorsAndL, ldA, tau, c, ldC);

            // Check MKL execution
            if (infoMult == 0) return Vector.CreateFromArray(c, false);
            else throw MKLUtilities.ProcessNegativeInfo(infoMult); // info < 0. This function does not return info > 0
        }
    }
}
