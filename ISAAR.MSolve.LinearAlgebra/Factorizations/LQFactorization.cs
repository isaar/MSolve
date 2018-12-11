using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;

// TODO: Handle the case when n < m. That would be a QL factorization.
// TODO: I think least squares only work if the columns of A are independent. This is not taken care of so far.
// TODO: Supposedly LAPACK's LQ is much slower than QR of transpose(A) (see 
// https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/505731). Since they both solve the same problem, it 
// might be worth looking into it.
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    /// <summary>
    /// LQ factorization of a matrix A: A = L * Q, where A is an m-by-n matrix with m &lt;= n, L is an m-by-n lower trapezoidal 
    /// matrix and Q is an n-by-n orthogonal matrix. Note that this is equivalent to taking the QR factorization of A^T:
    /// A^T = Qo * Ro => A = Ro^T * Qo^T => L = Ro^T and Q = Qo^T. Uses LAPACK. For more see:
    /// https://software.intel.com/en-us/mkl-developer-reference-c-orthogonal-factorizations-lapack-computational-routines#E832D468-0891-40EC-9468-
    /// Authors: Serafeim Bakalakos
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
        /// Calculates the LQ factorization of a matrix, such that A = L * Q. Requires an extra 
        /// min(<paramref name="numRows"/>, <paramref name="numCols"/>) available memory.
        /// </summary>
        /// <param name="numRows">The number of rows of the matrix.</param>
        /// <param name="numCols">The number of columns of the matrix.</param>
        /// <param name="matrix">The internal buffer storing the matrix entries in column major order. It will 
        ///     be overwritten with the factorization data.</param>
        /// <exception cref="NotImplementedException">Thrown if <paramref name="numCols"/> &lt; <paramref name="numRows"/>.
        ///     </exception>
        /// <exception cref="Exceptions.LapackException">Thrown if tha call to LAPACK fails due to an invalid 
        ///     <paramref name="matrix"/>.</exception>
        public static LQFactorization Factorize(int numRows, int numCols, double[] matrix)
        {
            if (numRows > numCols)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }

            // Prepare LAPACK arguments
            int m = numRows;
            int n = numCols;
            int ldA = m;
            int minDim = Math.Min(m, n); //TODO: this is known to be numRows (though it may change in the future)
            double[] tau = new double[minDim];

            // LAPACK query to find the optimum block size. TODO: Use a LAPACKE like wrapper to simplify this.
            int infoQuery = LapackUtilities.DefaultInfo;
            var workQuery = new double[1];
            int lWorkQuery = -1;
            Lapack.Dgelqf(m, n, matrix, 0, ldA, tau, 0, workQuery, 0, lWorkQuery, ref infoQuery);
            CheckLapackExecution(infoQuery);

            // Call LAPACK to actually perform the operation
            int lWorkOptim = (int)(workQuery[0]);
            if (lWorkOptim < 1) lWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lWorkOptim];
            int info = LapackUtilities.DefaultInfo;
            Lapack.Dgelqf(m, n, matrix, 0, ldA, tau, 0, workOptim, 0, lWorkOptim, ref info);
            CheckLapackExecution(info);

            return new LQFactorization(numRows, numCols, matrix, tau);
        }

        /// <summary>
        /// Explicitly creates the lower trapezoidal matrix L that resulted from the factorization: A = L * Q, where A is m-by-n, 
        /// L is m-by-n and Q is n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetFactorL()
        {
            if (NumRows > NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }
            double[] r = Conversions.RectColMajorToRectLowerColMajor(NumRows, NumColumns, reflectorsAndL);
            return Matrix.CreateFromArray(r, NumRows, NumColumns, false);
        }

        /// <summary>
        /// Explicitly creates the orthogonal matrix Q that resulted from the factorization: A = L * Q, where A is m-by-n, 
        /// L is m-by-n and Q is n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetFactorQ()
        {
            if (NumRows > NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }

            // Prepare LAPACK arguments
            int n = NumColumns;
            int m = n;
            int k = NumRows;
            int ldA = m;

            // We need a larger buffer for Q (n-by-n) > reflectorsAndL (p-by-n)
            double[] q = ArrayColMajor.IncreaseRows(NumRows, NumColumns, NumColumns, reflectorsAndL);

            // LAPACK query to find the optimum block size. TODO: Use a LAPACKE like wrapper to simplify this.
            int infoQuery = LapackUtilities.DefaultInfo;
            var workQuery = new double[1];
            int lWorkQuery = -1;
            Lapack.Dorglq(m, n, k, q, 0, ldA, tau, 0, workQuery, 0, lWorkQuery, ref infoQuery);
            CheckLapackExecution(infoQuery);

            // Call LAPACK to actually perform the operation
            int lWorkOptim = (int)(workQuery[0]);
            if (lWorkOptim < 1) lWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lWorkOptim];
            int info = LapackUtilities.DefaultInfo;
            Lapack.Dorglq(m, n, k, q, 0, ldA, tau, 0, workOptim, 0, lWorkOptim, ref info);
            CheckLapackExecution(info);

            return Matrix.CreateFromArray(q, NumColumns, NumColumns, false);
        }

        /// <summary>
        /// Solve the minimum norm problem A * x = b => transpose(Q) * inv(L) * b. The minimum nrom problem is defined as: 
        /// from the infinite solutions of A * x = b, find the one with the minimum norm2(x). Warning: the rows of the original 
        /// matrix must be independent for this to work: If A (m-by-n, with m &lt;= n) has full row rank, r = m => the column 
        /// space of A is an m dimensional subspace of R^m (every vector A*z is m-by-1) => the column space of A is the whole
        /// R^m and A*x=b, always has at least one solution. The nullspace of A is an (n-m) dimensional subspace of R^n. Thus 
        /// there are infinite solutions to A*x=b.
        /// </summary>
        /// <param name="rhsVector">The right hand side vector b. It may lie outside the column space of the original matrix. Its 
        ///     <see cref="IIndexable1D.Length"/> must be equal to this.<see cref="NumRows"/>.</param>
        /// <exception cref="Exceptions.LapackException">Thrown if tha call to LAPACK fails due to <paramref name="rhsVector"/> 
        ///     having a different <see cref="IIndexable1D.Length"/> than this.<see cref="NumRows"/>.</exception>
        public Vector SolveMinNorm(Vector rhsVector) //TODO: perhaps I should use the driver routines of LAPACKE
        {
            if (NumRows > NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be <= the number of columns");
            }
            Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);

            // Min norm: A * x = b => L * Q * x = b, where b is the right hand side vector. 
            // Step 1: L * c = b. L is m-by-n, b is m-by-1 => c is n-by-1. Reminder n>=m.
            // Decomposing L: [L1 0] * [c1; c2] = b => L1*c1 = b, where: 
            // L1 is m-by-m, lower triangular and stored in the factorized col major dat
            // c1 is m-by-1 and can be found by forward substitution
            // c2 is (n-m)-by-1 and can be any vector whatsoever.
            int ldA = NumRows; // Also k = min(NumRows, NumColumns) = NumRows = ldA
            int incC = 1; // step in rhs array
            double[] c = new double[NumColumns];
            rhsVector.CopyToArray(0, c, 0, NumRows);
            CBlas.Dtrsv(CBlasLayout.ColMajor, CBlasTriangular.Lower, CBlasTranspose.NoTranspose, CBlasDiagonal.NonUnit,
                NumRows, reflectorsAndL, 0, ldA, c, 0, incC);
            // TODO: Check output of BLAS somehow. E.g. Zero diagonal entries will result in NaN in the result vector.

            // Step 2: x = Q^T * c. Q is n-by-n, c is n-by-1, x is n-by-1

            // Prepare LAPACK arguments
            int m = NumColumns; // c has been grown past the limits of rhs.
            int nRhs = 1; // rhs = m-by-1
            int k = tau.Length;
            int ldC = m;

            // LAPACK query to find the optimum block size. TODO: Use a LAPACKE like wrapper to simplify this.
            int infoQuery = LapackUtilities.DefaultInfo;
            var workQuery = new double[1];
            int lWorkQuery = -1;
            Lapack.Dormlq("L", "T", m, nRhs, k, reflectorsAndL, 0, ldA, tau, 0, c, 0, ldC,
                workQuery, 0, lWorkQuery, ref infoQuery);
            CheckLapackExecution(infoQuery);

            // Call LAPACK to actually perform the operation
            int infoMult = LapackUtilities.DefaultInfo;
            int lWorkOptim = (int)(workQuery[0]);
            if (lWorkOptim < 1) lWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lWorkOptim];
            Lapack.Dormlq("L", "T", m, nRhs, k, reflectorsAndL, 0, ldA, tau, 0, c, 0, ldC,
                workOptim, 0, lWorkOptim, ref infoMult);
            CheckLapackExecution(infoMult);

            return Vector.CreateFromArray(c, false);
        }

        private static void CheckLapackExecution(int info)
        {
            if (info < 0) throw LapackUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
            else if (info > 0) throw new Exception(
                "This should not have happened. The corresponding LAPACK function does not return positive codes in info");
        }
    }
}
