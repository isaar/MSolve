using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;

// TODO: Perhaps I should use the QR with pivoting
// TODO: Add the option to specify if the diagonal entries of A are non diagonal.
// TODO: Handle the case when m < n. That would be a RQ factorization. Also update the documentation when necessary.
// TODO: I think least squares only work if the columns of A are independent. This is not taken care of so far.
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    /// <summary>
    /// QR factorization of a matrix A: A = Q * R, where A is an m-by-n matrix with m &gt;= n, Q is an m-by-m orthogonal matrix
    /// and R is an m-by-n upper trapezoidal matrix. Uses LAPACK.
    /// For more see:
    /// https://software.intel.com/en-us/mkl-developer-reference-c-orthogonal-factorizations-lapack-computational-routines#E832D468-0891-40EC-9468- ,
    /// http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga3766ea903391b5cf9008132f7440ec7b.html#ga3766ea903391b5cf9008132f7440ec7b.
    /// Authors: Serafeim Bakalakos
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
        /// Calculates the QR factorization of a matrix, such that A = Q * R. Requires an extra 
        /// min(<paramref name="numRows"/>, <paramref name="numCols"/>) available memory.
        /// </summary>
        /// <param name="numRows">The number of rows of the original matrix.</param>
        /// <param name="numCols">The number of columns of the original matrix.</param>
        /// <param name="matrix">The internal buffer storing the matrix entries in column major layout. It will 
        ///     be overwritten with the factorization data.</param>
        /// <exception cref="NotImplementedException">Thrown if <paramref name="numCols"/> &gt; <paramref name="numRows"/>.
        ///     </exception>
        /// <exception cref="Exceptions.LapackException">Thrown if tha call to LAPACK fails due to an invalid 
        ///     <paramref name="matrix"/>.</exception>
        public static QRFactorization Factorize(int numRows, int numCols, double[] matrix)
        {
            if (numRows < numCols)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }

            // Prepare LAPACK arguments
            int m = numRows;
            int n = numCols;
            int ldA = m;
            int minDim = Math.Min(m, n); //TODO: this is known to be numCols (though it may change in the future)
            double[] tau = new double[minDim];

            // LAPACK query to find the optimum block size. TODO: Use a LAPACKE like wrapper to simplify this.
            int infoQuery = LapackUtilities.DefaultInfo;
            var workQuery = new double[1];
            int lWorkQuery = -1;
            Lapack.Dgeqrf(m, n, matrix, 0, ldA, tau, 0, workQuery, 0, lWorkQuery, ref infoQuery);
            CheckLapackExecution(infoQuery);

            // Call LAPACK to actually perform the operation
            int lWorkOptim = (int)(workQuery[0]);
            if (lWorkOptim < 1) lWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lWorkOptim];
            int info = LapackUtilities.DefaultInfo;
            Lapack.Dgeqrf(m, n, matrix, 0, ldA, tau, 0, workOptim, 0, lWorkOptim, ref info);
            CheckLapackExecution(info);

            return new QRFactorization(numRows, numCols, matrix, tau);
        }

        /// <summary>
        /// Explicitly creates the matrix Q1, which consists of the first n columns of the orthogonal matrix Q that resulted  
        /// from the factorization: A =  Q * R = [Q1, Q2] * [R1; 0] (Matlab notation) = Q1 * R1, 
        /// where A is m-by-n, Q is m-by-m, R is m-by-n, Q1 is m-by-n and R1 is (n-by-n). 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetEconomyFactorQ()
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }
            int m = NumRows;
            int n = NumColumns;
            int k = NumColumns;
            int ldA = m;

            var q = new double[NumRows * NumColumns];
            Array.Copy(reflectorsAndR, q, q.Length);

            // LAPACK query to find the optimum block size. TODO: Use a LAPACKE like wrapper to simplify this.
            int infoQuery = LapackUtilities.DefaultInfo;
            var workQuery = new double[1];
            int lWorkQuery = -1;
            Lapack.Dorgqr(m, n, k, q, 0, ldA, tau, 0, workQuery, 0, lWorkQuery, ref infoQuery);
            CheckLapackExecution(infoQuery);

            // Call LAPACK to actually perform the operation
            int lWorkOptim = (int)(workQuery[0]);
            if (lWorkOptim < 1) lWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lWorkOptim];
            int info = LapackUtilities.DefaultInfo;
            Lapack.Dorgqr(m, n, k, q, 0, ldA, tau, 0, workOptim, 0, lWorkOptim, ref info);
            CheckLapackExecution(info);

            return Matrix.CreateFromArray(q, NumRows, NumColumns, false);
        }

        /// <summary>
        /// Explicitly creates the upper triangular matrix R1, which consists of the first n rows of the matrix R that resulted  
        /// from the factorization: A =  Q * R = [Q1, Q2] * [R1; 0] (Matlab notation) = Q1 * R1, 
        /// where A is m-by-n, Q is m-by-m, R is m-by-n, Q1 is m-by-n and R1 is (n-by-n). 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public TriangularUpper GetEconomyFactorR()
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }
            double[] r = Conversions.RectColMajorToSquarePackedUpperColMajor(NumRows, NumColumns, reflectorsAndR);
            return TriangularUpper.CreateFromArray(NumColumns, r, false);
        }

        /// <summary>
        /// Explicitly creates the orthogonal matrix Q that resulted from the factorization: A = Q * R, where A is m-by-n, 
        /// Q is m-by-m and R is m-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetFactorQ()
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }

            // Prepare LAPACK arguments
            int m = NumRows;
            int n = m;
            int k = NumColumns;
            int ldA = m;

            // We need a larger buffer for Q (m-by-m) > reflectorsAndR (m-by-p)
            double[] q = ArrayColMajor.ResizeCols(NumRows, NumColumns, NumRows, reflectorsAndR);

            // LAPACK query to find the optimum block size. TODO: Use a LAPACKE like wrapper to simplify this.
            int infoQuery = LapackUtilities.DefaultInfo;
            var workQuery = new double[1];
            int lWorkQuery = -1;
            Lapack.Dorgqr(m, n, k, q, 0, ldA, tau, 0, workQuery, 0, lWorkQuery, ref infoQuery);
            CheckLapackExecution(infoQuery);

            // Call LAPACK to actually perform the operation
            int lWorkOptim = (int)(workQuery[0]);
            if (lWorkOptim < 1) lWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lWorkOptim];
            int info = LapackUtilities.DefaultInfo;
            Lapack.Dorgqr(m, n, k, q, 0, ldA, tau, 0, workOptim, 0, lWorkOptim, ref info);
            CheckLapackExecution(info);

            return Matrix.CreateFromArray(q, NumRows, NumRows, false);
        }

        /// <summary>
        /// Explicitly creates the upper trapezoidal matrix R that resulted from the factorization: A = Q * R, where A is m-by-n, 
        /// Q is m-by-m and R is m-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetFactorR()
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }
            double[] r = Conversions.RectColMajorToRectUpperColMajor(NumRows, NumColumns, reflectorsAndR);
            return Matrix.CreateFromArray(r , NumRows, NumColumns, false);
        }

        /// <summary>
        /// Solve the linear least squares problem A * x = b => x = inv(R) * transpose(Q) * b.
        /// Warning: the columns of the original matrix A must be independent for this to work.
        /// </summary>
        /// <param name="rhsVector">The right hand side vector b. It may lie outside the column space of the original matrix. Its 
        ///     <see cref="IIndexable1D.Length"/> must be equal to this.<see cref="NumRows"/>.</param>
        /// <exception cref="Exceptions.LapackException">Thrown if the call to LAPACK fails due to <paramref name="rhsVector"/> 
        ///     having a different <see cref="IIndexable1D.Length"/> than this.<see cref="NumRows"/>.</exception>
        public Vector SolveLeastSquares(Vector rhsVector) //TODO: perhaps I should use the driver routines of LAPACKE
        {
            if (NumRows < NumColumns)
            {
                throw new NotImplementedException("For now, the number of rows must be >= the number of columns");
            }
            Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);

            // Least squares: x = inv(A^T * A) * A^T * b = inv(R) * Q^T * b, where b is the right hand side vector. 
            // Step 1: c = Q^T * b. Q is m-by-m, b is m-by-1 => c is m-by-1

            // Prepare LAPACK arguments
            double[] c = rhsVector.CopyToArray();
            int m = rhsVector.Length;
            int nRhs = 1; // rhs = m-by-1
            int k = tau.Length;
            int ldA = m;
            int ldC = m;

            // LAPACK query to find the optimum block size. TODO: Use a LAPACKE like wrapper to simplify this.
            int infoQuery = LapackUtilities.DefaultInfo;
            var workQuery = new double[1];
            int lWorkQuery = -1;
            Lapack.Dormqr("L", "T", m, nRhs, k, reflectorsAndR, 0, ldA, tau, 0, c, 0, ldC, 
                workQuery, 0, lWorkQuery, ref infoQuery);
            CheckLapackExecution(infoQuery);

            // Call LAPACK to actually perform the operation
            int infoMult = LapackUtilities.DefaultInfo;
            int lWorkOptim = (int)(workQuery[0]);
            if (lWorkOptim < 1) lWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lWorkOptim];
            Lapack.Dormqr("L", "T", m, nRhs, k, reflectorsAndR, 0, ldA, tau, 0, c, 0, ldC,
                workOptim, 0, lWorkOptim, ref infoMult);
            CheckLapackExecution(infoMult);

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
            CBlas.Dtrsv(CBlasLayout.ColMajor, CBlasTriangular.Upper, CBlasTranspose.NoTranspose, CBlasDiagonal.NonUnit,
                n, reflectorsAndR, 0, ldR, c, 0, incC);
            // TODO: Check output of BLAS somehow. E.g. Zero diagonal entries will result in NaN in the result vector.
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
