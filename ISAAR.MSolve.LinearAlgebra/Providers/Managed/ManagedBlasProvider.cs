using System;
using DotNumerics.LinearAlgebra.CSLapack;

//TODO: find a managed BLAS that supports the methods DotNumerics doesn't.
namespace ISAAR.MSolve.LinearAlgebra.Providers.Managed
{
    /// <summary>
    /// Provides managed C# implementations of the linear algebra operations defined by <see cref="IBlasProvider"/>. Uses the 
    /// library DotNumerics (see http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSBlas/Default.aspx) for the most
    /// part. For BLAS subroutines not provided by DotNumerics, custom C# implementations are used instead. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ManagedBlasProvider : IBlasProvider
    {
        //TODO: perhaps these should not be static.
        private static readonly DAXPY daxpy = new DAXPY();
        private static readonly DDOT ddot = new DDOT();
        private static readonly DGEMM dgemm = new DGEMM();
        private static readonly DGEMV dgemv = new DGEMV();
        private static readonly DNRM2 dnrm2 = new DNRM2();
        private static readonly DSCAL dscal = new DSCAL();
        private static readonly DTRSV dtrsv = new DTRSV();

        internal static ManagedBlasProvider UniqueInstance { get; } = new ManagedBlasProvider();

        private ManagedBlasProvider() { } // private constructor for singleton pattern

        #region BLAS Level 1
        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/daxpy.aspx
        /// </summary>
        public void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => daxpy.Run(n, alpha, x, offsetX, incX, ref y, offsetY, incY);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/ddot.aspx
        /// </summary>
        public double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY)
            => ddot.Run(n, x, offsetX, incX, y, offsetY, incY);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dnrm2.aspx
        /// </summary>
        public double Dnrm2(int n, double[] x, int offsetX, int incX)
            => dnrm2.Run(n, x, offsetX, incX);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dscal.aspx
        /// </summary>
        public void Dscal(int n, double alpha, double[] x, int offsetX, int incX)
            => dscal.Run(n, alpha, ref x, offsetX, incX);
        #endregion

        #region BLAS Level 2

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgemv.aspx
        /// </summary>
        public void Dgemv(TransposeMatrix transA, int m, int n,
            double alpha, double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX,
            double beta, double[] y, int offsetY, int incY)
            => dgemv.Run(transA.Translate(), m, n, alpha, a, offsetA, ldA, x, offsetX, incX, beta, ref y, offsetY, incY);

        public void Dspmv(StoredTriangle uplo, int n,
            double alpha, double[] a, int offsetA, double[] x, int offsetX, int incX,
            double beta, double[] y, int offsetY, int incY)
        {
            // y = alpha * L * x + beta * y 
            CblasLevel2Implementations.LowerTimesVectorPackedRowMajor(
                CblasLevel2Implementations.Diagonal.Regular, n, alpha, a, offsetA, x, offsetX, incX, beta, y, offsetY, incY);

            // y = alpha * U * x + y, where U has 0 diagonal
            CblasLevel2Implementations.UpperTimesVectorPackedColMajor(
                CblasLevel2Implementations.Diagonal.Zero, n, alpha, a, offsetA, x, offsetX, incX, 1.0, y, offsetY, incY);
        }

        public void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
        {
            // The copy may be avoidable in trangular operations, if we start the dot products from the bottom
            var input = new double[x.Length];
            Array.Copy(x, input, x.Length);

            CblasLevel2Implementations.Diagonal managedDiag = (diag == DiagonalValues.NonUnit) ?
                CblasLevel2Implementations.Diagonal.Regular : CblasLevel2Implementations.Diagonal.Unit;
            if (UseUpperImplementation(uplo, transA))
            {
                CblasLevel2Implementations.UpperTimesVectorPackedColMajor(
                    managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
            }
            else
            {
                CblasLevel2Implementations.LowerTimesVectorPackedRowMajor(
                    managedDiag, n, 1.0, a, offsetA, input, offsetX, incX, 0.0, x, offsetX, incX);
            }
        }

        public void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX)
        {
            bool unit = (diag == DiagonalValues.Unit) ? true : false;
            if (UseUpperImplementation(uplo, transA))
            {
                CblasLevel2Implementations.BackSubstitutionPackedColMajor(unit, n, a, offsetA, x, offsetX, incX);
            }
            else CblasLevel2Implementations.ForwardSubstitutionPackedRowMajor(unit, n, a, offsetA, x, offsetX, incX);
        }

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dtrsv.aspx
        /// </summary>
        public void Dtrsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX)
            => dtrsv.Run(uplo.Translate(), transA.Translate(), diag.Translate(), n, a, offsetA, ldA, ref x, offsetX, incX);
        #endregion

        #region BLAS Level 3

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgemm.aspx
        /// </summary>
        public void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha,
            double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC)
            => dgemm.Run(transA.Translate(), transB.Translate(), m, n, k, alpha, a, offsetA, ldA, b, offsetB, ldB,
                beta, ref c, offsetC, ldC);
        #endregion

        private static bool UseUpperImplementation(StoredTriangle uplo, TransposeMatrix transA)
        {
            //if (transA == TransposeMatrix.ConjugateTranspose)
            //    throw new ArgumentException("Cannot use conjugate transpose operations for double matrices and vectors.");

            if (uplo == StoredTriangle.Upper && transA == TransposeMatrix.NoTranspose) return true;
            if (uplo == StoredTriangle.Upper && transA == TransposeMatrix.Transpose) return false;
            if (uplo == StoredTriangle.Lower && transA == TransposeMatrix.NoTranspose) return true;
            if (uplo == StoredTriangle.Lower && transA == TransposeMatrix.Transpose) return false;
            throw new Exception("This code should not have been reached");
        }
    }
}
