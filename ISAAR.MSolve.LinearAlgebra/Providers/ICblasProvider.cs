//TODO: this and its implementations should be internal. The user should select a Provider that will then specify different BLAS,
//      SparseBLAS, LAPACK, etc. providers. The implementations should also be singletons or enums.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public interface ICblasProvider
    {
        #region BLAS Level 1
        /// <summary>
        /// y = alpha * x + y
        /// </summary>
        void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY);

        /// <summary>
        /// result = x * y
        /// </summary>
        double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY);

        /// <summary>
        /// result = sqrt(x * x)
        /// </summary>
        double Dnrm2(int n, double[] x, int offsetX, int incX);

        /// <summary>
        /// x = alpha * x
        /// </summary>
        void Dscal(int n, double alpha, double[] x, int offsetX, int incX);
        #endregion

        #region BLAS Level 2

        /// <summary>
        /// y = alpha * op(A) * x + beta * y, where op(A) = A or transpose(A). A is a general matrix, stored in full format.
        /// </summary>
        void Dgemv(CblasLayout layout, CblasTranspose trans, int m, int n, double alpha, double[] a, int offsetA, int ldA,
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);

        /// <summary>
        /// y = alpha * op(A) * x + beta * y, where op(A) = A or transpose(A). A is a symmetric matrix, stored in packed format.
        /// </summary>
        void Dspmv(CblasLayout layout, CblasTriangular uplo, int n, double alpha, double[] a, int offsetA,
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);

        /// <summary>
        /// x = op(A) * x, where op(A) = A or transpose(A). A is a triangular matrix, stored in packed format.
        /// </summary>
        void Dtpmv(CblasLayout layout, CblasTriangular uplo, CblasTranspose trans, CblasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX);

        /// <summary>
        /// x = inv(op(A)) * x, where op(A) = A or transpose(A). A is a triangular matrix, stored in packed format.
        /// </summary>
        void Dtpsv(CblasLayout layout, CblasTriangular uplo, CblasTranspose trans, CblasDiagonal diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX);
        #endregion

        #region BLAS Level 3
        /// <summary>
        /// C = alpha * op(A) * op(B) + beta * C, where op(A) = A or transpose(A), op(B) = B or transpose(B). A, B, C are general
        /// matrices, stored in full format.
        /// </summary>
        /// <param name="m">The number of rows of op(A) and C.</param>
        /// <param name="n">The number of columns of op(B) and C.</param>
        /// <param name="k">The number of columns of op(A), which must be equal to the number of rows of op(B).</param>
        void Dgemm(CblasLayout layout, CblasTranspose transA, CblasTranspose transB, int m, int n, int k, double alpha,
            double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC);
        #endregion

        #region BLAS-like extensions
        /// <summary>
        /// y = alpha * x + beta * y
        /// </summary>
        void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);
        #endregion
    }
}
