//TODO: I could also provide a CBLAS like wrapper. However it is not as high a priority, since BLAS is easy to use. If I do so,
//      the enums should be transfered to the wrapper. Here I would need to use strings only.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Provides linear algebra operations as defined by the BLAS (Basic Linear Algebra Subroutines) interface. These operations
    /// concern vector-vector operations (level 1), matrix-vector operations (level 2) and matrix-matrix operations (level 3). 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal interface IBlasProvider
    {
        #region BLAS Level 1

        /// <summary>
        /// y = alpha * x + y. See 
        /// http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga8f99d6a644d3396aa32db472e0cfc91c.html#ga8f99d6a644d3396aa32db472e0cfc91c
        /// </summary>
        void Daxpy(int n, double alpha, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY);

        /// <summary>
        /// result = x * y. See
        /// http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html#ga75066c4825cb6ff1c8ec4403ef8c843a
        /// </summary>
        double Ddot(int n, double[] x, int offsetX, int incX, double[] y, int offsetY, int incY);

        /// <summary>
        /// result = sqrt(x * x). See
        /// http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga4b25c539b862414e6f91ebb30b391d7c.html#ga4b25c539b862414e6f91ebb30b391d7c
        /// </summary>
        double Dnrm2(int n, double[] x, int offsetX, int incX);

        /// <summary>
        /// x = alpha * x. See
        /// http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga793bdd0739bbd0e0ec8655a0df08981a.html#ga793bdd0739bbd0e0ec8655a0df08981a
        /// </summary>
        void Dscal(int n, double alpha, double[] x, int offsetX, int incX);
        #endregion

        #region BLAS Level 2

        /// <summary>
        /// y = alpha * op(A) * x + beta * y, where op(A) = A or transpose(A). A is a general matrix, stored in full format. See
        /// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9
        /// </summary>
        void Dgemv(TransposeMatrix transA, int m, int n, double alpha, double[] a, int offsetA, int ldA,
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);

        /// <summary>
        /// y = alpha * op(A) * x + beta * y, where op(A) = A or transpose(A). A is a symmetric matrix, stored in packed format.
        /// See
        /// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gab746575c4f7dd4eec72e8110d42cefe9.html#gab746575c4f7dd4eec72e8110d42cefe9
        /// </summary>
        void Dspmv(StoredTriangle uplo, int n, double alpha, double[] a, int offsetA,
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);

        /// <summary>
        /// x = op(A) * x, where op(A) = A or transpose(A). A is a triangular matrix, stored in packed format.
        /// See
        /// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga1d9a8ecfddfea2c84e73e28e1ebb74cf.html#ga1d9a8ecfddfea2c84e73e28e1ebb74cf
        /// </summary>
        void Dtpmv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX);

        /// <summary>
        /// x = inv(op(A)) * x, where op(A) = A or transpose(A). A is a triangular matrix, stored in packed format. See
        /// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0fff73e765a7655a67da779c898863f1.html#ga0fff73e765a7655a67da779c898863f1
        /// </summary>
        void Dtpsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, double[] x, int offsetX, int incX);

        /// <summary>
        /// x = inv(op(A)) * x, where op(A) = A or transpose(A). A is a triangular matrix, stored in full format, but only the 
        /// entries corresponding to the designated triangle will be accessed. See
        /// http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_gad2a01dd62718b28e35b752dbad8474ab.html#gad2a01dd62718b28e35b752dbad8474ab
        /// </summary>
        void Dtrsv(StoredTriangle uplo, TransposeMatrix transA, DiagonalValues diag, int n,
            double[] a, int offsetA, int ldA, double[] x, int offsetX, int incX);
        #endregion

        #region BLAS Level 3

        /// <summary>
        /// C = alpha * op(A) * op(B) + beta * C, where op(A) = A or transpose(A), op(B) = B or transpose(B). A, B, C are general
        /// matrices, stored in full format. See
        /// http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html#gaeda3cbd99c8fb834a60a6412878226e1
        /// </summary>
        /// <param name="m">The number of rows of op(A) and C.</param>
        /// <param name="n">The number of columns of op(B) and C.</param>
        /// <param name="k">The number of columns of op(A), which must be equal to the number of rows of op(B).</param>
        void Dgemm(TransposeMatrix transA, TransposeMatrix transB, int m, int n, int k, double alpha,
            double[] a, int offsetA, int ldA, double[] b, int offsetB, int ldB, double beta, double[] c, int offsetC, int ldC);
        #endregion
    }
}
