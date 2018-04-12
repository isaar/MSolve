using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;

namespace ISAAR.MSolve.Numerical.MKL
{
    // This works because the dlls are copied due to Compute.NET. If that library is removed or different dlls need to be used, I 
    // must handle that. Ideally, this should be done using a nuget package for MKL bindings (actually 1 for x64 and 1 for x86).
    public static class LAPACKE
    {
        public const int LAPACK_ROW_MAJOR = 101;
        public const int LAPACK_COL_MAJOR = 102;
        public const char LAPACK_NO_TRANSPOSE = 'N';
        public const char LAPACK_TRANSPOSE = 'T';
        public const char LAPACK_HERMITIAN_TRANSPOSE = 'C';

        /// <summary>
        /// Full matrix-vector multiplication. Covered by Compute.NET. I just implemented it for practice and testing purposes.
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "cblas_dgemv")]
        internal static extern void Dgemv(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA, int m, int n, 
            double alpha, double[] A, int ldA, double[] X, int incX, double beta, double[] Y, int incY ); // CBLAS, rather than LAPACKE

        /// <summary>
        /// QR factorization of a full <paramref name="m"/>-by-<paramref name="n"/> matrix. LAPACKE interface. For the Fortran 77
        /// interface see https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgels.htm,
        /// http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga3766ea903391b5cf9008132f7440ec7b.html#ga3766ea903391b5cf9008132f7440ec7b
        /// </summary>
        /// <param name="matrixLayout">Specifies whether the matrix storage layout is row major (<see cref="LAPACK_ROW_MAJOR"/>) 
        ///     or column major (<see cref="LAPACK_COL_MAJOR"/>).</param>
        /// <param name="m">The number of rows in the matrix A (<paramref name="m"/> ≥ 0).</param>
        /// <param name="n">The number of columns in the matrix A (<paramref name="n"/> ≥ 0).</param>
        /// <param name="A"> Contains the original matrix. Array of size max(1, <paramref name="ldA"/> * <paramref name="n"/>)
        ///     for column major layout and max(1, <paramref name="ldA"/> * <paramref name="n"/>) for row major layout. 
        ///     Overwritten by the factorization data as follows: The elements on and above the diagonal contain the 
        ///     min(<paramref name="m"/>, <paramref name="n"/>)-by-<paramref name="n"/> upper trapezoidal matrix R (R is upper 
        ///     triangular if <paramref name="m"/> ≥ <paramref name="n"/>). The elements below the diagonal, with the array 
        ///     <paramref name="Tau"/>, present the orthogonal matrix Q as a product of min(<paramref name="m"/>, 
        ///     <paramref name="n"/>) elementary reflectors.</param>
        /// <param name="ldA">The leading dimension of <paramref name="A"/>; at least max(1, <paramref name="m"/>) for column 
        ///     major layout and max(1, <paramref name="n"/>) for row major layout.</param>
        /// <param name="Tau"> Array of size at least max(1, min(<paramref name="m"/>, <paramref name="n"/>)). Contains scalars 
        ///     that define elementary reflectors for the matrix Q in its decomposition in a product of elementary reflectors.
        ///     </param>
        /// <returns></returns>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dgeqrf")]
        internal static extern int Dgeqrf(int matrixLayout, int m, int n, double[] A, int ldA, double[] Tau);

        /// <summary>
        /// LU factorization of full matrix. Covered by Compute.NET. I just implemented it for practice and testing purposes.
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dgetrf")]
        internal static extern int Dgetrf(int matrixLayout, int m, int n, double[] A, int ldA, int[] iPiv);

        /// <summary>
        /// Linear system solution using the full LU factorization from <see cref="Dgetrf(int, int, int, double[], int, int[])"/>.
        /// Covered by Compute.NET. I just implemented it for practice and testing purposes.
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dgetrs")]
        internal static extern int Dgetrs(int matrixLayout, char transA, int n, int nRhs, 
            double[] A, int ldA, int[] iPiv, double[] B, int ldB);

        /// <summary>
        /// Generates the orthogonal matrix Q from the QR factorization created by 
        /// <see cref="Dgeqrf(int, int, int, double[], int, double[])"/>.
        /// </summary>
        /// <param name="matrixLayout">Specifies whether the matrix storage layout is row major (<see cref="LAPACK_ROW_MAJOR"/>) 
        ///     or column major (<see cref="LAPACK_COL_MAJOR"/>).</param>
        /// <param name="m">The order of the orthogonal matrix Q (<paramref name="m"/> ≥ 0).</param>
        /// <param name="n">The number of columns of Q to be computed (0 ≤ <paramref name="n"/> ≤ <paramref name="m"/>). 
        ///     To generate the whole matrix Q, use <paramref name="n"/> = <paramref name="m"/>.</param>
        /// <param name="k">The number of elementary reflectors whose product defines the matrix Q (0 ≤ <paramref name="k"/>
        ///     ≤ <paramref name="n"/>). To generate the whole matrix Q, use: <paramref name="k"/> = p, where p = number of
        ///     columns of the original matrix (if <paramref name="m"/> ≥ p).</param>
        /// <param name="A">Array returned by <see cref="Dgeqrf(int, int, int, double[], int, double[])"/>. The elements on and 
        ///     above the diagonal contain the min(<paramref name="m"/>,p)-by-p upper trapezoidal matrix R (R is upper triangular 
        ///     if <paramref name="m"/> ≥ p) The elements below the diagonal, with the array <paramref name="Tau"/>, present the 
        ///     orthogonal matrix Q as a product of min(<paramref name="m"/>, p) elementary reflectors. The size of 
        ///     <paramref name="A"/> is max(1, <paramref name="ldA"/> * <paramref name="n"/>) for column major layout and 
        ///     max(1, <paramref name="ldA"/> * <paramref name="m"/>) for row major layout.</param>
        /// <param name="ldA">The leading dimension of <paramref name="A"/>; at least max(1, <paramref name="m"/>) for column 
        ///     major layout and max(1, <paramref name="n"/>) for row major layout.</param>
        /// <param name="Tau">Array returned by <see cref="Dgeqrf(int, int, int, double[], int, double[])"/>. Contains scalars 
        ///     that define elementary reflectors for the matrix Q in its decomposition in a product of elementary reflectors.
        ///     The size of <paramref name="Tau"/> must be at least max(1, <paramref name="k"/>).</param>
        /// <returns></returns>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dorgqr")]
        internal static extern int Dorgqr(int matrixLayout, int m, int n, int k, double[] A, int ldA, double[] Tau);

    }
}
