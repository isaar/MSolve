using System.Runtime.InteropServices;
using IntelMKL.LP64;

//TODO: Replace Compute.NET functions with th LAPACKE interface for triangulations (~trf) and system solutions (~trs)
//TODO: Wrap these methods with ones using enums, which will then call MKL with the correct arguments.
//TODO: This class works because the dlls are copied due to Compute.NET. If that library is removed or different dlls need to 
//      be used, I must handle that. Ideally, this should be done using a nuget package for MKL bindings (actually 1 for x64 
//      and 1 for x86).
namespace ISAAR.MSolve.LinearAlgebra.Providers.MKL
{
    /// <summary>
    /// Platform Invoke methods for Intel MKL's LAPACKE (C interface of LAPACK). These are not covered by any nuget packages.
    /// Also see the included "lapacke.h" C header file and MKL's C user guide.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class LapackePInvokes
    {
        public const int LAPACK_ROW_MAJOR = 101;
        public const int LAPACK_COL_MAJOR = 102;
        public const char LAPACK_NO_TRANSPOSE = 'N';
        public const char LAPACK_TRANSPOSE = 'T';
        public const char LAPACK_HERMITIAN_TRANSPOSE = 'C';
        public const char LAPACK_SIDE_LEFT = 'L';
        public const char LAPACK_SIDE_RIGHT = 'R';
        public const char LAPACK_UPPER = 'U';
        public const char LAPACK_LOWER = 'L';

        /// <summary>
        /// LQ factorization of a full <paramref name="m"/>-by-<paramref name="n"/> matrix.
        /// </summary>
        /// <param name="matrixLayout"></param>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <param name="A"></param>
        /// <param name="ldA"></param>
        /// <param name="Tau"></param>
        /// <returns></returns>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dgelqf")]
        internal static extern int Dgelqf(int matrixLayout, int m, int n, double[] A, int ldA, double[] Tau);

        /// <summary>
        /// Full matrix-vector multiplication. Covered by Compute.NET. I just implemented it for practice and testing purposes.
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "cblas_dgemv")]
        public static extern void Dgemv(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA, int m, int n, 
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
        /// <param name="A">Array containing the original matrix, with size: max(1, <paramref name="ldA"/> * <paramref name="n"/>)
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
        public static extern int Dgetrf(int matrixLayout, int m, int n, double[] A, int ldA, int[] iPiv);

        /// <summary>
        /// Inversion of matrix in full format, that underwent LU factorization by 
        /// <see cref="Lapack.Dgetrf(ref int, ref int, ref double, ref int, ref int, ref int)"/>. 
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dgetri")]
        internal static extern int Dgetri(int matrixLayout, int n, double[] A, int ldA, int[] iPiv);

        /// <summary>
        /// Linear system solution using the full LU factorization from <see cref="Dgetrf(int, int, int, double[], int, int[])"/>.
        /// Covered by Compute.NET. I just implemented it for practice and testing purposes.
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dgetrs")]
        internal static extern int Dgetrs(int matrixLayout, char transA, int n, int nRhs, 
            double[] A, int ldA, int[] iPiv, double[] B, int ldB);

        /// <summary>
        /// Generates the orthogonal matrix Q from the LQ factorization created by 
        /// <see cref="Dgelqf(int, int, int, double[], int, double[])"/>.
        /// </summary>
        /// <param name="matrixLayout"></param>
        /// <param name="n">The order of the orthogonal matrix Q (<paramref name="n"/> ≥ 0). WARNING: equal to the number of 
        ///     columns of the original matrix, rather than the number of rows.</param>
        /// <param name="m">The number of rows of Q to be computed (0 ≤ <paramref name="m"/> ≤ <paramref name="n"/>). 
        ///     To generate the whole matrix Q, use <paramref name="m"/> = <paramref name="n"/>.</param>
        /// <param name="k">The number of elementary reflectors whose product defines the matrix Q (0 ≤ <paramref name="k"/>
        ///     ≤ <paramref name="m"/>). To generate the whole matrix Q, use: <paramref name="k"/> = p, where p = number of
        ///     rows of the original matrix (if <paramref name="n"/> ≥ p).</param>
        /// <param name="A">WARNING: To get the whole matrix Q, you must allocate enough space = 
        ///     <paramref name="n"/>-by-<paramref name="n"/>, although the factorized data only need 
        ///     p-by-<paramref name="n"/>.</param>
        /// <param name="ldA">The leading dimension of <paramref name="A"/>; at least max(1, <paramref name="m"/>) for column 
        ///     major layout and max(1, <paramref name="n"/>) for row major layout.</param>
        /// <param name="Tau">Array returned by <see cref="Dgelqf(int, int, int, double[], int, double[])"/>. Contains scalars 
        ///     that define elementary reflectors for the matrix Q in its decomposition in a product of elementary reflectors.
        ///     The size of <paramref name="Tau"/> must be at least max(1, <paramref name="k"/>)</param>
        /// <returns></returns>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dorglq")]
        internal static extern int Dorglq(int matrixLayout, int m, int n, int k, double[] A, int ldA, double[] Tau);

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
        ///     max(1, <paramref name="ldA"/> * <paramref name="m"/>) for row major layout. WARNING: To get the whole matrix Q, 
        ///     you must allocate enough space = <paramref name="m"/>-by-<paramref name="m"/>, although the factorized data only
        ///     need <paramref name="m"/>-by-p.</param>
        /// <param name="ldA">The leading dimension of <paramref name="A"/>; at least max(1, <paramref name="m"/>) for column 
        ///     major layout and max(1, <paramref name="n"/>) for row major layout.</param>
        /// <param name="Tau">Array returned by <see cref="Dgeqrf(int, int, int, double[], int, double[])"/>. Contains scalars 
        ///     that define elementary reflectors for the matrix Q in its decomposition in a product of elementary reflectors.
        ///     The size of <paramref name="Tau"/> must be at least max(1, <paramref name="k"/>).</param>
        /// <returns></returns>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dorgqr")]
        internal static extern int Dorgqr(int matrixLayout, int m, int n, int k, double[] A, int ldA, double[] Tau);

        /// <summary>
        /// Multiply an <paramref name="m"/>-by-<paramref name="n"/> (or vector) by the orthogonal matrix Q created by 
        /// <see cref="Dgelqf(int, int, int, double[], int, double[])"/> or its transpose.
        /// </summary>
        /// <param name="matrixLayout">Specifies whether the matrix storage layout is row major (<see cref="LAPACK_ROW_MAJOR"/>) 
        ///     or column major (<see cref="LAPACK_COL_MAJOR"/>).</param>
        /// <param name="side"> Must be either <see cref="LAPACK_SIDE_LEFT"/> or <see cref="LAPACK_SIDE_RIGHT"/>. If 
        ///     <paramref name="side"/> = <see cref="LAPACK_SIDE_LEFT"/>, Q or Q^T is multiplied to <paramref name="C"/> from the
        ///     left. If <paramref name="side"/> = <see cref="LAPACK_SIDE_RIGHT"/>, Q or Q^T is multiplied to 
        ///     <paramref name="C"/> from the right.</param>
        /// <param name="transQ">Must be either <see cref="LAPACK_NO_TRANSPOSE"/> or <see cref="LAPACK_TRANSPOSE"/>. If 
        ///     <paramref name="transQ"/> = <see cref="LAPACK_NO_TRANSPOSE"/>, the routine multiplies <paramref name="C"/> by Q. 
        ///     If trans='T'<paramref name="transQ"/> = <see cref="LAPACK_TRANSPOSE"/>, the routine multiplies 
        ///     <paramref name="C"/> by Q^T.</param>
        /// <param name="m">The number of rows in the matrix <paramref name="C"/> (<paramref name="m"/> ≥ 0).</param>
        /// <param name="n">The number of columns in the matrix <paramref name="C"/> (<paramref name="n"/> ≥ 0).</param>
        /// <param name="k">The number of elementary reflectors whose product defines the matrix Q. Constraints:        
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_LEFT"/> : 0 ≤ <paramref name="k"/> ≤ <paramref name="m"/>.
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_RIGHT"/> : 0 ≤ <paramref name="k"/> ≤ <paramref name="n"/>.</param>
        /// <param name="A"></param>
        /// <param name="ldA">The leading dimension of <paramref name="A"/>. Constraints:
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_LEFT"/> : <paramref name="ldA"/> ≥ 
        ///     max(1, <paramref name="k"/>) for column major layout and max(1, <paramref name="m"/>) for row major layout.
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_RIGHT"/> : <paramref name="ldA"/> ≥ 
        ///     max(1, <paramref name="k"/>) for column major layout and max(1, <paramref name="n"/>) for row major layout.
        ///     </param>
        /// <param name="Tau">Array returned by <see cref="Dgelqf(int, int, int, double[], int, double[])"/>. Contains scalars 
        ///     that define elementary reflectors for the matrix Q in its decomposition in a product of elementary reflectors. 
        ///     The size of <paramref name="Tau"/> must be at least max(1, <paramref name="k"/>).</param>
        /// <param name="C">Array containing the <paramref name="m"/>-by<paramref name="n"/> right hand side matrix. Its size is
        ///     max(1, <paramref name="ldC"/> * <paramref name="n"/>) for column major layout and 
        ///     max(1, <paramref name="ldC"/> * <paramref name="m"/>) for row major layout. Overwritten by the product Q * C, 
        ///     Q^T * C, C * Q, or C * Q^T (specified by <paramref name="side"/> and <paramref name="transQ"/>).</param>
        /// <param name="ldC">The leading dimension of <paramref name="C"/>. Constraints:  <paramref name="ldC"/> ≥
        ///     max(1, <paramref name="m"/>) for column major layout and max(1, <paramref name="n"/>) for row major layout.</param>
        /// <returns></returns>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dormlq")]
        internal static extern int Dormlq(int matrixLayout, char side, char transQ,
            int m, int n, int k, double[] A, int ldA, double[] Tau, double[] C, int ldC);

        /// <summary>
        /// Multiply an <paramref name="m"/>-by-<paramref name="n"/> matrix (or vector) by the orthogonal matrix Q created by 
        /// <see cref="Dgeqrf(int, int, int, double[], int, double[])"/> or its transpose.
        /// </summary>
        /// <param name="matrixLayout">Specifies whether the matrix storage layout is row major (<see cref="LAPACK_ROW_MAJOR"/>) 
        ///     or column major (<see cref="LAPACK_COL_MAJOR"/>).</param>
        /// <param name="side"> Must be either <see cref="LAPACK_SIDE_LEFT"/> or <see cref="LAPACK_SIDE_RIGHT"/>. If 
        ///     <paramref name="side"/> = <see cref="LAPACK_SIDE_LEFT"/>, Q or Q^T is multiplied to <paramref name="C"/> from the
        ///     left. If <paramref name="side"/> = <see cref="LAPACK_SIDE_RIGHT"/>, Q or Q^T is multiplied to 
        ///     <paramref name="C"/> from the right.</param>
        /// <param name="transQ">Must be either <see cref="LAPACK_NO_TRANSPOSE"/> or <see cref="LAPACK_TRANSPOSE"/>. If 
        ///     <paramref name="transQ"/> = <see cref="LAPACK_NO_TRANSPOSE"/>, the routine multiplies <paramref name="C"/> by Q. 
        ///     If trans='T'<paramref name="transQ"/> = <see cref="LAPACK_TRANSPOSE"/>, the routine multiplies 
        ///     <paramref name="C"/> by Q^T.</param>
        /// <param name="m">The number of rows in the matrix <paramref name="C"/> (<paramref name="m"/> ≥ 0).</param>
        /// <param name="n">The number of columns in the matrix <paramref name="C"/> (<paramref name="n"/> ≥ 0).</param>
        /// <param name="k">The number of elementary reflectors whose product defines the matrix Q. Constraints:        
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_LEFT"/> : 0 ≤ <paramref name="k"/> ≤ <paramref name="m"/>.
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_RIGHT"/> : 0 ≤ <paramref name="k"/> ≤ <paramref name="n"/>.</param>
        /// <param name="A">Array returned by <see cref="Dgeqrf(int, int, int, double[], int, double[])"/>. The elements on and 
        ///     above the diagonal contain upper trapezoidal matrix R. The elements below the diagonal, with the array 
        ///     <paramref name="Tau"/>, present the orthogonal matrix Q as a product of elementary reflectors. The size of 
        ///     <paramref name="A"/> is max(1, <paramref name="ldA"/> * <paramref name="k"/>) for column major layout, 
        ///     max(1, <paramref name="ldA"/> * <paramref name="m"/>) for row major layout and <paramref name="side"/> = 
        ///     <see cref="LAPACK_SIDE_LEFT"/> and max(1, <paramref name="ldA"/> * <paramref name="n"/>) for row major layout and
        ///     <paramref name="side"/> = <see cref="LAPACK_SIDE_RIGHT"/>.</param>
        /// <param name="ldA">The leading dimension of <paramref name="A"/>. Constraints:
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_LEFT"/> : <paramref name="ldA"/> ≥ 
        ///     max(1, <paramref name="m"/>) for column major layout and max(1, <paramref name="k"/>) for row major layout.
        ///     If <paramref name="side"/> = <see cref="LAPACK_SIDE_RIGHT"/> : <paramref name="ldA"/> ≥ 
        ///     max(1, <paramref name="n"/>) for column major layout and max(1, <paramref name="k"/>) for row major layout.
        ///     </param>
        /// <param name="Tau">Array returned by <see cref="Dgeqrf(int, int, int, double[], int, double[])"/>. Contains scalars 
        ///     that define elementary reflectors for the matrix Q in its decomposition in a product of elementary reflectors. 
        ///     The size of <paramref name="Tau"/> must be at least max(1, <paramref name="k"/>).</param>
        /// <param name="C">Array containing the <paramref name="m"/>-by<paramref name="n"/> right hand side matrix. Its size is
        ///     max(1, <paramref name="ldC"/> * <paramref name="n"/>) for column major layout and 
        ///     max(1, <paramref name="ldC"/> * <paramref name="m"/>) for row major layout. Overwritten by the product Q * C, 
        ///     Q^T * C, C * Q, or C * Q^T (specified by <paramref name="side"/> and <paramref name="transQ"/>).</param>
        /// <param name="ldC">The leading dimension of <paramref name="C"/>. Constraints:  <paramref name="ldC"/> ≥
        ///     max(1, <paramref name="m"/>) for column major layout and max(1, <paramref name="n"/>) for row major layout.
        ///     </param>
        /// <returns></returns>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dormqr")]
        internal static extern int Dormqr(int matrixLayout, char side, char transQ,
            int m, int n, int k, double[] A, int ldA, double[] Tau, double[] C, int ldC);

        /// <summary>
        /// Inversion of positive definite matrix in full format, that underwent cholesky factorization by 
        /// <see cref="Lapack.Dpotrf(string, ref int, ref double, ref int, ref int)"/>. 
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dpotri")]
        internal static extern int Dpotri(int matrixLayout, char uplo, int n, double[] A, int ldA);

        /// <summary>
        /// Inversion of positive definite matrix in packed format, that underwent cholesky factorization by 
        /// <see cref="Lapack.Dpptrf(string, ref int, ref double, ref int)"/>. 
        /// </summary>
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "LAPACKE_dpptri")]
        internal static extern int Dpptri(int matrixLayout, char uplo, int n, double[] Ap);
    }
}
