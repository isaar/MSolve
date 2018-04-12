using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.SuiteSparse
{
    internal static class SuiteSparseUtilities
    {
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_create_common")]
        internal static extern IntPtr CreateCommon();

        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_destroy_common")]
        internal static extern void DestroyCommon(ref IntPtr common);

        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_destroy_factor")]
        internal static extern void DestroyFactor(ref IntPtr factorizedMatrix, IntPtr common);

        /// <summary>
        /// Factorize a symmetric matrix using cholesky algorithm. The matrix is in csc form, with only the upper triangle stored.
        /// If cholesky is successful: -1 is returned and factorizedMatrix points to the factorized upper triangle.
        /// If the matrix is not positive definite: the index(0-based) of the column where cholesky failed
        /// and factorizedMatrix = <see cref="IntPtr.Zero"/>.
        /// If the something another failure occurs, such as memory not being sufficient due to excessive fill-in: -2 is returned
        /// and factorizedMatrix = <see cref="IntPtr.Zero"/>.
        /// </summary>
        /// <param name="order"></param>
        /// <param name="nnz"></param>
        /// <param name="values"></param>
        /// <param name="rowIndices"></param>
        /// <param name="colOffsets"></param>
        /// <param name="factorizedMatrix"></param>
        /// <param name="common"></param>
        /// <returns></returns>
        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_factorize_cscupper")]
        internal static extern int FactorizeCSCUpper(int order, int nnz, double[] values, int[] rowIndices, int[] colOffsets, 
            out IntPtr factorizedMatrix, IntPtr common);

        [DllImport("suitesparse_utilities.dll", EntryPoint = "util_solve")]
        internal static extern void Solve(int order, IntPtr factorizedMatrix, double[] rhs, double[] outSolution, IntPtr common);
    }
}
