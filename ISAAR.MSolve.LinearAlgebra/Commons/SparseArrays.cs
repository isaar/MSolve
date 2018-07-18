//TODO: These should be delegated to C .dlls or to MKL if possible.
namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Low level array operations for sparse matrix formats.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class SparseArrays
    {
        /// <summary>
        /// Copied from 
        /// https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/csr.h#L376.
        /// Compute B = A for CSR matrix A, CSC matrix B.
        /// Also, with the appropriate arguments can also be used to:
        ///   - compute B = A ^ t for CSR matrix A, CSR matrix B
        ///   - compute B = A ^ t for CSC matrix A, CSC matrix B
        ///   - convert CSC->CSR
        /// Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
        /// </summary>
        /// <param name="n_row">Number of rows in A.</param>
        /// <param name="n_col">Number of columns in A.</param>
        /// <param name="Ap">Row pointers. Size = n_row+1.</param>
        /// <param name="Aj">Column indices. Size = nnz(A). They are not assumed to be in sorted order.</param>
        /// <param name="Ax">Non-zero values. Size = nnz(A).</param>
        /// <param name="Bp">Preallocated ouput argument. Column pointers. Size = n_col+1</param>
        /// <param name="Bi">Preallocated ouput argument. Row indices. Size = nnz(A). They will be in sorted order.</param>
        /// <param name="Bx">Preallocated ouput argument. Non-zero values. Size = nnz(A).</param>
        internal static void CsrToCsc(int n_row, int n_col, int[] Ap, int[] Aj, double[] Ax, int[] Bp, int[] Bi, double[] Bx)
        {
            int nnz = Ap[n_row];

            //compute number of non-zero entries per column of A 
            //std::fill(Bp, Bp + n_col, 0); //In C# the array is initilized to 0.0 by default

            for (int n = 0; n < nnz; n++)
            {
                Bp[Aj[n]]++;
            }

            //cumsum the nnz per column to get Bp[]
            for (int col = 0, cumsum = 0; col < n_col; col++)
            {
                int temp = Bp[col];
                Bp[col] = cumsum;
                cumsum += temp;
            }
            Bp[n_col] = nnz;

            for (int row = 0; row < n_row; row++)
            {
                for (int jj = Ap[row]; jj < Ap[row + 1]; jj++)
                {
                    int col = Aj[jj];
                    int dest = Bp[col];

                    Bi[dest] = row;
                    Bx[dest] = Ax[jj];

                    Bp[col]++;
                }
            }

            for (int col = 0, last = 0; col <= n_col; col++)
            {
                int temp = Bp[col];
                Bp[col] = last;
                last = temp;
            }
        }
    }
}
