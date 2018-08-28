//WARNING: IGeneralMatrixBuilder.AddSubmatrix() appears to be the same as ISymmetricMatrixBuilder.AddSubmatrix(). However,  
// their semantics are different. ISymmetricMatrixBuilder.AddSubmatrix() cannot be used for adding symmetric matrices, as it 
// will add the off-diagonal entries twice. This is due to the implementations of ISymmetricMatrixBuilder. On the other hand 
// the implementations of IGeneralMatrixBuilder make it safe to pass in symmetric matrices to AddSubmatrix(), although  
// there are more efficient methods.
namespace ISAAR.MSolve.LinearAlgebra.Matrices.Builders
{
    /// <summary>
    /// Facilitates the construction of large sparse matrices, usually by adding smaller submatrices.  
    /// The large matrices and their properties will be characterized as "global" in this namespace.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IMatrixBuilder: IIndexable2D
    {
        /// <summary>
        /// Sets the entry (<paramref name="rowIdx"/>, <paramref name="colIdx"/>) to the provided value.
        /// </summary>
        /// <param name="rowIdx">The row index of the entry to set. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="colIdx">The column index of the entry to set. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="IIndexable2D.NumColumns"/>.</param>
        /// <returns></returns>
        new double this[int rowIdx, int colIdx] { set; }

        /// <summary>
        /// Adds the provided <paramref name="value"/> to the entry (<paramref name="rowIdx"/>, <paramref name="colIdx"/>). 
        /// </summary>
        /// <param name="rowIdx">The row index of the entry to modify. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="colIdx">The column index of the entry to modify. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; this.<see cref="IIndexable2D.NumColumns"/>.</param>
        /// <param name="value">The value that will be added to the entry (<paramref name="colIdx"/>, <paramref name="colIdx"/>).
        ///     </param>
        void AddToEntry(int rowIdx, int colIdx, double value);
    }
}
