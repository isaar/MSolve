using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public interface IIndexable2D
    {
        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        int NumRows { get; }

        /// <summary>
        /// The entry with row index = i and column index = j. 
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= rowIdx &lt; <see cref="NumRows"/></param>
        /// <param name="colIdx">The column index: 0 &lt;= colIdx &lt; <see cref="NumColumns"/></param>
        /// <returns>The entry with indices i, j</returns>
        double this[int rowIdx, int colIdx] { get; }

        bool Equals(IIndexable2D other, double tolerance = 1e-13);
    }
}
