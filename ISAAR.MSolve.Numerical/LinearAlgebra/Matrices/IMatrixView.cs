using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public interface IMatrixView
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
        /// <param name="i">The row index: 0 &lt;= i &lt; <see cref="NumRows"/></param>
        /// <param name="j">The column index: 0 &lt;= j &lt; <see cref="NumColumns"/></param>
        /// <returns>The entry with indices i, j</returns>
        double this[int i, int j] { get; }

        /// <summary>
        /// Only structural non zeros
        /// </summary>
        /// <returns></returns>
        int NumNonZeros { get; }

        IMatrixView DoPointwise(IMatrixView other, Func<double, double, double> binaryOperation);
        IMatrixView DoToAllEntries(Func<double, double> unaryOperation);

        bool Equals(IMatrixView other, ValueComparer comparer = null);

        IMatrixView MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        IMatrixView MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        IVectorView MultiplyRight(IVectorView vector, bool transposeThis = false);

        //Perhaps these should return dense matrices directly
        IMatrixView Slice(int[] rowIndices, int[] colIndices);
        IMatrixView Slice(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive);

        IMatrixView Transpose();

        void WriteToConsole(Array2DFormatting format = null);
        void WriteToFile(string path, bool append = false, Array2DFormatting format = null);
    }
}
