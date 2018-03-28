using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public interface IMatrixView: IIndexable2D
    {
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
        void WriteToFile(string path, bool append = false, Array2DFormatting format = null); //Perhaps the fromater must be a static field, since it depends on each matrix type
    }
}
