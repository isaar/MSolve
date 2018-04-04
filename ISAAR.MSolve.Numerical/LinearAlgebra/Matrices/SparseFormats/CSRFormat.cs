using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.SparseFormats
{
    public class CSRFormat: ISparseFormat
    {
        public CSRFormat(double[] values, int[] columnIndices, int[] rowOffsets)
        {
            RawValuesArray = values;
            RawIndexArrays = new Dictionary<string, IReadOnlyList<int>>(2)
            {
                { "Column indices", columnIndices },
                { "Row offsets", rowOffsets }
            };
        }

        public IReadOnlyList<double> RawValuesArray { get; }
        public Dictionary<string, IReadOnlyList<int>> RawIndexArrays { get; }
    }
}
