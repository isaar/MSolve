using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.SparseFormats
{
    public class CSCFormat : ISparseFormat
    {
        public CSCFormat(double[] values, int[] rowIndices, int[] columnOffsets)
        {
            RawValuesArray = values;
            RawIndexArrays = new Dictionary<string, IReadOnlyList<int>>(2)
            {
                { "Row indices", rowIndices },
                { "Column offsets", columnOffsets }
            };
        }

        public IReadOnlyList<double> RawValuesArray { get; }
        public Dictionary<string, IReadOnlyList<int>> RawIndexArrays { get; }
    }
}
