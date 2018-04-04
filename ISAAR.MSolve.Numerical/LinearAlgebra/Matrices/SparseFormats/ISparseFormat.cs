using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.SparseFormats
{
    public interface ISparseFormat
    {
        IReadOnlyList<double> RawValuesArray { get; }
        Dictionary<string, IReadOnlyList<int>> RawIndexArrays { get; }
    }
}
