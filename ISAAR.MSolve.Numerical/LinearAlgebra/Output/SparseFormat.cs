using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    // TODO: While this works for matrices that use arrays, it doesn't work for DOKs. Ideally the sparse format object should 
    // know how to print itself
    // TODO: I should also have a triangular format class for Symmetric, Lower and Upper classes.
    public class SparseFormat
    {
        public SparseFormat()
        {
            this.RawIndexArrays = new Dictionary<string, IReadOnlyList<int>>();
        }

        public Dictionary<string, IReadOnlyList<int>> RawIndexArrays { get; }
        public IReadOnlyList<double> RawValuesArray { get; set; }
        public string RawValuesTitle { get; set; } = "Values";
    }
}
