using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public interface ISliceable1D
    {
        VectorMKL Slice(int[] indices);
        VectorMKL Slice(int startInclusive, int endExclusive);
    }
}
