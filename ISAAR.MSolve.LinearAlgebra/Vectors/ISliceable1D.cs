using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public interface ISliceable1D
    {
        Vector Slice(int[] indices);
        Vector Slice(int startInclusive, int endExclusive);
    }
}
