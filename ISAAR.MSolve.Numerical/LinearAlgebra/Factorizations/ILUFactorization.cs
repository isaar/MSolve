using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public interface ILUFactorization
    {
        double CalcDeterminant();
        IVector Solve(IVectorView rhs);
    }
}
