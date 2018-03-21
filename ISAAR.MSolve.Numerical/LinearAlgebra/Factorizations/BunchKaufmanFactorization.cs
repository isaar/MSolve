using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public class BunchKaufmanFactorization : IFactorization
    {
        double IFactorization.CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        VectorMKL IFactorization.SolveLinearSystem(VectorMKL rhs)
        {
            throw new NotImplementedException();
        }
    }
}
