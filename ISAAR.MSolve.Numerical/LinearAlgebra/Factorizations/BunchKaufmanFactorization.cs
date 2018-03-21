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

        Vector IFactorization.SolveLinearSystem(Vectors.VectorMKL rhs)
        {
            throw new NotImplementedException();
        }
    }
}
