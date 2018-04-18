using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    public class BunchKaufmanFactorization : IFactorization
    {
        double IFactorization.CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        Vector IFactorization.SolveLinearSystem(Vector rhs)
        {
            throw new NotImplementedException();
        }
    }
}
