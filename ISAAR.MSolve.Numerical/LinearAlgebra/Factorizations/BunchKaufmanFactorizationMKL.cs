using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public class BunchKaufmanFactorizationMKL : IFactorizationMKL
    {
        double IFactorizationMKL.CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        DenseVector IFactorizationMKL.SolveLinearSystem(DenseVector rhs)
        {
            throw new NotImplementedException();
        }
    }
}
