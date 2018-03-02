using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public class BunchKaufmanFactorizationMKL : IFactorizationMKL
    {
        double IFactorizationMKL.CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        VectorMKL IFactorizationMKL.SolveLinearSystem(VectorMKL rhs)
        {
            throw new NotImplementedException();
        }
    }
}
