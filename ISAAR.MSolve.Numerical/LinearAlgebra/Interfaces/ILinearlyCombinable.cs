using System.Collections.Generic;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces
{
    public interface ILinearlyCombinable<T> where T : IMatrix2D
    {
        void LinearCombination(IList<double> coefficients, IList<T> matrices);
    }
}
