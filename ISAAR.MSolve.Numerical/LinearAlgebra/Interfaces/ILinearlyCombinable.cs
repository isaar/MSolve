using System.Collections.Generic;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces
{
    public interface ILinearlyCombinable<T> where T : IMatrix2D
    {
        void LinearCombination(IList<double> coefficients, IList<T> matrices);
    }

    public interface ILinearlyCombinable
    {
        void LinearCombination(IList<double> coefficients, IList<IMatrix2D> matrices);
    }
}
