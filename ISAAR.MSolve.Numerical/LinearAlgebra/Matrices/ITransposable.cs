using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public interface ITransposable: IIndexable2D
    {
        ITransposable Transpose();
    }
}
