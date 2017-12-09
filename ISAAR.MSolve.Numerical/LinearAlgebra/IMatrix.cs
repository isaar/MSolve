using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public interface IMatrix: IMatrixView
    {
        new double this[int row, int col] { get; set; } // Concrete classes only need to implement the getter once.

        void DoPointwiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation);
        void DoToAllEntriesIntoThis(Func<double, double> unaryOperation);
        void SetAll(double value); // Awkward with DoToAllEntriesIntoThis
    }
}
