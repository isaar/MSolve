using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    //TODO: Perhaps Addition, Subtraction and Scaling must be done without using delegates, for performance
    public interface IEntrywiseOperable: IIndexable2D
    {
        IEntrywiseOperable DoEntrywise(IEntrywiseOperable other, Func<double, double, double> binaryOperation);
        IEntrywiseOperable DoToAllEntries(Func<double, double> unaryOperation);
    }
}
