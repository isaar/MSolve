using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public interface IVector: IVectorView
    {
        new double this[int index] { get; set; } // Concrete classes only need to implement the getter once.
        
        void AddIntoThis(IVectorView vector);
        void DoPointwiseIntoThis(IVectorView vector, Func<double, double, double> operation);
        void DoToAllEntriesIntoThis(Func<double, double> operation);
        void MultiplyPointwiseIntoThis(IVectorView vector);
        void MultiplyScalarIntoThis(double scalar);
        void SetAll(double value);
        void SubtractIntoThis(IVectorView vector);
    }
}
