using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Optimization.Problems
{
    //TODO: use a design class for this.
    public delegate double EqualityConstraint(Vector x);

    public interface IConstraintFunction
    {
        double Evaluate(double[] x);
    }
}
