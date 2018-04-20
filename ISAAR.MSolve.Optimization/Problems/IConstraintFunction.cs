using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Problems
{

    public interface IConstraintFunction
    {
        double Evaluate(double[] x);
    }
}
