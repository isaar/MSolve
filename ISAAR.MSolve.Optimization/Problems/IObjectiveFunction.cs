using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Problems
{
    /// <summary>
    /// Represents an objective function of an optimization problem. It will be called be the selected optimization algorithm. 
    /// The user hooks the model he wants to optimize by implementing this interface. 
    /// All optimizers will try to find the global minimum of the IObjectiveFunction. 
    /// To maximize a function the user needs to provide its negative -f(x).
    /// </summary>
    public interface IObjectiveFunction
    {

        /// <summary>
        /// Calculates the value of the IObjectiveFunction object. 
        /// This call may be time consuming (e.g. launching a FEM simulation). 
        /// Unless it is executed in a different thread it may block the current iteration of the optimization algorithm. 
        /// </summary>
        /// <param name="x">The continuous decision variable vector.</param>
        /// <returns>The value of the IObjectiveFunction object.</returns>
        double Evaluate(double[] x);
    }
}
