using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Commons.Exceptions
{
    /// <summary>
    /// The exception that is thrown when the user selects parameters or strategies for implementing the "operators" of an 
    /// optimization algorithm that cannot work effectively together. As many of these strategies are stochastic, this exception 
    /// is not guarranteed to always be thrown during initialization or the first iteration of the optimization algorithm. It 
    /// might be thrown later or not even at all.
    /// </summary>
    class IncompatibleStrategiesException: Exception
    {
        public IncompatibleStrategiesException(): base()
        {
        }

        public IncompatibleStrategiesException(string message) : base(message)
        {
        }

        public IncompatibleStrategiesException(string message, Exception inner) : base(message, inner)
        {
        }
    }
}
