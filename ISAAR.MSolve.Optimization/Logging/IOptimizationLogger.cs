using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Logging
{
    /// <summary>
    /// Common interface for all classes that log data during the optimization procedure.
    /// </summary>
    /// <remarks>
    /// Logging during the optimization procedure is performed using the Observer design pattern. Classes that implement
    /// <see cref="IOptimizationLogger"/> are the Observers, while the class that implements <see cref="IOptimizationAlgorithm"/>
    /// is the Subject. Data is transfered to the observers by pulling. <see cref="IOptimizationAlgorithm"/> provides properties 
    /// (or methods) that expose its current state. Once every iteration it calls <see cref="Log(IOptimizationAlgorithm)"/> on  
    /// its registered Observers and passes its instance. The Observers can then read only the state values they need by using 
    /// the provided properties.</remarks>
    public interface IOptimizationLogger
    {
        /// <summary>
        /// Instructs the <see cref="IOptimizationLogger"/> to log the current state of the <see cref="IOptimizationAlgorithm"/>.
        /// </summary>
        /// <param name="algorithm">The object that is observed (logged) by this <see cref="IOptimizationLogger"/></param>
        void Log(IOptimizationAlgorithm algorithm);
    }
}
