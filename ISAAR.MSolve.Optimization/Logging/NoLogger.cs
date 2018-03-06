using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Optimization.Logging
{
    /// <summary>
    /// An <see cref="IOptimizationLogger"/> that does nothing.
    /// </summary>
    public class NoLogger : IOptimizationLogger
    {
        public void Log(IOptimizationAlgorithm algorithm)
        {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void PrintToConsole()
        {
        }
    }
}
