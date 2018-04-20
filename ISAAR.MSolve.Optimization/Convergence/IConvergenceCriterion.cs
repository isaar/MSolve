using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Optimization.Convergence
{
    // Convergence criteria must be checked at the END of each iteration of the optimization process
    public interface IConvergenceCriterion
    {
        bool HasConverged(IOptimizationAlgorithm algorithm);
    }
}
