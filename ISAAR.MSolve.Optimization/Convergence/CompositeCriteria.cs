using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Convergence
{
    public static class CompositeCriteria
    {
        #region static methods
        public static IConvergenceCriterion NOT(IConvergenceCriterion originalCriterion)
        {
            return new NotCriterion(originalCriterion);
        }

        public static IConvergenceCriterion AND(IConvergenceCriterion criterion1, IConvergenceCriterion criterion2)
        {
            return new AndCriterion(criterion1, criterion2);
        }

        public static IConvergenceCriterion OR(IConvergenceCriterion criterion1, IConvergenceCriterion criterion2)
        {
            return new OrCriterion(criterion1, criterion2);
        }
        #endregion

        #region private classes
        private class NotCriterion : IConvergenceCriterion
        {
            private readonly IConvergenceCriterion originalCriterion;

            internal NotCriterion(IConvergenceCriterion originalCriterion)
            {
                this.originalCriterion = originalCriterion;
            }

            public bool HasConverged(IOptimizationAlgorithm algorithm)
            {
                return !originalCriterion.HasConverged(algorithm);
            }
        }

        private class AndCriterion : IConvergenceCriterion
        {
            private readonly IConvergenceCriterion criterion1;
            private readonly IConvergenceCriterion criterion2;

            internal AndCriterion(IConvergenceCriterion criterion1, IConvergenceCriterion criterion2)
            {
                this.criterion1 = criterion1;
                this.criterion2 = criterion2;
            }

            public bool HasConverged(IOptimizationAlgorithm algorithm)
            {
                return criterion1.HasConverged(algorithm) && criterion2.HasConverged(algorithm);
            }
        }

        private class OrCriterion : IConvergenceCriterion 
        {
            private readonly IConvergenceCriterion criterion1;
            private readonly IConvergenceCriterion criterion2;

            internal OrCriterion(IConvergenceCriterion criterion1, IConvergenceCriterion criterion2)
            {
                this.criterion1 = criterion1;
                this.criterion2 = criterion2;
            }

            public bool HasConverged(IOptimizationAlgorithm algorithm)
            {
                return criterion1.HasConverged(algorithm) || criterion2.HasConverged(algorithm);
            }
        }
        #endregion
    }


}
