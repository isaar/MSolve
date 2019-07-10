using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.OC.Convergence
{
    public static class CompositeCriteria
    {
        #region static methods
        public static IOptimalityCriteriaConvergence NOT(IOptimalityCriteriaConvergence originalCriterion)
            => new NotCriterion(originalCriterion);

        public static IOptimalityCriteriaConvergence AND(IOptimalityCriteriaConvergence criterion1, 
            IOptimalityCriteriaConvergence criterion2)
            => new AndCriterion(criterion1, criterion2);

        public static IOptimalityCriteriaConvergence OR(IOptimalityCriteriaConvergence criterion1, 
            IOptimalityCriteriaConvergence criterion2)
            => new OrCriterion(criterion1, criterion2);
        #endregion

        #region private classes
        private class NotCriterion : IOptimalityCriteriaConvergence
        {
            private readonly IOptimalityCriteriaConvergence originalCriterion;

            internal NotCriterion(IOptimalityCriteriaConvergence originalCriterion)
            {
                this.originalCriterion = originalCriterion;
            }

            public bool HasConverged(int currentIteration, double currentObjectiveFunction, IVectorView nextDesignVariables)
                => !originalCriterion.HasConverged(currentIteration, currentObjectiveFunction, nextDesignVariables);
        }

        private class AndCriterion : IOptimalityCriteriaConvergence
        {
            private readonly IOptimalityCriteriaConvergence criterion1;
            private readonly IOptimalityCriteriaConvergence criterion2;

            internal AndCriterion(IOptimalityCriteriaConvergence criterion1, IOptimalityCriteriaConvergence criterion2)
            {
                this.criterion1 = criterion1;
                this.criterion2 = criterion2;
            }

            public bool HasConverged(int currentIteration, double currentObjectiveFunction, IVectorView nextDesignVariables)
                => criterion1.HasConverged(currentIteration, currentObjectiveFunction, nextDesignVariables) 
                    && criterion2.HasConverged(currentIteration, currentObjectiveFunction, nextDesignVariables);
        }

        private class OrCriterion : IOptimalityCriteriaConvergence
        {
            private readonly IOptimalityCriteriaConvergence criterion1;
            private readonly IOptimalityCriteriaConvergence criterion2;

            internal OrCriterion(IOptimalityCriteriaConvergence criterion1, IOptimalityCriteriaConvergence criterion2)
            {
                this.criterion1 = criterion1;
                this.criterion2 = criterion2;
            }

            public bool HasConverged(int currentIteration, double currentObjectiveFunction, IVectorView nextDesignVariables)
                => criterion1.HasConverged(currentIteration, currentObjectiveFunction, nextDesignVariables)
                    || criterion2.HasConverged(currentIteration, currentObjectiveFunction, nextDesignVariables);
        }
        #endregion
    }


}
