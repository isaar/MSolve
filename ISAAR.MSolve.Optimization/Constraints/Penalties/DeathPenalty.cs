using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Constraints.Penalties
{
    public class DeathPenalty : IPenaltyStatic
    {
        //private readonly IConstraintFunction[] inequalityConstraints;

        //public DeathPenalty(IConstraintFunction[] inequalityConstraints)
        //{
        //    this.inequalityConstraints = inequalityConstraints;
        //}

        public double Evaluate(double fitness, IDesign design)
        {
            foreach (var constraintValue in design.ConstraintValues)
            {
                if (constraintValue > 0) return double.MaxValue;
            }
            return fitness;
        }


    }
}
