using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Problems
{
    //public delegate double RealFunction(double[] x);

    public class SingleObjectiveUnconstrained: OptimizationProblem
    {
        /// <summary>
        /// Single objective unconstrained mathematical benchmark problems need only to set the objective function 
        /// (and the bounds)
        /// </summary>
        public Func<double[], double> ObjectiveFunction
        {
            set { base.DesignFactory = new Factory(value); }
        }

        /// <summary>
        /// Do not call this property. Its only purpose is to prevent subclasses of <see cref="SingleObjectiveUnconstrained"/>
        /// from setting the <see cref="OptimizationProblem.DesignFactory"/> property directly, although it is still 
        /// possible.
        /// </summary>
        public new IDesignFactory DesignFactory
        {
            get { return base.DesignFactory; }
        }

        private class Factory: IDesignFactory
        {
            private readonly Func<double[], double> objectiveFunction;

            internal Factory(Func<double[], double> objectiveFunction)
            {
                this.objectiveFunction = objectiveFunction;
            }

            public IDesign CreateDesign(double[] x)
            {
                return new Design(objectiveFunction(x));
            }
        }

        private class Design : IDesign
        {
            public double[] ObjectiveValues { get; }

            public double[] ConstraintValues
            {
                get { throw new InvalidOperationException("There are no constraints"); }
            }

            internal Design(double objectiveValue)
            {
                ObjectiveValues = new double[] { objectiveValue };
            }
        }
    }
}
