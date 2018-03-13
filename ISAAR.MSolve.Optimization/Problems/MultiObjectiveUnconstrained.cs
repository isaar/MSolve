using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Problems
{
    class MultiObjectiveUnconstrained: OptimizationProblem
    {
        private readonly Factory factory;

        public MultiObjectiveUnconstrained()
        {
            this.factory = new Factory();
            base.DesignFactory = factory;
        }

        public void AddObjectiveFunction(Func<double[], double> objectiveFunction)
        {
            factory.ObjectiveFunctions.Add(objectiveFunction);
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

        private class Factory : IDesignFactory
        {
            internal List<Func<double[], double>> ObjectiveFunctions { get; }

            internal Factory()
            {
                this.ObjectiveFunctions = new List<Func<double[], double>>();
            }

            public IDesign CreateDesign(double[] x)
            {
                double[] objectiveValues = new double[ObjectiveFunctions.Count];
                for (int i = 0; i < ObjectiveFunctions.Count; ++i)
                {
                    objectiveValues[i] = ObjectiveFunctions[i](x);
                }
                return new Design(objectiveValues);
            }
        }

        private class Design : IDesign
        {
            public double[] ObjectiveValues { get; }

            public double[] ConstraintValues
            {
                get { throw new InvalidOperationException("There are no constraints"); }
            }

            internal Design(double[] objectiveValues)
            {
                ObjectiveValues = objectiveValues;
            }
        }
    }
}
