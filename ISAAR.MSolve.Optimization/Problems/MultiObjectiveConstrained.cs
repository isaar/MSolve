using System;
using System.Collections.Generic;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Problems
{
    class MultiObjectiveConstrained: OptimizationProblem
    {
        private readonly Factory factory;

        public MultiObjectiveConstrained()
        {
            this.factory = new Factory();
            base.DesignFactory = factory;
        }

        public void AddObjectiveFunction(Func<double[], double> objectiveFunction)
        {
            factory.ObjectiveFunctions.Add(objectiveFunction);
        }

        public void AddConstraintFunction(Func<double[], double> constraintFunction)
        {
            factory.ConstraintFunctions.Add(constraintFunction);
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
            internal List<Func<double[], double>> ConstraintFunctions { get; }

            internal Factory()
            {
                this.ObjectiveFunctions = new List<Func<double[], double>>();
                this.ConstraintFunctions = new List<Func<double[], double>>();
            }

            public IDesign CreateDesign(double[] x)
            {
                double[] objectiveValues = new double[ObjectiveFunctions.Count];
                for (int i = 0; i < ObjectiveFunctions.Count; ++i)
                {
                    objectiveValues[i] = ObjectiveFunctions[i](x);
                }

                double[] constraintValues = new double[ConstraintFunctions.Count];
                for (int i = 0; i < ConstraintFunctions.Count; ++i)
                {
                    constraintValues[i] = ConstraintFunctions[i](x);
                }

                return new Design(objectiveValues, constraintValues);
            }
        }

        private class Design : IDesign
        {
            public double[] ObjectiveValues { get; }
            public double[] ConstraintValues { get; }

            internal Design(double[] objectiveValues, double[] constraintValues)
            {
                ObjectiveValues = objectiveValues;
                ConstraintValues = constraintValues;
            }
        }
    }
}
