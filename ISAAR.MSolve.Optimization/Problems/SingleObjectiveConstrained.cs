using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Problems;

namespace ISAAR.MSolve.Optimization.Problems
{
    public class SingleObjectiveConstrained : OptimizationProblem
    {
        private readonly Factory factory;

        public SingleObjectiveConstrained()
        {
            this.factory = new Factory();
            base.DesignFactory = factory;
        }

        public Func<double[], double> ObjectiveFunction
        {
            set { factory.ObjectiveFunction = value; }
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
            internal Func<double[], double> ObjectiveFunction { get; set; }
            internal List<Func<double[], double>> ConstraintFunctions { get; }

            internal Factory()
            {
                this.ConstraintFunctions = new List<Func<double[], double>>();
            }

            public IDesign CreateDesign(double[] x)
            {
                double objectiveValue = ObjectiveFunction(x);
                double[] constraintValues = new double[ConstraintFunctions.Count];
                for (int i = 0; i < ConstraintFunctions.Count; ++i)
                {
                    constraintValues[i] = ConstraintFunctions[i](x);
                }
                return new SingleObjectiveConstrainedDesign(objectiveValue, constraintValues);
            }
        }

        private class SingleObjectiveConstrainedDesign : IDesign
        {
            public double[] ObjectiveValues { get; }
            public double[] ConstraintValues { get; }

            internal SingleObjectiveConstrainedDesign(double objectiveValue, double[] constraintValues)
            {
                ObjectiveValues = new double[] { objectiveValue };
                ConstraintValues = constraintValues;
            }
        }
    }
}