using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Optimization.Logging
{
    /// <summary>
    /// Logs the value of the objective function in each iteration and stores it for later use. If the optimization algorithm is
    /// population based, then the best objective value over all the members is logged.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ObjectiveFunctionLogger //: IOptimizationLogger
    {
        private readonly List<double> objectiveValues = new List<double>();

        public IReadOnlyList<double> ObjectiveFunctionValues => objectiveValues;

        public void Log(double objectiveValue) => objectiveValues.Add(objectiveValue);
        //public void Log(IOptimizationAlgorithm algorithm) => objectiveValues.Add(algorithm.BestFitness);
    }
}
