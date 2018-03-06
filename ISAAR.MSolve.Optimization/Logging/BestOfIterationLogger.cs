using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

    namespace ISAAR.MSolve.Optimization.Logging
{
    /// <summary>
    /// An <see cref="IOptimizationLogger"/> that logs the best values of the objective function and the corresponding design 
    /// variables found at each iteration of the <see cref="IOptimizationAlgorithm"/>.
    /// </summary>
    public class BestOfIterationLogger: IOptimizationLogger
    {
        private List<double[]> bestContinuousVariables;
        private List<int[]> bestIntegerVariables;
        private List<double> bestObjectives;

        /// <summary>
        /// Creates a new instance of <see cref="BestOfIterationLogger"/>.
        /// </summary>
        public BestOfIterationLogger()
        {
            bestContinuousVariables = new List<double[]>();
            bestIntegerVariables = new List<int[]>();
            bestObjectives = new List<double>();
        }

        public void Log(IOptimizationAlgorithm algorithm)
        {
            bestContinuousVariables.Add(algorithm.BestPosition);
            //bestIntegerVariables.Add();
            bestObjectives.Add(algorithm.BestFitness);
        }

        /// <summary>
        /// Writes the logged values of the objective function and the corresponding design variables for each iteration 
        /// to the console.
        /// </summary>
        public void PrintToConsole()
        {
            Console.Write("Initialization: ");
            PrintEntryToConsole(0);
            for (int iteration = 0; iteration < bestObjectives.Count-1; ++iteration)
            {
                Console.Write("Iteration " + iteration + ": ");
                PrintEntryToConsole(iteration + 1);
            }
        }

        /// <summary>
        /// Writes the logged values of the objective function and the corresponding design variables of the iteration denoted by
        /// index to the console.
        /// </summary>
        /// <param name="index">The index assigned to each iteration. Initialization is assigned index = 0.</param>
        private void PrintEntryToConsole(int index)
        {
            Console.Write("continuous variables = " + ArrayToString<double>(bestContinuousVariables[index]));
            Console.Write(" , integer variables = " + ArrayToString<int>(bestIntegerVariables[index]));
            Console.WriteLine(" , objective value = " + bestObjectives[index]);
        }

        private static string ArrayToString<T>(T[] array)
        {
            StringBuilder builder = new StringBuilder();
            builder.Append("{ ");
            foreach (var entry in array)
            {
                builder.Append(entry);
                builder.Append(' ');
            }
            builder.Append("}");
            return builder.ToString();
        }
        
    }
}
