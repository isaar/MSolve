using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

//TODO: Use enums instead of strings for the solver task and dof category. Or use interfaces & enum classes, to adhere to 
//      open-closed principle.
namespace ISAAR.MSolve.Solvers
{
    public class SolverLogger
    {
        private readonly string solverName;
        private readonly List<(int iterations, double residualNormRatio)> iterativeAlgorithmData = new List<(int, double)>();
        private readonly List<SortedDictionary<string, long>> taskDurations = new List<SortedDictionary<string, long>>();
        private readonly List<SortedDictionary<string, int>> numDofsPerCategory = new List<SortedDictionary<string, int>>();
        private int currentIteration;

        public SolverLogger(string solverName)
        {
            this.solverName = solverName;
            currentIteration = 0;
            taskDurations.Add(new SortedDictionary<string, long>());
            numDofsPerCategory.Add(new SortedDictionary<string, int>());
        }

        public void LogTaskDuration(string task, long duration)
            => taskDurations[currentIteration][task] = duration;

        public void LogNumDofs(string category, int numDofs)
            => numDofsPerCategory[currentIteration][category] = numDofs;

        public void LogIterativeAlgorithm(int iterations, double residualNormRatio)
            => iterativeAlgorithmData[currentIteration] = (iterations, residualNormRatio);

        /// <summary>
        /// Each iteration is defined by the solution phase of ISolver. Dof ordering and matrix assembly may also be included, 
        /// but they are not necessarily repeated in all analyses. Thus call it at the end of the Solve() method.
        /// </summary>
        public void IncrementIteration() //TODO: Forgetting to call this is easy. A better design is needed. 
        {
            ++currentIteration;
            taskDurations.Add(new SortedDictionary<string, long>());
            numDofsPerCategory.Add(new SortedDictionary<string, int>());
        }

        public void WriteToFile(string path, string title, bool append)
        {
            if (append && !File.Exists(path)) append = false;
            using (var writer = new StreamWriter(path, append))
            {
                WriteHeader(writer, title);

                // Dofs
                for (int i = 0; i <= currentIteration; ++i)
                {
                    foreach (var categoryDofs in numDofsPerCategory[i])
                    {
                        writer.Write($"Analysis iteration {i}: ");
                        writer.WriteLine($"Dof category = {categoryDofs.Key} - num dofs = {categoryDofs.Value}");
                    }
                }

                // Durations
                for (int i = 0; i <= currentIteration; ++i)
                {
                    foreach (var taskDuration in taskDurations[i])
                    {
                        writer.Write($"Analysis iteration {i}: ");
                        writer.WriteLine($"Task = {taskDuration.Key} - duration = {taskDuration.Value} ms");
                    }
                }

                // Iterative algorithm data
                for (int i = 0; i <= currentIteration; ++i)
                {
                    (int iter, double res) = iterativeAlgorithmData[i];
                    writer.Write($"Analysis iteration {i}: ");
                    writer.WriteLine($"Iterative algorithm: num iterations = {iter} - redisual norm ratio = {res} ms");
                }

                writer.WriteLine();
                writer.WriteLine();
            }
        }

        public void WriteAggregatesToFile(string path, string title, bool append)
        {
            if (append && !File.Exists(path)) append = false;
            using (var writer = new StreamWriter(path, append))
            {
                WriteHeader(writer, title);

                // Dofs
                var numDofsRange = new SortedDictionary<string, (int min, int max)>();
                for (int i = 0; i <= currentIteration; ++i)
                {
                    foreach (var categoryDofs in numDofsPerCategory[i])
                    {
                        string category = categoryDofs.Key;
                        int numDofs = categoryDofs.Value;

                        bool exists = numDofsRange.TryGetValue(category, out (int min, int max) rangeSoFar); 
                        if (!exists)
                        {
                            rangeSoFar.min = int.MaxValue;
                            rangeSoFar.max = int.MinValue;
                        }
                        numDofsRange[category] =
                            (
                                (numDofs < rangeSoFar.min) ? numDofs : rangeSoFar.min,
                                (numDofs > rangeSoFar.max) ? numDofs : rangeSoFar.max
                            );
                    }
                }
                foreach (var categoryRange in numDofsRange)
                {
                    string category = categoryRange.Key;
                    (int min, int max) = categoryRange.Value;
                    writer.WriteLine($"Dof category = {category}: min = {min} - max = {max}");
                }

                // Durations
                long totalDuration = 0;
                var taskTotalDurations = new SortedDictionary<string, long>();
                for (int i = 0; i <= currentIteration; ++i)
                {
                    foreach (var taskDuration in taskDurations[i])
                    {
                        string task = taskDuration.Key;
                        long duration = taskDuration.Value;
                        totalDuration += duration;
                        taskTotalDurations.TryGetValue(task, out long timeSoFar);
                        taskTotalDurations[task] = timeSoFar + duration;
                    }
                }
                foreach (var taskTime in taskTotalDurations)
                {
                    writer.WriteLine($"Task = {taskTime.Key}: total duration = {taskTime.Value} ms");
                }
                writer.WriteLine($"All tasks: total duration = {totalDuration} ms");

                // Iterative algorithm data
                int minIterations = int.MaxValue, maxIterations = int.MinValue;
                double minResNorm = double.MaxValue, maxResNorm = double.MinValue;
                for (int i = 0; i <= currentIteration; ++i)
                {
                    (int iter, double res) = iterativeAlgorithmData[i];
                    minIterations = (iter < minIterations) ? iter : minIterations;
                    maxIterations = (iter > maxIterations) ? iter : maxIterations;
                    minResNorm = (res < minResNorm) ? res : minResNorm;
                    maxResNorm = (res < maxResNorm) ? res : maxResNorm;
                }
                Console.WriteLine($"Iterative algorithm iterations: min = {minIterations} - max = {maxIterations}");
                Console.WriteLine($"Iterative algorithm residual norm ratio: min = {minResNorm} - max = {maxResNorm}");

                writer.WriteLine();
                writer.WriteLine();
            }
        }

        private void WriteHeader(StreamWriter writer, string title)
        {
            writer.WriteLine();
            writer.WriteLine("*********************************************************************");
            writer.WriteLine(title);
            var culture = new CultureInfo("fr-FR");
            var date = DateTime.Now;
            writer.WriteLine("Date: " + date.ToString(culture));
            writer.WriteLine("Solver: " + solverName);
            writer.WriteLine("*********************************************************************");
        }
    }
}
