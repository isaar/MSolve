using ISAAR.MSolve.XFEM.Utilities;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

//TODO: make a sorted table
namespace ISAAR.MSolve.XFEM.Solvers
{
    class SolverLogger
    {
        private readonly string solverName;
        private readonly SortedDictionary<int, long> freedomDegrees;
        private readonly SortedDictionary<int, SortedDictionary<string, long>> durations;

        public SolverLogger(string solverName = "") //TODO: remove the default parameter value
        {
            this.solverName = solverName;
            this.freedomDegrees = new SortedDictionary<int, long>();
            this.durations = new SortedDictionary<int, SortedDictionary<string, long>>();
        }

        public void LogDofs(int iteration, long numDofs)
        {
            freedomDegrees.Add(iteration, numDofs);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="iteration">Iteration 0 = initialization, before any linear system solutions take place.</param>
        /// <param name="task"></param>
        /// <param name="duration"></param>
        public void LogDuration(int iteration, string task, long duration)
        {
            if (durations.ContainsKey(iteration)) durations[iteration].Add(task, duration);
            else
            {
                var newEntry = new SortedDictionary<string, long>();
                newEntry.Add(task, duration);
                durations.Add(iteration, newEntry);
            }
        }

        public void WriteToFile(string directory, bool append)
        {
            if (append && !File.Exists(directory)) append = false;
            using (var writer = new StreamWriter(directory, append))
            {
                writer.WriteLine("*********************************************************************");
                writer.WriteLine(solverName);
                writer.WriteLine("*********************************************************************");
                foreach (var wholeIteration in freedomDegrees)
                {
                    writer.Write($"Iteration {wholeIteration.Key}: ");
                    writer.WriteLine($"dofs = {wholeIteration.Value}");
                }
                foreach (var wholeIteration in durations)
                {
                    foreach (var taskTime in wholeIteration.Value)
                    {
                        writer.Write($"Iteration {wholeIteration.Key}: ");
                        writer.WriteLine($"Task = {taskTime.Key} - duration = {taskTime.Value} ms");
                    }
                }
                writer.WriteLine();
                writer.WriteLine();
            }
        }

    }
}
