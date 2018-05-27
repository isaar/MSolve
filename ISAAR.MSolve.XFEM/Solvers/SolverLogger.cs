using ISAAR.MSolve.XFEM.Utilities;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

//TODO: make a sorted table
namespace ISAAR.MSolve.XFEM.Solvers
{
    class SolverLogger
    {
        private readonly SortedDictionary<int, long> freedomDegrees;
        private readonly SortedDictionary<int, SortedDictionary<string, long>> durations;

        public SolverLogger(string solverName = "") //TODO: remove the default parameter value
        {
            this.SolverName = solverName;
            this.freedomDegrees = new SortedDictionary<int, long>();
            this.durations = new SortedDictionary<int, SortedDictionary<string, long>>();
        }

        public string SolverName { get; }

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

        public void WriteToFile(string path, string title, bool append)
        {
            if (append && !File.Exists(path)) append = false;
            using (var writer = new StreamWriter(path, append))
            {
                WriteHeader(writer, title);

                // Dofs
                foreach (var wholeIteration in freedomDegrees)
                {
                    writer.Write($"Iteration {wholeIteration.Key}: ");
                    writer.WriteLine($"dofs = {wholeIteration.Value}");
                }

                // Durations
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

        public void WriteSumsToFile(string path, string title, bool append)
        {
            if (append && !File.Exists(path)) append = false;
            using (var writer = new StreamWriter(path, append))
            {
                WriteHeader(writer, title);

                // Dofs
                long minNumDofs = long.MaxValue;
                long maxNumDofs = long.MinValue;
                foreach (var numDofs in freedomDegrees.Values)
                {
                    if (numDofs < minNumDofs) minNumDofs = numDofs;
                    if (numDofs > maxNumDofs) maxNumDofs = numDofs;
                }
                writer.WriteLine($"Number of dofs: min = {minNumDofs}, max = {maxNumDofs}");

                // Durations
                var taskDurations = new SortedDictionary<string, long>();
                foreach (var wholeIteration in durations)
                {
                    int iteration = wholeIteration.Key;
                    foreach (var taskTime in wholeIteration.Value)
                    {
                        string task = taskTime.Key;
                        long duration = taskTime.Value;
                        taskDurations.TryGetValue(task, out long timeSoFar);
                        taskDurations[task] = timeSoFar + duration;
                    }
                }
                foreach (var taskTime in taskDurations)
                {
                    writer.WriteLine($"Task = {taskTime.Key}: total duration = {taskTime.Value} ms");
                }

                writer.WriteLine();
                writer.WriteLine();
            }
        }

        private void WriteHeader(StreamWriter writer, string title)
        {
            writer.WriteLine("*********************************************************************");
            writer.WriteLine(title);
            var culture = new CultureInfo("fr-FR");
            var date = DateTime.Now;
            writer.WriteLine("Date: " + date.ToString(culture));
            writer.WriteLine("Solver: " + SolverName);
            writer.WriteLine("*********************************************************************");
        }

    }
}
