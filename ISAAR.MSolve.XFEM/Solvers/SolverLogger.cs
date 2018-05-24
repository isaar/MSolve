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
        private readonly SortedDictionary<int, SortedDictionary<string, long>> durations;

        public SolverLogger()
        {
            this.durations = new SortedDictionary<int, SortedDictionary<string, long>>();
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

        public void WriteToFile(string path, string header, bool append)
        {
            using (var writer = new StreamWriter(path, append))
            {
                writer.WriteLine("*********************************************************************");
                writer.WriteLine(header);
                writer.WriteLine("*********************************************************************");
                writer.WriteLine();

                foreach (var wholeIteration in durations)
                {
                    writer.Write($"Iteration {wholeIteration.Key}: ");
                    foreach (var taskTime in wholeIteration.Value)
                    {
                        writer.WriteLine($"Task = {taskTime.Key} - duration = {taskTime.Value} ms");
                    }
                }
                writer.WriteLine();
                writer.WriteLine();
            }
        }

    }
}
