using ISAAR.MSolve.XFEM.Utilities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: Separate dof enumeration, global matrix assembly and linear system solution times
//TODO: use something more sophisticated than Stopwatch.ElapsedMilliseconds
//TODO: encapsulate the time measurement in methods of this class, rather than letting each solver do it.
namespace ISAAR.MSolve.XFEM.Solvers
{
    class SolverLoggerOLD
    {
        private Table<string, int, long> durations;

        public SolverLoggerOLD()
        {
            this.SolutionTimes = new List<long>();
        }

        public void LogDuration(string task, int iteration, long duration)
        {
            durations[task, iteration] = duration;
        }

        /// <summary>
        /// The duration in milliseconds of the initialization.
        /// </summary>
        public long InitializationTime { get; set; }

        /// <summary>
        /// The durations in milliseconds of all linear system solutions, if the solver is called more than once.
        /// </summary>
        public List<long> SolutionTimes { get; }

        /// <summary>
        /// The total duration of initialization and all linear system solutions
        /// </summary>
        /// <returns></returns>
        public long CalcTotalTime()
        {
            // TODO: perhaps I should use a data type that will not overflow
            long sum = InitializationTime;
            for (int i = 0; i < SolutionTimes.Count; ++i) sum += SolutionTimes[i];
            return sum;
        }

        public void Clear()
        {
            InitializationTime = 0;
            SolutionTimes.Clear();
        }
    }
}
