using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging
{
    public class TotalDisplacementsPerIterationLog: IAnalyzerLog
    {
        private readonly Dictionary<int, int[]> watchDofs;
        private readonly List<Dictionary<int, Dictionary<int, double>>> dofDisplacementsPerIter;

        /// <summary>
        /// Initializes a new instance of <see cref="TotalDisplacementsPerIterationLog"/>.
        /// </summary>
        /// <param name="watchDofs">Which freedom degrees to track for each subdomain.</param>
        public TotalDisplacementsPerIterationLog(Dictionary<int, int[]> watchDofs)
        {
            this.watchDofs = watchDofs;
            this.dofDisplacementsPerIter = new List<Dictionary<int, Dictionary<int, double>>>();
        }

        /// <summary>
        /// Stores the total displacements = u_converged + du, for a new iteration.
        /// </summary>
        /// <param name="totalDisplacements">The total displacements for each subdomain.</param>
        public void StoreDisplacements(Dictionary<int, Vector> totalDisplacements)
        {
            var currentIterDisplacements = new Dictionary<int, Dictionary<int, double>>();
            foreach (var subdomainDofsPair in watchDofs)
            {
                int subdomainID = subdomainDofsPair.Key;
                var subdomainDisplacements = new Dictionary<int, double>();
                foreach (int dof in subdomainDofsPair.Value)
                {
                    subdomainDisplacements[dof] = totalDisplacements[subdomainID][dof];
                }
                currentIterDisplacements[subdomainID] = subdomainDisplacements;
            }
            dofDisplacementsPerIter.Add(currentIterDisplacements);
        }

        /// <summary>
        /// Stores the total displacements = u_converged + du, for a new iteration.
        /// </summary>
        /// <param name="totalDisplacements">The total displacements for each subdomain.</param>
        public void StoreDisplacements_v2(Dictionary<int, LinearAlgebra.Vectors.IVector> totalDisplacements)
        {
            var currentIterDisplacements = new Dictionary<int, Dictionary<int, double>>();
            foreach (var subdomainDofsPair in watchDofs)
            {
                int subdomainID = subdomainDofsPair.Key;
                var subdomainDisplacements = new Dictionary<int, double>();
                foreach (int dof in subdomainDofsPair.Value)
                {
                    subdomainDisplacements[dof] = totalDisplacements[subdomainID][dof];
                }
                currentIterDisplacements[subdomainID] = subdomainDisplacements;
            }
            dofDisplacementsPerIter.Add(currentIterDisplacements);
        }

        public double GetTotalDisplacement(int iteration, int subdomainID, int dof) => dofDisplacementsPerIter[iteration][subdomainID][dof];

        public void StoreResults(DateTime startTime, DateTime endTime, IVector solution)
        {
            throw new NotImplementedException();
        }
    }
}
