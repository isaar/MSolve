using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging
{
    public class IncrementalDisplacementsLog
    {
        private readonly Dictionary<int, int[]> watchDofs;
        private readonly List<Dictionary<int, Dictionary<int, double>>> dofDisplacementsPerIter;

        /// <summary>
        /// Initializes a new instance of <see cref="IncrementalDisplacementsLog"/>.
        /// </summary>
        /// <param name="watchDofs">Which freedom degrees to track for each subdomain.</param>
        public IncrementalDisplacementsLog(Dictionary<int, int[]> watchDofs)
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

        public double GetTotalDisplacement(int iteration, int subdomainID, int dof) => dofDisplacementsPerIter[iteration][subdomainID][dof];
    }
}
