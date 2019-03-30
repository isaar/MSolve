using System;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;

//TODO: we need direct access to the observers themselves. Perhaps Builders to assist the user, but the analyzer should be 
//      passed the logger directly.
namespace ISAAR.MSolve.Logging
{
    public class TotalDisplacementsLog : IAnalyzerLog_v2
    {
        private readonly ISubdomain_v2 subdomain;
        private readonly DofTable watchedDofs;
        private double[] displacements;

        public TotalDisplacementsLog(ISubdomain_v2 subdomain, int numWatchedDofs, DofTable watchedDofs)
        {
            this.subdomain = subdomain;
            this.watchedDofs = watchedDofs;
            this.displacements = new double[numWatchedDofs];
        }

        public DateTime EndTime { get; private set; }

        public double GetDisplacementAt(INode node, DOFType dofType) => displacements[watchedDofs[node, dofType]];

        public DateTime StartTime { get; private set; }

        public void StoreResults(DateTime startTime, DateTime endTime, IVectorView solution)
        {
            StartTime = startTime;
            EndTime = endTime;

            foreach ((INode node, DOFType dofType, int dofIdx) in watchedDofs)
            {
                int globalDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                displacements[dofIdx] = solution[globalDofIdx];
            }
        }

        public class Factory : ILogFactory_v2
        {
            private readonly ISubdomain_v2 subdomain;
            private readonly DofTable watchedDofs = new DofTable();
            private int numWatchedDofs = 0;

            public Factory(ISubdomain_v2 subdomain) => this.subdomain = subdomain;

            public IAnalyzerLog_v2[] CreateLogs()
                => new IAnalyzerLog_v2[] { new TotalDisplacementsLog(subdomain, numWatchedDofs, watchedDofs) };

            public void WatchDof(INode node, DOFType dofType)
            {
                bool isAdded = watchedDofs.TryAdd(node, dofType, numWatchedDofs++);
                if (!isAdded) throw new ArgumentException(
                    $"At node {node.ID}, {dofType} is set to be watched more than once");
            }
        }
    }
}
