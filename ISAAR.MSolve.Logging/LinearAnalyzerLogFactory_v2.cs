using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.Logging
{
    public class LinearAnalyzerLogFactory_v2 : ILogFactory_v2
    {
        private readonly int[] dofs;
        private readonly Element_v2[] stressElements, forceElements;

        public LinearAnalyzerLogFactory_v2(int[] dofs, Element_v2[] stressElements, Element_v2[] dofElements)
        {
            this.dofs = dofs;
            this.stressElements = stressElements;
            this.forceElements = dofElements;
        }

        public LinearAnalyzerLogFactory_v2(int[] dofs) : this(dofs, new Element_v2[0], new Element_v2[0])
        {
        }

        public IAnalyzerLog_v2[] CreateLogs()
        {
            var l = new List<IAnalyzerLog_v2>();
            l.Add(new DOFSLog_v2(dofs));
            if (stressElements.Length > 0)
                l.Add(new StressesLog_v2(stressElements));
            if (forceElements.Length > 0)
                l.Add(new ForcesLog_v2(forceElements));
            return l.ToArray();
        }
    }
}
