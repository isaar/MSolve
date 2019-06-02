using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Logging.DomainDecomposition
{
    public class DomainDecompositionLogger : IDomainDecompositionLogger
    {
        private readonly string plotDirectoryPath;
        private int analysisStep;

        public DomainDecompositionLogger(string plotDirectoryPath) //TODO: make sure the path does not end in "\"
        {
            this.plotDirectoryPath = plotDirectoryPath;
            analysisStep = 0;
        }

        public void PlotSubdomains(IStructuralModel model)
        {
            var writer = new MeshPartitionWriter();
            writer.WriteSubdomainElements($"{plotDirectoryPath}\\subdomains_{analysisStep}.vtk", model);
            writer.WriteBoundaryNodes($"{plotDirectoryPath}\\boundary_nodes_{analysisStep}.vtk", model);
            ++analysisStep;
        }
    }
}
