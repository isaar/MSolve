using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Logging.DomainDecomposition
{
    public class DomainDecompositionLogger : IDomainDecompositionLogger
    {
        private readonly string plotDirectoryPath;
        private readonly bool shuffleSubdomainColors;
        private int analysisStep;

        //TODO: make sure the path does not end in "\"
        public DomainDecompositionLogger(string plotDirectoryPath, bool shuffleSubdomainColors = false) 
        {
            this.plotDirectoryPath = plotDirectoryPath;
            this.shuffleSubdomainColors = shuffleSubdomainColors;
            analysisStep = 0;
        }

        public void PlotSubdomains(IStructuralModel model)
        {
            var writer = new MeshPartitionWriter(shuffleSubdomainColors);
            writer.WriteSubdomainElements($"{plotDirectoryPath}\\subdomains_{analysisStep}.vtk", model);
            writer.WriteBoundaryNodes($"{plotDirectoryPath}\\boundary_nodes_{analysisStep}.vtk", model);
            ++analysisStep;
        }
    }
}
