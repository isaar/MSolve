using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Logging.DomainDecomposition
{
    public interface IDomainDecompositionLogger
    {
        void PlotSubdomains(IStructuralModel model);
    }
}
