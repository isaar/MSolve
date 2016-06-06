using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IStaticProvider : IAnalyzerProvider
    {
        void CalculateMatrix(ISolverSubdomain subdomain);
        //void CalculateMatrices();
    }
}
