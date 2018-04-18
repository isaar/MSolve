using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Analysis
{
    // TODO: The Analyzer should be separate from the solver.
    interface ILinearStaticAnalysis
    {
        Vector Solution { get; }

        void Initialize();
        void Solve();
        void PrintSolution();
    }
}
