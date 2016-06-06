using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolver
    {
        void Initialize();
        void Solve();
    }
}
