using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IPCCoefficientsProvider : IStochasticCoefficientsProvider
    {
        int OrderM { get; }
        int OrderP { get; }
        int CurrentOrder { get; set; }
        int ExpansionOrder { get; }
        int NoOfMatrices { get; }
        IPolynomialChaosCoefficients Calculator { get; }
    }
}
