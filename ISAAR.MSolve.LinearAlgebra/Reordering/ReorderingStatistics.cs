using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    public class ReorderingStatistics
    {
        public ReorderingStatistics(int factorizedNumNonZeros)
        {
            this.FactorizedNumNonZeros = factorizedNumNonZeros;
        }

        public int FactorizedNumNonZeros { get; }
    }
}
