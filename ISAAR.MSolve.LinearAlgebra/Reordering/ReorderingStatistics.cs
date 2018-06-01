using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    public class ReorderingStatistics
    {
        public ReorderingStatistics(int supFactorizedNumNonZeros, int numMovedDenseRows)
        {
            this.SupFactorizedNumNonZeros = supFactorizedNumNonZeros;
            this.NumMovedDenseRows = numMovedDenseRows;
        }

        public int SupFactorizedNumNonZeros { get; }
        public int NumMovedDenseRows { get; }
    }
}
