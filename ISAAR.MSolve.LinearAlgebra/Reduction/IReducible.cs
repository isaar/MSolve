using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Reduction
{

    public delegate double ProcessEntry(double entry, double aggregator);
    public delegate double ProcessZeros(int numZeros, double aggregator);
    public delegate double Finalize(double aggregator);

    public interface IReducible
    {
        double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize);
    }
}