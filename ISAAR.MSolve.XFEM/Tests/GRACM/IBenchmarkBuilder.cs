using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    interface IBenchmarkBuilder
    {
        string TimingPath { get; }
        IBenchmark BuildBenchmark();
    }
}
