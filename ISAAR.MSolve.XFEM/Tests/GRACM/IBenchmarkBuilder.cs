using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    interface IBenchmarkBuilder
    {
        int MaxIterations { get; set; }
        string TimingOutputDirectory { get; }
        IBenchmark BuildBenchmark();
    }
}
