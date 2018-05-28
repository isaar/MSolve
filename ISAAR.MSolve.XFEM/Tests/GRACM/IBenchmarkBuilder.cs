using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    interface IBenchmarkBuilder
    {
        string TimingOutputDirectory { get; }
        IBenchmark BuildBenchmark();
    }
}
