using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    interface IBenchmarkBuilder
    {
        string LsmOutputDirectory { get; set; }
        string TimingPath { get; }
        IBenchmark BuildBenchmark();
    }
}
