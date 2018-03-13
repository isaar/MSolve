using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Optimization;

namespace ISAAR.MSolve.OptimizationTests
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestDE.Run();
            //TestDEConstrained.Run();
            //TestGA.Run();
            //TestPSO.Run();

            TestTruss10Benchmark.Run();
        }
    }
}
