using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestLibs;

namespace ISAAR.MSolve.SamplesConsole.Linear_Algebra
{
    class TestMKLInstallation
    {
        public static void Main()
        {
            ComputeNETexample.RunBlasExample1();
            BLASlvl1Tests.Run();
        }
    }
}
