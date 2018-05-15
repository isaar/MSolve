using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Tests.DofOrdering;
using ISAAR.MSolve.XFEM.Tests.GRACM;
using ISAAR.MSolve.XFEM.Tests.MenkBordas;

namespace ISAAR.MSolve.XFEM.Tests
{
    public class XFEMTestSuite
    {
        public static void RunAll()
        {
            DCBSolvers.Run();
            //ReanalysisDebugging.Run();
            //ReorderingTests.Run();
            //SubdomainTest1.Run();
            //SubdomainTest2.Run();
            //TestMenkBordasCG.Run();
        }
    }
}
