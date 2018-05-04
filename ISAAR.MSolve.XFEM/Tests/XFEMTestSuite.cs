using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Tests.DofOrdering;
using ISAAR.MSolve.XFEM.Tests.GRACM;

namespace ISAAR.MSolve.XFEM.Tests
{
    public class XFEMTestSuite
    {
        public static void RunAll()
        {
            ReorderingTests.Run();
            //DCBSolvers.Run();
        }
    }
}
