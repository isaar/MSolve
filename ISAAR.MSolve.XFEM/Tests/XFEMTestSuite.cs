using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Tests.DofOrdering;
using ISAAR.MSolve.XFEM.Tests.GRACM;
using ISAAR.MSolve.XFEM.Tests.Khoei;
using ISAAR.MSolve.XFEM.Tests.Subdomains;
using ISAAR.MSolve.XFEM.Tests.Visualization;

namespace ISAAR.MSolve.XFEM.Tests
{
    public class XFEMTestSuite
    {
        public static void RunAll()
        {
            //DCB3x1.Run();
            //DCBSolvers.Run();
            //BenchmarkSolvers.Run();
            IntersectedMeshCantilever.Run();

            //ReanalysisDebugging.Run();
            //ReorderingTests.Run();
            //SubdomainTest1.Run();
            //SubdomainTest2.Run();
            //AutomaticDecompositionTest.Run();
            //TestMenkBordasSolver.Run();
        }
    }
}
