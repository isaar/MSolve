using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests
{
    public static class FilletTests
    {
        [Fact]
        public static void TestPropagation()
        {
            string meshPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Meshes\fillet.msh";
            //string propagationPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Propagation\crack_growth.txt";
            //string plotPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Plots";
            //string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Timing";

            double growthLength = 5; // mm. Must be sufficiently larger than the element size.
            var builder = new FilletBenchmark.Builder(growthLength, meshPath);
            //builder.WritePropagation = writePropagationPath;
            builder.HeavisideEnrichmentTolerance = 0.01;
            builder.RigidBCs = true;
            builder.NumSubdomains = 5;

            //builder.LsmPlotDirectory = plotLSM ? plotPath : null;
            builder.MaxIterations = 13;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            // Run the benchmark
            FilletBenchmark benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            SkylineSolver solver = new SkylineSolver.Builder().BuildSolver(benchmark.Model);
            benchmark.Analyze(solver);

            CheckPropagationPath(benchmark.Crack.CrackPath);
        }

        private static void CheckPropagationPath(IReadOnlyList<CartesianPoint> crackPath)
        {
            var expectedPath = new CartesianPoint[15];
            expectedPath[0] = new CartesianPoint(150.000000000000, 95.000000000000);
            expectedPath[1] = new CartesianPoint(155.000000000000, 95.000000000000);
            expectedPath[2] = new CartesianPoint(159.972400076204, 94.475369194419);
            expectedPath[3] = new CartesianPoint(164.969005754588, 94.2911636170193);
            expectedPath[4] = new CartesianPoint(169.968291638482, 94.3756659903413);
            expectedPath[5] = new CartesianPoint(174.962888574803, 94.6080480974199);
            expectedPath[6] = new CartesianPoint(179.951921636336, 94.9390301376132);
            expectedPath[7] = new CartesianPoint(184.933954583460, 95.3625247819367);
            expectedPath[8] = new CartesianPoint(189.916060436371, 95.7851608788293);
            expectedPath[9] = new CartesianPoint(194.909363596072, 96.043856753013);
            expectedPath[10] = new CartesianPoint(199.907624295387, 96.1757278327639);
            expectedPath[11] = new CartesianPoint(204.907160604311, 96.1076345456372);
            expectedPath[12] = new CartesianPoint(209.904818476420, 95.9546123467587);
            expectedPath[13] = new CartesianPoint(214.901423131894, 95.7703790248536);
            expectedPath[14] = new CartesianPoint(219.900347693701, 95.6666812652064);

            Assert.Equal(expectedPath.Length, crackPath.Count);
            int precision = 4;
            for (int i = 0; i < expectedPath.Length; ++i)
            {
                Assert.Equal(expectedPath[i].X, crackPath[i].X, precision);
                Assert.Equal(expectedPath[i].Y, crackPath[i].Y, precision);
            }
        }
    }
}
