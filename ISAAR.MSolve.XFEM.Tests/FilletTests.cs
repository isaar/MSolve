using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests
{
    public static class FilletTests
    {
        [Theory]
        [InlineData(true)]
        [InlineData(false)]
        public static void TestPropagation(bool rigidBoundaryConditions)
        {
            string meshPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\fillet_1272dofs.msh";
            //string propagationPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Propagation\crack_growth.txt";
            //string plotPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Plots";
            //string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Timing";

            double growthLength = 5; // mm. Must be sufficiently larger than the element size.
            var builder = new FilletBenchmark.Builder(meshPath, growthLength);
            //builder.WritePropagation = writePropagationPath;
            builder.HeavisideEnrichmentTolerance = 0.01;
            builder.RigidBCs = rigidBoundaryConditions;

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
            //SuiteSparseSolver solver = new SuiteSparseSolver.Builder().BuildSolver(benchmark.Model);
            SkylineSolver solver = new SkylineSolver.Builder().BuildSolver(benchmark.Model);
            benchmark.Analyze(solver);

            CheckPropagationPath(benchmark.Crack.CrackPath, rigidBoundaryConditions);
        }

        //TODO: The reference crack paths were taken from solving this benchmark with previous versions of my code. Although
        //      they match visually the expected crack paths in publications, the exact ones must be searched.
        private static void CheckPropagationPath(IReadOnlyList<CartesianPoint> crackPath, bool rigidBCs)
        {
            var expectedPath = new CartesianPoint[15];
            if (rigidBCs)
            {
                expectedPath[0] = new CartesianPoint(150.000000000000, 95.000000000000);
                expectedPath[1] = new CartesianPoint(155.000000000000, 95.000000000000);
                expectedPath[2] = new CartesianPoint(159.928845175266, 94.1594732376293);
                expectedPath[3] = new CartesianPoint(164.895272288676, 93.5810255501270);
                expectedPath[4] = new CartesianPoint(169.873369876570, 93.1135382183093);
                expectedPath[5] = new CartesianPoint(174.865946489105, 92.8411804242900);
                expectedPath[6] = new CartesianPoint(179.863172733028, 92.6746575518206);
                expectedPath[7] = new CartesianPoint(184.863046509768, 92.6391298631669);
                expectedPath[8] = new CartesianPoint(189.862887217122, 92.6790410281461);
                expectedPath[9] = new CartesianPoint(194.861153856214, 92.8106867821075);
                expectedPath[10] = new CartesianPoint(199.854887645723, 93.0609323405507);
                expectedPath[11] = new CartesianPoint(204.846200988297, 93.3555356607026);
                expectedPath[12] = new CartesianPoint(209.835693764484, 93.6795138025205);
                expectedPath[13] = new CartesianPoint(214.822110231868, 94.0478218445237);
                expectedPath[14] = new CartesianPoint(219.810510335299, 94.3882104578577);
            }
            else
            {
                expectedPath[0] =  new CartesianPoint(150.000000000000, 95.0000000000000);
                expectedPath[1] =  new CartesianPoint(155.000000000000, 95.0000000000000);
                expectedPath[2] =  new CartesianPoint(159.546278695233, 92.9188104302414);
                expectedPath[3] =  new CartesianPoint(163.855055966478, 90.3821909422607);
                expectedPath[4] =  new CartesianPoint(167.842709039709, 87.3657985660029);
                expectedPath[5] =  new CartesianPoint(171.568844458286, 84.0317447508364);
                expectedPath[6] =  new CartesianPoint(175.075157758010, 80.4672297706211);
                expectedPath[7] =  new CartesianPoint(178.428998912345, 76.7588946304959);
                expectedPath[8] =  new CartesianPoint(181.635027717824, 72.9220496638023);
                expectedPath[9] =  new CartesianPoint(184.762601328988, 69.0209878459438);
                expectedPath[10] = new CartesianPoint(187.806841672205, 65.0545535867224);
                expectedPath[11] = new CartesianPoint(190.768299071507, 61.0259347077364);
                expectedPath[12] = new CartesianPoint(193.591488900770, 56.8992412143998);
                expectedPath[13] = new CartesianPoint(196.150616188081, 52.6037939120413);
                expectedPath[14] = new CartesianPoint(198.453603879810, 48.1657484310051);
            }

            Assert.Equal(expectedPath.Length, crackPath.Count);
            int precision = 2;
            for (int i = 0; i < expectedPath.Length; ++i)
            {
                Assert.Equal(expectedPath[i].X, crackPath[i].X, precision);
                Assert.Equal(expectedPath[i].Y, crackPath[i].Y, precision);
            }
        }
    }
}
