using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests
{
    public static class HolesTests
    {
        [Fact]
        public static void TestPropagation()
        {
            string meshPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\holes_4442dofs.msh";
            //string propagationPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Propagation\crack_growth.txt";
            //string plotPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Plots";
            //string timingPath = @"C:\Users\Serafeim\Desktop\GRACM\Benchmark_Fillet\Timing";

            double growthLength = 1.0; // mm. Must be sufficiently larger than the element size.
            var builder = new HolesBenchmark.Builder(meshPath, growthLength);
            builder.HeavisideEnrichmentTolerance = 0.12;

            // Usually should be in [1.5, 2.5). The J-integral radius must be large enough to at least include elements around
            // the element that contains the crack tip. However it must not be so large that an element intersected by the 
            // J-integral contour is containes the previous crack tip. Thus the J-integral radius must be sufficiently smaller
            // than the crack growth length.
            builder.JintegralRadiusOverElementSize = 2.0;

            // If you modify the following two parameters significantly, then you will need to redefine which nodes are expected 
            // to be enriched.
            builder.TipEnrichmentRadius = 0.5;
            builder.BC = HolesBenchmark.BoundaryConditions.BottomConstrainXDisplacementY_TopConstrainXDisplacementY;

            builder.MaxIterations = 12;
            //builder.LeftLsmPlotDirectory = plotLSM ? plotPathLeft : null;
            //builder.RightLsmPlotDirectory = plotLSM ? plotPathRight : null;

            // Run the benchmark
            HolesBenchmark benchmark = builder.BuildBenchmark();
            benchmark.InitializeModel();
            //SuiteSparseSolver solver = new SuiteSparseSolver.Builder().BuildSolver(benchmark.Model);
            var solverBuilder = new SkylineSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithCSparseAmd());
            SkylineSolver solver = solverBuilder.BuildSolver(benchmark.Model);
            benchmark.Analyze(solver);

            CheckCrackPropagationPath(benchmark);
        }

        //TODO: The reference crack paths were taken from solving this benchmark with previous versions of my code. Although
        //      they match visually the expected crack paths in publications, the exact ones must be searched.
        private static void CheckCrackPropagationPath(HolesBenchmark benchmark)
        {
            // Left crack
            var expectedPathLeft = new CartesianPoint[14];
            expectedPathLeft[0]  = new CartesianPoint(0, 2.85);
            expectedPathLeft[1]  = new CartesianPoint(1, 2.85);
            expectedPathLeft[2]  = new CartesianPoint(1.99995199459598, 2.85979839290530);
            expectedPathLeft[3]  = new CartesianPoint(2.99929868004322, 2.89593981354163);
            expectedPathLeft[4]  = new CartesianPoint(3.97604603958355, 2.68154626039464);
            expectedPathLeft[5]  = new CartesianPoint(4.90503431488344, 2.31143710844710);
            expectedPathLeft[6]  = new CartesianPoint(5.86904738124316, 2.04558224599815);
            expectedPathLeft[7]  = new CartesianPoint(6.84852018861039, 1.84400590902343);
            expectedPathLeft[8]  = new CartesianPoint(7.84391778239331, 1.74817464955572);
            expectedPathLeft[9]  = new CartesianPoint(8.84344037927620, 1.71727839864243);
            expectedPathLeft[10] = new CartesianPoint(9.83984063952885, 1.80205175158095);
            expectedPathLeft[11] = new CartesianPoint(10.8276024105833, 1.95802188860506);
            expectedPathLeft[12] = new CartesianPoint(11.7939950361318, 2.21509248841279);
            expectedPathLeft[13] = new CartesianPoint(12.7922235636750, 2.27458876391790);

            Assert.Equal(expectedPathLeft.Length, benchmark.LeftCrack.CrackPath.Count);
            int precisionLeft = 2;
            for (int i = 0; i < expectedPathLeft.Length; ++i)
            {
                Assert.Equal(expectedPathLeft[i].X, benchmark.LeftCrack.CrackPath[i].X, precisionLeft);
                Assert.Equal(expectedPathLeft[i].Y, benchmark.LeftCrack.CrackPath[i].Y, precisionLeft);
            }

            // Right crack
            var expectedPathRight = new CartesianPoint[14];
            expectedPathRight[0]  = new CartesianPoint(20, 7.15);
            expectedPathRight[1]  = new CartesianPoint(19, 7.15);
            expectedPathRight[2]  = new CartesianPoint(18.0108334835354, 7.00320217063868);
            expectedPathRight[3]  = new CartesianPoint(17.0296410112692, 6.81017022414205);
            expectedPathRight[4]  = new CartesianPoint(16.0310132864364, 6.8625407043302);
            expectedPathRight[5]  = new CartesianPoint(15.1290705380152, 7.29439632669779);
            expectedPathRight[6]  = new CartesianPoint(14.1597696631581, 7.54027396700138);
            expectedPathRight[7]  = new CartesianPoint(13.1820518706078, 7.75019756770084);
            expectedPathRight[8]  = new CartesianPoint(12.1871075395841, 7.85062554270676);
            expectedPathRight[9]  = new CartesianPoint(11.1878028931389, 7.81333984044864);
            expectedPathRight[10] = new CartesianPoint(10.1894240274647, 7.75642199987822);
            expectedPathRight[11] = new CartesianPoint(9.20526987140180, 7.57910703405718);
            expectedPathRight[12] = new CartesianPoint(8.22923499887833, 7.36149286260210);
            expectedPathRight[13] = new CartesianPoint(7.24370720476692, 7.19197894861886);

            Assert.Equal(expectedPathRight.Length, benchmark.RightCrack.CrackPath.Count);
            int precisionRight = 2;
            for (int i = 0; i < expectedPathRight.Length; ++i)
            {
                Assert.Equal(expectedPathRight[i].X, benchmark.RightCrack.CrackPath[i].X, precisionRight);
                Assert.Equal(expectedPathRight[i].Y, benchmark.RightCrack.CrackPath[i].Y, precisionRight);
            }
        }
    }
}
