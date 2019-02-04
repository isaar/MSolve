using System.IO;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Logging;
using ISAAR.MSolve.Optimization.Structural.Topology.SIMP;
using Xunit;

namespace ISAAR.MSolve.Optimization.Structural.Tests.Topology.SIMP
{
    public static class TopologySimp99LinesTests
    {
        [Fact]
        private static void TestMbbBeam()
        {
            // Run the optimization
            (double compliance, Matrix densities, ObjectiveFunctionLogger logger) =
                TopologySimp99Lines.TopologyOptimization(60, 20, 0.5, 3.0, 1.5);

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(new double[]
            {
                1007.0221077739, 579.5597996395, 412.5077544968, 343.9883378073, 322.0344413300,
                 308.4706441855, 298.5663741370, 289.0954086264, 280.6621712319, 272.6552677735,
                 265.3879215834, 257.3357660875, 248.5522595584, 239.6383605602, 232.2492910000,
                 226.8881259681, 222.7741557552, 219.1685979726, 215.7808219590, 212.4871275855,
                 209.8967509607, 208.4453911762, 207.6186619678, 206.8286073982, 206.2165460520,
                 205.6797624859, 205.3860265680, 205.2792282918, 205.2036531370, 205.1364112880,
                 205.0504324839, 204.9721946517, 204.8797879878, 204.7946202261, 204.6801287848,
                 204.5780757746, 204.5108212174, 204.4584196387, 204.4346922851, 204.4153156607,
                 204.4003966915, 204.3651643531, 204.3433804827, 204.3182686879, 204.2806455407,
                 204.2370937897, 204.2080892859, 204.1580543438, 204.1192294874, 204.0831242101,
                 204.0464827613, 204.0069076953, 203.9660031619, 203.9233722869, 203.8840008611,
                 203.8439420170, 203.8088419596, 203.7829012778, 203.7561308102, 203.7315369931,
                 203.7126382407, 203.6960699182, 203.6704306112, 203.6601799424, 203.6451914557,
                 203.6217600500, 203.5993675987, 203.5830253410, 203.5567458841, 203.5422949406,
                 203.5355981337, 203.5207013382, 203.5124652445, 203.5005096659, 203.4899429205,
                 203.4860402865, 203.4764533125, 203.4670241037, 203.4591931267, 203.4502780887,
                 203.4481317111, 203.4370046665, 203.4344554090, 203.4178128204, 203.4131747941,
                 203.3943491471, 203.3892304549, 203.3808469626, 203.3713921064, 203.3578766315,
                 203.3419123647, 203.3334938592, 203.3237563148, 203.3060616186
            });
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-10));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_MBBbeam_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Matrix densitiesExpected = reader.ReadFile(densitiesPath);
            Assert.True(densitiesExpected.Equals(densities, 1E-9));
        }
    }
}
