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
    /// <summary>
    /// Tests for <see cref="TopologySimp99Lines"/>. They cover all extensions presented in the original paper
    /// "A 99 line topology optimization code written in Matlab, O. Sigmund, 1991".
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TopologySimp99LinesTests
    {
        internal static double[] ComplianceHistoryForCantilever => new double[]
        {
            429.2330474336, 254.3009164179, 178.5784401160, 140.4462442231, 122.9623005009,
            107.8651217652,  95.8505642454,  85.8912909662,  78.6566669333,  73.7832456086,
             69.9984973735,  66.6419165142,  63.8484320311,  61.9206601960,  60.7705805103,
             60.0935472538,  59.6830408718,  59.3694615068,  59.1097304745,  58.8994639921,
             58.7353645277,  58.5985347639,  58.4910112988,  58.3942898852,  58.3058184270,
             58.2017305010,  58.1012619809,  57.9937574872,  57.8831828827,  57.7912833926,
             57.7206119641,  57.6657539093,  57.6283391372,  57.6079401395,  57.5923734421,
             57.5806239245,  57.5681345244,  57.5605926804,  57.5503857773,  57.5472472631,
             57.5424407787,  57.5366686348,  57.5291346966,  57.5252901746,  57.5156589365,
             57.5063466667,  57.4946198375,  57.4805897757,  57.4609883087,  57.4463886235,
             57.4340110171,  57.4156397500,  57.3999272954,  57.3816337042,  57.3605463691,
             57.3461120267,  57.3290164141,  57.3231741894,  57.3170543899,  57.3171800351,
             57.3240324531,  57.3318729880,  57.3343124161,  57.3391126470,  57.3407987904,
             57.3437826605,  57.3483968256,  57.3489822796,  57.3477706358,  57.3518184734,
             57.3491630698
        };

        internal static double[] ComplianceHistoryForCantileverWith2LoadCases => new double[]
        {
            445.4931968905, 269.4398266551, 206.2288018075, 173.3365748825, 148.7277106959,
            121.7785005410,  96.7196190100,  80.1278584879,  71.6231799623,  67.5616293945,
             65.4659885024,  64.4160866042,  63.7772338552,  63.4050771480,  63.1379921169,
             62.9448629611,  62.7768440453,  62.6313387398,  62.4844129855,  62.3384360551,
             62.1844813279,  62.0415951767,  61.8999732289,  61.8008502754,  61.7433075384,
             61.7126519863,  61.6694558665,  61.6351251085,  61.6049798860,  61.5603070009,
             61.5210031765,  61.4866237060,  61.4521455022,  61.4056481800,  61.3716597821,
             61.3423384638,  61.3069785434,  61.2741048168,  61.2475521102,  61.2210118854,
             61.2132250082,  61.2058428092,  61.2087307437,  61.2092545297,  61.2132370930,
             61.2208462397,  61.2236218767,  61.2353478498,  61.2422782128,  61.2479828092,
             61.2550029789,  61.2621224736,  61.2620514681,  61.2674655751,  61.2665835469,
             61.2683492413,  61.2696793883,  61.2712922629,  61.2659908078,  61.2667439256,
             61.2605340692,  61.2574086575,  61.2553721768,  61.2572265868,  61.2612208612,
             61.2609843020,  61.2640822865,  61.2637731797,  61.2675474850,  61.2688982144,
             61.2735008253,  61.2769082666,  61.2785478443,  61.2787461935,  61.2827553924,
             61.2879997329
        };

        internal static double[] ComplianceHistoryForCantileverWithPassiveElements => new double[]
        {
            370.8613052581, 136.5205784542, 79.2534024658, 61.1468154986, 56.0035596005,
                 53.9558359252,  53.0247359134, 52.6095598968, 52.4033536890, 52.3002773050,
                 52.2566965195,  52.2282416664, 52.2138028187, 52.2035282075, 52.1885945575,
                 52.1825970120,  52.1689168100, 52.1615513604, 52.1553488042, 52.1499139184,
                 52.1455376426,  52.1417267820, 52.1382801344, 52.1348974348, 52.1228165997,
                 52.1256729715,  52.1242030289, 52.1137737132, 52.1176546700, 52.1084199673,
                 52.1130325742,  52.1126367052, 52.1034540016, 52.1079769128
        };

        internal static double[] ComplianceHistoryForMbbBeam => new double[]
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
        };

        [Fact]
        private static void TestCantileverBeam()
        {
            // Run the optimization
            var simp = new TopologySimp99Lines(32, 20, 0.4, 3.0, 1.2, TopologySimp99Lines.BoundaryConditions.ShortCantilever,
                TopologySimp99Lines.PassiveElements.No, TopologySimp99Lines.OptimAlgorithm.OC);
            (double compliance, Matrix densities, ObjectiveFunctionLogger logger) = simp.TopologyOptimization();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(ComplianceHistoryForCantilever);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-11));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_shortcantilever_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Matrix densitiesExpected = reader.ReadFile(densitiesPath);
            Assert.True(densitiesExpected.Equals(densities, 1E-10));
        }

        [Fact]
        private static void TestCantileverBeamWith2LoadCases()
        {
            // Run the optimization
            var simp = new TopologySimp99Lines(30, 30, 0.4, 3.0, 1.2,
                TopologySimp99Lines.BoundaryConditions.Cantilever2LoadCases, TopologySimp99Lines.PassiveElements.No, 
                TopologySimp99Lines.OptimAlgorithm.OC);
            (double compliance, Matrix densities, ObjectiveFunctionLogger logger) = simp.TopologyOptimization();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(ComplianceHistoryForCantileverWith2LoadCases);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-11));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_cantilever_loadcases_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Matrix densitiesExpected = reader.ReadFile(densitiesPath);
            Assert.True(densitiesExpected.Equals(densities, 1E-11));
        }

        [Fact]
        private static void TestCantileverBeamWithPassiveElements()
        {
            // Run the optimization
            var simp = new TopologySimp99Lines(45, 30, 0.5, 3.0, 1.5, TopologySimp99Lines.BoundaryConditions.ShortCantilever,
                TopologySimp99Lines.PassiveElements.HoleInCantilever, TopologySimp99Lines.OptimAlgorithm.OC);
            (double compliance, Matrix densities, ObjectiveFunctionLogger logger) = simp.TopologyOptimization();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(ComplianceHistoryForCantileverWithPassiveElements);
            Assert.True(expectedCompliances.Equals(
                Vector.CreateFromArray(logger.ObjectiveFunctionValues.ToArray()),
                1E-11));

            // Check the optimum element densities (design variables).
            string densitiesPath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\topology99lines_cantilever_passive_densities.txt";
            var reader = new FullMatrixReader(false, ',');
            Matrix densitiesExpected = reader.ReadFile(densitiesPath);
            Assert.True(densitiesExpected.Equals(densities, 1E-10));
        }

        [Fact]
        private static void TestMbbBeam()
        {
            // Run the optimization
            var simp = new TopologySimp99Lines(60, 20, 0.5, 3.0, 1.5, TopologySimp99Lines.BoundaryConditions.MbbBeam,
                TopologySimp99Lines.PassiveElements.No, TopologySimp99Lines.OptimAlgorithm.OC);
            (double compliance, Matrix densities, ObjectiveFunctionLogger logger) = simp.TopologyOptimization();

            // Check the history of the compliance (objective function)
            var expectedCompliances = Vector.CreateFromArray(ComplianceHistoryForMbbBeam);
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
