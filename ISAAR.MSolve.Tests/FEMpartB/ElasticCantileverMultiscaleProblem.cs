using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Tests.FEMpartB.SeparationBenchmarks1;
using Xunit;


namespace ISAAR.MSolve.Tests.FEMpartB
{
    public static class ElasticCantileverMultiscaleProblem
    {
        [Fact]
        public static void CheckElasticMultiscaleCantilever()
        {
            IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
            TotalDisplacementsPerIterationLog computedDisplacements = IntegrationElasticCantileverBenchmark.RunExample();
            Assert.True( AreDisplacementsSame(expectedDisplacements, computedDisplacements));

        }



        private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements()
        {
            var expectedDisplacements = new Dictionary<int, double>[9]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES            

            expectedDisplacements[0] = new Dictionary<int, double> {
    { 0,3.335700883958350000e-02 }, {11,-2.632079302287960000e-02 }, {23,-4.942856730039380000e-02 }, {35,-6.269433016162970200e-02 }, {47,-6.765615287120199700e-02 }};
            expectedDisplacements[1] = new Dictionary<int, double> {
    { 0,3.450997936135530300e-02 }, {11,-2.739535658161109900e-02 }, {23,-5.634615891128919700e-02 }, {35,-8.131311540906950600e-02 }, {47,-1.019128163215870000e-01 }};
            expectedDisplacements[2] = new Dictionary<int, double> {
    { 0,3.433216808183409800e-02 }, {11,-2.726099954481620000e-02 }, {23,-5.629518934457999900e-02 }, {35,-8.199981263488670400e-02 }, {47,-1.039303808027040000e-01 }};
            expectedDisplacements[3] = new Dictionary<int, double> {
    { 0,3.431257880330890200e-02 }, {11,-2.724173809701519900e-02 }, {23,-5.624825754899259700e-02 }, {35,-8.192386243981529500e-02 }, {47,-1.038312844223340000e-01 }};
            expectedDisplacements[4] = new Dictionary<int, double> {
    { 0,6.894482083872839600e-02 }, {11,-5.475067582655650200e-02 }, {23,-1.173163323056880000e-01 }, {35,-1.790583221645240000e-01 }, {47,-2.376077026742320100e-01 }};
            expectedDisplacements[5] = new Dictionary<int, double> {
    { 0,6.972463521590920000e-02 }, {11,-5.540512678394360300e-02 }, {23,-1.220341983649240000e-01 }, {35,-1.920604720743059900e-01 }, {47,-2.620585115820520100e-01 }};
            expectedDisplacements[6] = new Dictionary<int, double> {
    { 0,6.858059919522070700e-02 }, {11,-5.432730995597089700e-02 }, {23,-1.192782599647590100e-01 }, {35,-1.873563500914020000e-01 }, {47,-2.554697448169410100e-01 }};
            expectedDisplacements[7] = new Dictionary<int, double> {
    { 0,6.835175024769160600e-02 }, {11,-5.410392477418309700e-02 }, {23,-1.186661258178350100e-01 }, {35,-1.861855064114160100e-01 }, {47,-2.536413732588089800e-01 }};
            expectedDisplacements[8] = new Dictionary<int, double> {
    { 0,6.834878138258780600e-02 }, {11,-5.410102312471470200e-02 }, {23,-1.186583648525040000e-01 }, {35,-1.861711431280629900e-01 }, {47,-2.536196564358649800e-01 }};


            return expectedDisplacements;
        }

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements,
            TotalDisplacementsPerIterationLog computedDisplacements)
        {
            int subdomainID = 1;
            var comparer = new ValueComparer(1E-13);
            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
    }
}
