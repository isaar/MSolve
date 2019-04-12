using ISAAR.MSolve.Tests.SupportiveClasses;
using MGroup.Stochastic;
using MGroup.Stochastic.Structural;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public static class CantileverBeamWithStochasticMaterial
    {
        [Fact]
        public static void Solve()
        {
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var domainMapper = new CantileverStochasticDomainMapper(new[] { 0d, 0d, 0d });
            var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var m = new MonteCarlo(iterations, evaluator, evaluator);
            m.Evaluate();
            double[] expectedMeanValueResponse = new double[]
            {
                1.5999674517697445E-08, -2.399309224401548E-08
            };

            for (int i = 0; i < expectedMeanValueResponse.Length; i++)
                Assert.Equal(expectedMeanValueResponse[i], m.MeanValueResponse[i], 10);
            double[] expectedStandardDeviationResponse = new double[]
            {
                4.3180461976960273E-12, 2.6374139668652252E-12
            };

            for (int i = 0; i < expectedStandardDeviationResponse.Length; i++)
                Assert.Equal(expectedStandardDeviationResponse[i], m.StandardDeviationResponse[i], 10);

        }
    }
}
