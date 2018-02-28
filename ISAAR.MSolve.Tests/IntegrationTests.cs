using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class IntegrationTests
    {
        [Fact]
        public void TestSolveHexaCantileverBeam()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            HexaSimpleCantileverBeam.MakeCantileverBeam(model, 0, 0, 0, model.NodesDictionary.Count + 1, model.ElementsDictionary.Count + 1, 1);
            
            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[16], DOF = DOFType.Z });
            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[17], DOF = DOFType.Z });
            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[18], DOF = DOFType.Z });
            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[19], DOF = DOFType.Z });

            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            //NewtonRaphsonNonLinearAnalyzer analyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems, provider, 10, 48);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            double[] expectedDisplacements = new double[]
            {
                -0.0000025899520106, -0.0000004898560318, -0.0000031099520106, -0.0000025899520106, 0.0000004898560318,
                -0.0000031099520106, 0.0000025899520106, 0.0000004898560318, -0.0000031099520106, 0.0000025899520106,
                -0.0000004898560318, -0.0000031099520106, -0.0000045673419128, -0.0000002423136749, -0.0000107872459340,
                -0.0000045673419128, 0.0000002423136749, -0.0000107872459340, 0.0000045673419128, 0.0000002423136749,
                -0.0000107872459340, 0.0000045673419128, -0.0000002423136749, -0.0000107872459340, -0.0000057299058132,
                -0.0000001253780263, -0.0000216044936601, -0.0000057299058132, 0.0000001253780263, -0.0000216044936601,
                0.0000057299058132, 0.0000001253780263, -0.0000216044936601, 0.0000057299058132, -0.0000001253780263,
                -0.0000216044936601, -0.0000061325564473, -0.0000000425738760, -0.0000339869559207, -0.0000061325564473,
                0.0000000425738760, -0.0000339869559207, 0.0000061325564473, 0.0000000425738760, -0.0000339869559207,
                0.0000061325564473, -0.0000000425738760, -0.0000339869559207
            };

            for (int i = 0; i < expectedDisplacements.Length; i++)
                Assert.Equal(expectedDisplacements[i], linearSystems[1].Solution[i],10);
        }
    }
}
