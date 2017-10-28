using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Problems.Structural;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.SamplesConsole
{
    class Program
    {
        private static void SolveBuildingInNoSoilSmall()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, 1, 4, false, false);
            model.Loads.Add(new Load() { Amount = -100, Node = model.Nodes[21], DOF = DOFType.X });
            model.ConnectDataStructures();

            var linearSystem = new SkylineLinearSystem(1, model.SubdomainsDictionary[1].Forces);
            SolverFBSubstitution solver = new SolverFBSubstitution(linearSystem);
            ProblemStructural provider = new ProblemStructural(model);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystem);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystem);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        private static void SolveBuildingInNoSoilSmallDynamic()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, 1, 4, false, false);
            model.ConnectDataStructures();

            var linearSystem = new SkylineLinearSystem(1, model.SubdomainsDictionary[1].Forces);
            SolverFBSubstitution solver = new SolverFBSubstitution(linearSystem);
            ProblemStructural provider = new ProblemStructural(model);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystem);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, linearSystem, 0.5, 0.25, 0.01, 0.1);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }


        static void Main(string[] args)
        {
            SolveBuildingInNoSoilSmall();
        }
    }
}
