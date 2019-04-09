using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
//using ISAAR.MSolve.SamplesConsole.DdmBenchmarks1;
using ISAAR.MSolve.SamplesConsole.Solvers;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Tests.FEMpartB;
using MGroup.Stochastic;
using MGroup.Stochastic.Structural;
using MGroup.Stochastic.Structural.Example;


namespace ISAAR.MSolve.SamplesConsole
{
    class Program
    {
        private const int subdomainID = 0;

        static void Main(string[] args)
        {
            //SolveBuildingInNoSoilSmall();
            //TrussExample.Run();
            //FEM.Cantilever2D.Run();
            //FEM.Cantilever2DPreprocessor.Run();
            //FEM.WallWithOpenings.Run();
            //SeparateCodeCheckingClass.Check06();
            //SolveBuildingInNoSoilSmall_v2();
            //SolveBuildingInNoSoilSmallDynamic_v2();
            //SolveStochasticMaterialBeam2DWithBruteForceMonteCarlo();
            //CNTExamples.CNT_4_4_DisplacementControl_v2();
            //CNTExamples.CNT_4_4_NewtonRaphson_v2();
            //Tests.FEM.Shell8andCohesiveNonLinear.RunTest_v2();
            //AppliedDisplacementExample.Run();

            //Logging.PrintForceDisplacementCurve.CantileverBeam2DCorotationalLoadControl();

            //SuiteSparseBenchmarks.MemoryConsumptionDebugging();
            //SolverBenchmarks.SuiteSparseMemoryConsumptionDebugging();
            //NRNLAnalyzerDevelopTest_v2.SolveDisplLoadsExample();
            //SeparateCodeCheckingClass4.Check05bStressIntegrationObje_v2_Integration();
            //SeparateCodeCheckingClass4.Check_Graphene_rve_Obje_v2_Integration();
            //IntegrationElasticCantileverBenchmark.RunExample();
            //OneRveExample.Check_Graphene_rve_serial();
            //BondSlipTest.CheckStressStrainBonSlipMaterial();
            //OneRveExample.Check_Graphene_rve_parallel();
            //LinearRves.CheckShellScaleTransitionsAndMicrostructure();
            SolveCantileverWithStochasticMaterial();
        }

        private static void SolveBuildingInNoSoilSmall()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = subdomainID });
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);
            model.Loads.Add(new Load() { Amount = -100, Node = model.Nodes[21], DOF = DOFType.X });
            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[subdomainID]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            analyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            DOFSLog log = (DOFSLog)analyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {420}, u = {log.DOFValues[420]}");
        }

        private static void SolveBuildingInNoSoilSmall_v2()
        {
            var model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding_v2(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);
            model.Loads.Add(new Load_v2() { Amount = -100, Node = model.Nodes[21], DOF = DOFType.X });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural_v2(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory_v2(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog_v2 log = (DOFSLog_v2)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");
        }

        private static void SolveBuildingInNoSoilSmallDynamic()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = subdomainID });
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);
            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[subdomainID]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, linearSystems, 0.25, 0.5, 0.01, 0.1);

            analyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            DOFSLog log = (DOFSLog)analyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {420}, u = {log.DOFValues[420]}");
        }

        private static void SolveBuildingInNoSoilSmallDynamic_v2()
        {
            var model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding_v2(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural_v2(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer_v2.Builder(model, solver, provider, childAnalyzer, 0.01, 0.1);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            NewmarkDynamicAnalyzer_v2 parentAnalyzer = parentAnalyzerBuilder.Build();

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory_v2(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog_v2 log = (DOFSLog_v2)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");

            //TODO: No loads have been defined so the result is bound to be 0.
        }

        private static void SolveStochasticMaterialBeam2DWithBruteForceMonteCarlo()
        {
            #region Beam2D Geometry Data
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 2.0e08;
            double poissonRatio = 0.3;
            double nodalLoad = 10.0;

            IStochasticMaterialCoefficientsProvider coefficientProvider = new PowerSpectrumTargetEvaluatorCoefficientsProvider(10, 0.1, .05, 20, 200, DOFType.X, 0.1, 200, 1e-10);
            StochasticElasticMaterial material = new StochasticElasticMaterial(coefficientProvider)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < 11; i++)
            {
                model.NodesDictionary.Add(i, new Node
                {
                    ID = i,
                    X = i * 1,
                    Y = 0,
                    Z = 0
                });
            }

            // Fix cantilever left end node of the model
            model.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[0].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            for (int i = 0; i < model.NodesDictionary.Count - 1; i++)
            {
                var element = new Element()
                {
                    ID = i,
                    ElementType = new Beam2DWithStochasticMaterial(material)
                    {
                        SectionArea = 1,
                        MomentOfInertia = 0.1
                    }
                };
                element.AddNode(model.NodesDictionary[i]);
                element.AddNode(model.NodesDictionary[i + 1]);
                model.ElementsDictionary.Add(i, element);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, element);
            }


            // Add nodal load values at the right end of the cantilever
            model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[model.NodesDictionary.Count - 1], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();
            #endregion

            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            MonteCarloAnalyzerWithStochasticMaterial stohasticAnalyzer =
            new MonteCarloAnalyzerWithStochasticMaterial(model, provider, parentAnalyzer, linearSystems,
                    coefficientProvider, 1, 10000);
            stohasticAnalyzer.Initialize();
            stohasticAnalyzer.Solve();

            //Assert.Equal(-2.08333333333333333e-5, stohasticAnalyzer.MonteCarloMeanValue, 8);
        }


        private static void SolveCantileverWithStochasticMaterial()
        {
            VectorExtensions.AssignTotalAffinityCount();
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var domainMapper = new CantileverStochasticDomainMapper(new[] { 0d, 0d, 0d });
            var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var m = new MonteCarlo(iterations, evaluator, evaluator);
            m.Evaluate();
        }
    }
}
