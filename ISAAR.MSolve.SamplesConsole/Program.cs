using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials.Interfaces;

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

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

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

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, linearSystems, 0.25, 0.5, 0.01, 0.1);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
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
            model.NodesDictionary[0].Constraints.Add(DOFType.X);
            model.NodesDictionary[0].Constraints.Add(DOFType.Y);
            model.NodesDictionary[0].Constraints.Add(DOFType.RotZ);

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
        
        static void Main(string[] args)
        {
            //SolveBuildingInNoSoilSmall();
            //TrussExample.Run();
            //FEM.Cantilever2D.Run();
            //FEM.Cantilever2DPreprocessor.Run();
            FEM.WallWithOpenings.Run();
        }
    }
}
