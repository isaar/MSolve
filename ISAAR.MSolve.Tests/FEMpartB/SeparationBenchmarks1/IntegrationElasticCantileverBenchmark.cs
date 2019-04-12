using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;

namespace ISAAR.MSolve.Tests.FEMpartB.SeparationBenchmarks1
{
    public static class IntegrationElasticCantileverBenchmark //checked
    {
        //checked: apotelesmata C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\APOTELESMATA_MS_hexa8_cantilever_nea\IntegrationElasticCantileverBenchmark RunExample
        
        //Origin opou htan checked branch example/ms_development_nl_elements_merge
        //modifications: egine v2
        public static TotalDisplacementsPerIterationLog RunExample()
        {
            //VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            int subdomainID = 1;  model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
            HexaCantileverBuilder_copyMS_222(model, 0.00219881744271988174427);

            //model.ConnectDataStructures();

            //var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            //linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.Subdomains[0].Forces);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            //var solver = new SolverSkyline(linearSystems[subdomainID]);
            //var linearSystemsArray = new[] { linearSystems[subdomainID] };

            //var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            //var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.MaxIterationsPerIncrement = 100;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
            var watchDofs = new Dictionary<int, int[]>();
            watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 47 });
            var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
            childAnalyzer.TotalDisplacementsPerIterationLog = log1;
           
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            return log1;
        }

        public static void HexaCantileverBuilder_copyMS_222(Model model, double load_value)
        {
            //Origin: ParadeigmataElegxwnBuilder.HexaCantileverBuilder_copyMS_222(Model model, double load_value)

            IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderNonLinear();

            IContinuumMaterial3DDefGrad material1 = new MicrostructureDefGrad3D(homogeneousRveBuilder1,
                m => (new SkylineSolver.Builder()).BuildSolver(m), false, 1);

            double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
            {0.250000,-0.250000,-1.000000},
            {-0.250000,0.250000,-1.000000},
            {0.250000,0.250000,-1.000000},
            {-0.250000,-0.250000,-0.500000},
            {0.250000,-0.250000,-0.500000},
            {-0.250000,0.250000,-0.500000},
            {0.250000,0.250000,-0.500000},
            {-0.250000,-0.250000,0.000000},
            {0.250000,-0.250000,0.000000},
            {-0.250000,0.250000,0.000000},
            {0.250000,0.250000,0.000000},
            {-0.250000,-0.250000,0.500000},
            {0.250000,-0.250000,0.500000},
            {-0.250000,0.250000,0.500000},
            {0.250000,0.250000,0.500000},
            {-0.250000,-0.250000,1.000000},
            {0.250000,-0.250000,1.000000},
            {-0.250000,0.250000,1.000000},
            {0.250000,0.250000,1.000000}};

            int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
            {2,12,11,9,10,8,7,5,6},
            {3,16,15,13,14,12,11,9,10},
            {4,20,19,17,18,16,15,13,14}, };

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

            }

            // orismos elements 
            Element e1;
            int subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinearDefGrad(material1, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
            }

            // constraint vashh opou z=-1
            for (int k = 1; k < 5; k++)
            {
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Z });
            }

            // fortish korufhs
            Load load1;
            for (int k = 17; k < 21; k++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[k],
                    DOF = DOFType.X,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }
    }
}
