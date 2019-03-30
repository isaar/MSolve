using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Problems.Structural.Elements;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;

namespace ISAAR.MSolve.SamplesConsole
{
    public class TrussExample
    {
        private const int subdomainID = 0;

        public static IList<Node_v2> CreateNodes()
        {
            IList<Node_v2> nodes = new List<Node_v2>();
            Node_v2 node1 = new Node_v2 { ID = 1, X = 0, Y = 0 };
            Node_v2 node2 = new Node_v2 { ID = 2, X = 0, Y = 40 };
            Node_v2 node3 = new Node_v2 { ID = 3, X = 40, Y = 40 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            return nodes;
        }

        public static void Run()
        {
            double youngMod = 10e6;
            //double poisson = 0.3;
            double loadX = 500;
            double loadY = 300;
            double sectionArea = 1.5;

            IList<Node_v2> nodes = TrussExample.CreateNodes();

            var model = new Model_v2();

            model.SubdomainsDictionary.Add(0, new Subdomain_v2(subdomainID));

            for (int i = 0; i < nodes.Count; i++)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = DOFType.Y });


            var element1 = new Element_v2()
            {
                ID = 1,
                ElementType = new Rod2D_v2(youngMod) { Density = 1, SectionArea = sectionArea }
            };

            var element2 = new Element_v2()
            {
                ID = 2,
                ElementType = new Rod2D_v2(youngMod) { Density = 1, SectionArea = sectionArea }
            };

            element1.AddNode(model.NodesDictionary[1]);
            element1.AddNode(model.NodesDictionary[3]);

            element2.AddNode(model.NodesDictionary[2]);
            element2.AddNode(model.NodesDictionary[3]);

            model.ElementsDictionary.Add(element1.ID, element1);
            model.ElementsDictionary.Add(element2.ID, element2);

            model.SubdomainsDictionary[subdomainID].Elements.Add(element1);
            model.SubdomainsDictionary[subdomainID].Elements.Add(element2);

            model.Loads.Add(new Load_v2() { Amount = loadX, Node = model.NodesDictionary[3], DOF = DOFType.X });
            model.Loads.Add(new Load_v2() { Amount = loadY, Node = model.NodesDictionary[3], DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural_v2(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Output requests
            var logFactory = new TotalDisplacementsLog.Factory(model.SubdomainsDictionary[subdomainID]);
            logFactory.WatchDof(model.NodesDictionary[3], DOFType.X);
            logFactory.WatchDof(model.NodesDictionary[3], DOFType.Y);
            childAnalyzer.LogFactories[subdomainID] = logFactory;

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Print output
            var logger = (TotalDisplacementsLog)(childAnalyzer.Logs[subdomainID][0]); //There is a list of logs for each subdomain and we want the first one
            double ux = logger.GetDisplacementAt(model.NodesDictionary[3], DOFType.X);
            double uy = logger.GetDisplacementAt(model.NodesDictionary[3], DOFType.Y);
            Console.WriteLine($"Displacements of Node 3: Ux = {ux}, Uy = {uy}");
        }
    }
}
