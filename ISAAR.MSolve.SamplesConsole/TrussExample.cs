using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.FEM.Problems.Structural.Elements;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.SamplesConsole
{
    public class TrussExample
    {
        public static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0, Y = 0 };
            Node node2 = new Node { ID = 2, X = 0, Y = 40 };
            Node node3 = new Node { ID = 3, X = 40, Y = 40 };
          
            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            return nodes;
        }

        public static void Run()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngMod = 10e6;
            double poisson = 0.3;
            double loadX = 500;
            double loadY = 300;
            double sectionArea = 1.5;

            IList<Node> nodes = TrussExample.CreateNodes();

            Model trussModel = new Model();

            trussModel.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            for (int i = 0; i < nodes.Count; i++)
            {
                trussModel.NodesDictionary.Add(i + 1, nodes[i]);
            }

            trussModel.NodesDictionary[1].Constraints.Add(DOFType.X);
            trussModel.NodesDictionary[1].Constraints.Add(DOFType.Y);
            trussModel.NodesDictionary[2].Constraints.Add(DOFType.X);
            trussModel.NodesDictionary[2].Constraints.Add(DOFType.Y);
            
            var element1 = new Element() { ID = 1, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea} };
            var element2 = new Element() { ID = 2, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea} };

            element1.AddNode(trussModel.NodesDictionary[1]);
            element1.AddNode(trussModel.NodesDictionary[3]);

            element2.AddNode(trussModel.NodesDictionary[2]);
            element2.AddNode(trussModel.NodesDictionary[3]);

            trussModel.ElementsDictionary.Add(element1.ID, element1);
            trussModel.ElementsDictionary.Add(element2.ID, element2);

            trussModel.SubdomainsDictionary[0].ElementsDictionary.Add(element1.ID, element1);
            trussModel.SubdomainsDictionary[0].ElementsDictionary.Add(element2.ID, element2);

            trussModel.Loads.Add(new Load() { Amount = loadX, Node = trussModel.NodesDictionary[3], DOF = DOFType.X });
            trussModel.Loads.Add(new Load() { Amount = loadY, Node = trussModel.NodesDictionary[3], DOF = DOFType.Y });

            trussModel.ConnectDataStructures();


            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[0] = new SkylineLinearSystem(0, trussModel.SubdomainsDictionary[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            ProblemStructural provider = new ProblemStructural(trussModel, linearSystems);

            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] {
                trussModel.NodalDOFsDictionary[3][DOFType.X],
                trussModel.NodalDOFsDictionary[3][DOFType.Y],
                });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Console.WriteLine("Displacements of Node 3, along axes X, Y:");
            Console.WriteLine(childAnalyzer.Logs[0][0]);
        }
    }
}
