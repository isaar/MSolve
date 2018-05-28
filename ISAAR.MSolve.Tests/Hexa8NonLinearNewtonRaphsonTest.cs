using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class Hexa8NonLinearNewtonRaphsonTest
    {
        [Fact]
        public void TestHexa8NonLinearNewtonRaphsonExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 20000.0;
            int nNodes = 12;
            int nElems = 5;
            int monitorNode = 12;

            // Create new 3D material
            ElasticMaterial3D_v2 material = new ElasticMaterial3D_v2
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 100.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 200.0, Y = 0.0, Z = 0.0 };
            Node node4 = new Node { ID = 4, X = 300.0, Y = 0.0, Z = 0.0 };
            Node node5 = new Node { ID = 5, X = 400.0, Y = 0.0, Z = 0.0 };
            Node node6 = new Node { ID = 6, X = 500.0, Y = 0.0, Z = 0.0 };
            Node node7 = new Node { ID = 7, X = 0.0, Y = 100.0, Z = 0.0 };
            Node node8 = new Node { ID = 8, X = 100.0, Y = 100.0, Z = 0.0 };
            Node node9 = new Node { ID = 9, X = 200.0, Y = 100.0, Z = 0.0 };
            Node node10 = new Node { ID = 10, X = 300.0, Y = 100.0, Z = 0.0 };
            Node node11 = new Node { ID = 11, X = 400.0, Y = 100.0, Z = 0.0 };
            Node node12 = new Node { ID = 12, X = 500.0, Y = 100.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);
            nodes.Add(node5);
            nodes.Add(node6);
            nodes.Add(node7);
            nodes.Add(node8);
            nodes.Add(node9);
            nodes.Add(node10);
            nodes.Add(node11);
            nodes.Add(node12);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[1].Constraints.Add(DOFType.Z);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);
            // Generate elements of the structure
            int iNode = 0;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // Create new Beam2D section and element
                //var beam = new EulerBeam2D(youngModulus)
                //{
                //    Density = 7.85,
                //    SectionArea = 91.04,
                //    MomentOfInertia = 8091.00,
                //};

                // Create new Beam3D section and element
                var hexa = new Hexa8NL_v2();

                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = hexa
                };
                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                var a = hexa.StiffnessMatrix(element);

                // Add Hexa element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 10;
            int totalDOFs = model.TotalDOFs;
            int maximumIteration = 120;
            int iterationStepsForMatrixRebuild = 500;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(148.93590375922295, linearSystems[1].Solution[7], 12);
        }
    }
}
