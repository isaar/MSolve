using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
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
            double youngModulus = 3.76;
            double poissonRatio = 0.3779;
            double nodalLoad = 500.0;
            int nNodes = 8;
            int nElems = 1;
            int monitorNode = 8;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 10.0, Z = 10.0 };
            Node node2 = new Node { ID = 2, X = 0.0, Y = 0.0, Z = 10.0 };
            Node node3 = new Node { ID = 3, X = 10.0, Y = 0.0, Z = 10.0 };
            Node node4 = new Node { ID = 4, X = 10.0, Y = 10.0, Z = 10.0 };
            Node node5 = new Node { ID = 5, X = 0.0, Y = 10.0, Z = 0.0 };
            Node node6 = new Node { ID = 6, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node7 = new Node { ID = 7, X = 10.0, Y = 0.0, Z = 0.0 };
            Node node8 = new Node { ID = 8, X = 10.0, Y = 10.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);
            nodes.Add(node5);
            nodes.Add(node6);
            nodes.Add(node7);
            nodes.Add(node8);

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

            model.NodesDictionary[2].Constraints.Add(DOFType.X);
            model.NodesDictionary[2].Constraints.Add(DOFType.Y);
            model.NodesDictionary[2].Constraints.Add(DOFType.Z);

            model.NodesDictionary[5].Constraints.Add(DOFType.X);
            model.NodesDictionary[5].Constraints.Add(DOFType.Y);
            model.NodesDictionary[5].Constraints.Add(DOFType.Z);

            model.NodesDictionary[6].Constraints.Add(DOFType.X);
            model.NodesDictionary[6].Constraints.Add(DOFType.Y);
            model.NodesDictionary[6].Constraints.Add(DOFType.Z);
            
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
                var hexa = new Hexa8NL_v2(material, 2, 2, 2);

                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = hexa
                };
                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[1]);
                element.AddNode(model.NodesDictionary[2]);
                element.AddNode(model.NodesDictionary[3]);
                element.AddNode(model.NodesDictionary[4]);
                element.AddNode(model.NodesDictionary[5]);
                element.AddNode(model.NodesDictionary[6]);
                element.AddNode(model.NodesDictionary[7]);
                element.AddNode(model.NodesDictionary[8]);

                var a = hexa.StiffnessMatrix(element);

                // Add Hexa element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[3], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[4], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[7], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[8], DOF = DOFType.X });

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
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(10.927956254399865, linearSystems[1].Solution[6], 2);
        }
    }
}
