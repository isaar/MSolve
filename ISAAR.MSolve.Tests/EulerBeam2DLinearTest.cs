using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class EulerBeam2DLinearTest
    {
        [Fact]
        public void TestEulerBeam2DLinearBendingExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 2000.0;
            int nElems = 2;
            int monitorNode = 3;

            // Create new material
            ElasticMaterial material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 100.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 200.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

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
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // Create new Beam2D section and element
                var beam = new EulerBeam2D(youngModulus)
                {
                    Density = 7.85,
                    SectionArea = 91.04,
                    MomentOfInertia = 8091.00,
                };

                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };
                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                var a = beam.StiffnessMatrix(element);

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

            // Choose child analyzer -> Child: LinearAnalyzer
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(31.388982074929341, linearSystems[1].Solution[4], 12);
        }
    }
}
