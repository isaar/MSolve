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
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class Beam3DNonlinearTest
    {
        [Fact]
        public void TestBeam3DNonlinearExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 5000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            int nNodes = 6;
            int nElems = 5;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0,   Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 100.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 200.0, Y = 0.0, Z = 0.0 };
            Node node4 = new Node { ID = 4, X = 300.0, Y = 0.0, Z = 0.0 };
            Node node5 = new Node { ID = 5, X = 400.0, Y = 0.0, Z = 0.0 };
            Node node6 = new Node { ID = 6, X = 500.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);
            nodes.Add(node5);
            nodes.Add(node6);

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
                // Create new Beam3D section
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

                // Create a new Beam3D element
                var beam = new Beam3DCorotationalQuaternion(nodes, material, 1.0, beamSection);

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
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[6], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            //LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 50;
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

            Assert.Equal(49.1812791359803, linearSystems[1].Solution[31-6], 12);
        }
    }
}
