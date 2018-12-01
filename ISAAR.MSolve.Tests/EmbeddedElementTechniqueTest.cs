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
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.FEM.Embedding;
using System.Linq;
using ISAAR.MSolve.Discretization.Integration.Quadratures;

namespace ISAAR.MSolve.Tests
{
    public class EmbeddedElementTechniqueTest
    {
        [Fact]
        public void EmbeddedElementTechniqueExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Choose model
            EmbeddedExamplesBuilder.ExampleWithEmbedded(model);
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);

            // Skyline Solver
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int totalDOFs = model.TotalDOFs;
            int increments = 10;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[5][DOFType.X],
            model.NodalDOFsDictionary[5][DOFType.Y],
            model.NodalDOFsDictionary[5][DOFType.Z]});

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            var expectedValue = linearSystems[1].Solution[11];
            Assert.Equal(11.584726466617692, expectedValue, 3);
        }

        public static class EmbeddedExamplesBuilder
        {
            public static void HostElementsBuilder(Model model)
            {
                // Nodes Geometry
                model.NodesDictionary.Add(1, new Node() { ID = 1, X = 10.00, Y = 2.50, Z = 2.50 });
                model.NodesDictionary.Add(2, new Node() { ID = 2, X = 0.00, Y = 2.50, Z = 2.50 });
                model.NodesDictionary.Add(3, new Node() { ID = 3, X = 0.00, Y = -2.50, Z = 2.50 });
                model.NodesDictionary.Add(4, new Node() { ID = 4, X = 10.00, Y = -2.50, Z = 2.50 });
                model.NodesDictionary.Add(5, new Node() { ID = 5, X = 10.00, Y = 2.50, Z = -2.50 });
                model.NodesDictionary.Add(6, new Node() { ID = 6, X = 0.00, Y = 2.50, Z = -2.50 });
                model.NodesDictionary.Add(7, new Node() { ID = 7, X = 0.00, Y = -2.50, Z = -2.50 });
                model.NodesDictionary.Add(8, new Node() { ID = 8, X = 10.00, Y = -2.50, Z = -2.50 });

                // Boundary Conditions
                model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = DOFType.Z });
                model.NodesDictionary[3].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[3].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[3].Constraints.Add(new Constraint { DOF = DOFType.Z });
                model.NodesDictionary[6].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[6].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[6].Constraints.Add(new Constraint { DOF = DOFType.Z });
                model.NodesDictionary[7].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[7].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[7].Constraints.Add(new Constraint { DOF = DOFType.Z });

                // Create Material
                ElasticMaterial3D solidMaterial = new ElasticMaterial3D()
                {
                    YoungModulus = 3.76,
                    PoissonRatio = 0.3779,
                };

                // Hexa8NL element definition
                Element hexa8NLelement = new Element()
                {
                    ID = 1,
                    ElementType = new Hexa8NonLinear(solidMaterial, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
                };

                // Add nodes to the created element
                hexa8NLelement.AddNode(model.NodesDictionary[1]);
                hexa8NLelement.AddNode(model.NodesDictionary[2]);
                hexa8NLelement.AddNode(model.NodesDictionary[3]);
                hexa8NLelement.AddNode(model.NodesDictionary[4]);
                hexa8NLelement.AddNode(model.NodesDictionary[5]);
                hexa8NLelement.AddNode(model.NodesDictionary[6]);
                hexa8NLelement.AddNode(model.NodesDictionary[7]);
                hexa8NLelement.AddNode(model.NodesDictionary[8]);

                // Add Hexa element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);

                // Add nodal load values at the top nodes of the model
                model.Loads.Add(new Load() { Amount = 25, Node = model.NodesDictionary[1], DOF = DOFType.Z });
                model.Loads.Add(new Load() { Amount = 25, Node = model.NodesDictionary[4], DOF = DOFType.Z });
                model.Loads.Add(new Load() { Amount = 25, Node = model.NodesDictionary[5], DOF = DOFType.Z });
                model.Loads.Add(new Load() { Amount = 25, Node = model.NodesDictionary[8], DOF = DOFType.Z });
            }

            public static void EmbeddedElementsBuilder(Model model)
            {
                // define mechanical properties
                double youngModulus = 1.0;
                double shearModulus = 1.0;
                double poissonRatio = (youngModulus / (2 * shearModulus)) - 1;
                double area = 1776.65;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
                double inertiaY = 1058.55;
                double inertiaZ = 1058.55;
                double torsionalInertia = 496.38;
                double effectiveAreaY = area;
                double effectiveAreaZ = area;

                // Geometry
                model.NodesDictionary.Add(9, new Node() { ID = 9, X = 0.00, Y = 0.00, Z = 0.00 });
                model.NodesDictionary.Add(10, new Node() { ID = 10, X = 10.00, Y = 0.00, Z = 0.00 });

                // Create new 3D material
                ElasticMaterial3D beamMaterial = new ElasticMaterial3D
                {
                    YoungModulus = youngModulus,
                    PoissonRatio = poissonRatio,
                };

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                // element nodes
                IList<Node> elementNodes = new List<Node>();
                elementNodes.Add(model.NodesDictionary[9]);
                elementNodes.Add(model.NodesDictionary[10]);
                var beam = new Beam3DCorotationalQuaternion(elementNodes, beamMaterial, 7.85, beamSection);
                var beamElement = new Element { ID = 2, ElementType = beam };

                beamElement.NodesDictionary.Add(9, model.NodesDictionary[9]);
                beamElement.NodesDictionary.Add(10, model.NodesDictionary[10]);

                model.ElementsDictionary.Add(beamElement.ID, beamElement);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(beamElement.ID, beamElement);
            }

            public static void ExampleWithEmbedded(Model model)
            {
                HostElementsBuilder(model);
                EmbeddedElementsBuilder(model);
                var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key == 1).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key == 2).Select(kv => kv.Value), true);
            }
        }
    }
}
