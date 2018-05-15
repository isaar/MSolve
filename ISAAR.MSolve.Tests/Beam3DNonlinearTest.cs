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
        private readonly INonLinearSubdomainUpdater[] subdomainUpdaters;
        private readonly ISubdomainGlobalMapping[] mappings;

        [Fact]
        public void TestBeam3DNonlinearExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 2.0e08;
            double poissonRatio = 0.3;
            double nodalLoad = 10.0;
            double area = 1.0;
            double inertiaY = 1.0;
            double inertiaZ = 1.0;
            double torsionalInertia = 1.0;
            double effectiveAreaY = 1.0;
            double effectiveAreaZ = 1.0;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0,   Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 500.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);

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

            // Create new Beam3D section
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

            // Create a new Beam3D element
            var beam = new Beam3DCorotationalQuaternion(nodes, material, 1.0, beamSection);

            var element = new Element()
            {
                ID = 1,
                ElementType = beam
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);

            var a = beam.StiffnessMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);

            //int increments = 10;
            //int totalDOFs = model.TotalDOFs;
            //int maximumIteration = 120;
            //int iterationStepsForMatrixRebuild = 500;
            //NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems[1], subdomainUpdaters, mappings,
            //provider, increments, totalDOFs, maximumIteration, iterationStepsForMatrixRebuild);

            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(-2.08333333333333333e-5, linearSystems[1].Solution[1], 10);
        }
    }
}
