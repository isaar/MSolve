using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.FEM.Meshes;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.VonMisesStress;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ContinuumElementQuad4LinearTest
    {
        [Fact]
        public void TestQuad4Cantilever()
        {
            const string workingDirectory = @"C:\Users\George\Desktop";
            double youngModulus = 3.76;
            double poissonRatio = 0.3779;
            double thickness = 1.0;
            double nodalLoad = 10.0;
            int nNodes = 8;
            int nElems = 3;
            int monitorNode = 4;
            VectorExtensions.AssignTotalAffinityCount();

            // Nodal Geometry
            IReadOnlyList<Node2D> nodes = new Node2D[]
            {
                new Node2D(1, 0.0, 0.0),
                new Node2D(2, 10, 0.0),
                new Node2D(3, 20.0, 0.0),
                new Node2D(4, 30.0, 0.0),
                new Node2D(5, 0.0, 10.0),
                new Node2D(6, 10.0, 10.0),
                new Node2D(7, 20.0, 10.0),
                new Node2D(8, 30.0, 10.0)
            };

            // Create new 2D material
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Subdomains
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add Nodes to Model
            for (int i = 0; i < nNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Connectivity
            CellConnectivity2D[] elements =
            {
                new CellConnectivity2D(CellType2D.Quad4, new Node2D[] { nodes[0], nodes[1], nodes[5], nodes[4]}),
                new CellConnectivity2D(CellType2D.Quad4, new Node2D[] { nodes[1], nodes[2], nodes[6], nodes[5]}),
                new CellConnectivity2D(CellType2D.Quad4, new Node2D[] { nodes[2], nodes[3], nodes[7], nodes[6]})
            };

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            for (int i = 0; i < nElems; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, elementWrapper);
            }

            // Constraints
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[5].Constraints.Add(DOFType.X);
            model.NodesDictionary[5].Constraints.Add(DOFType.Y);

            // Loads
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = nodes[3], DOF = DOFType.X });

            // Finalize
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            // Logging displacement, strain, and stress fields.
            //string outputDirectory = workingDirectory + "\\Plots";
            //childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            //{
            //    LogDisplacements = true,
            //    LogStrains = true,
            //    LogStresses = true,
            //    VonMisesStressCalculator = new PlaneStressVonMises()
            //};

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            
            Assert.Equal(10.927956254399865, linearSystems[1].Solution[6], 2);
        }
    }
}
