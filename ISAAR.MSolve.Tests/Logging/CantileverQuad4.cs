using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Tests.Logging
{
    public class CantileverQuad4
    {
        private const double length = 4.0;
        private const double height = 20.0;
        private const double thickness = 0.1;
        private const double youngModulus = 2E6;
        private const double poissonRatio = 0.3;
        private const double maxLoad = 1000.0; // TODO: this should be triangular

        [Fact]
        public static void Run()
        {
            SolveLinearStatic(CreateModelManually());
        }

        private static Model CreateModelManually()
        {
            VectorExtensions.AssignTotalAffinityCount();

            ElasticMaterial2D material = new ElasticMaterial2D()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
                StressState = "plstress"
            };

            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            #region Nodes
            Node2D[] nodes =
            {
                new Node2D(0, 0.0, 0.0),
                new Node2D(1, length, 0.0),
                new Node2D(2, 0.0, 0.25 * height),
                new Node2D(3, length, 0.25 * height),
                new Node2D(4, 0.0, 0.50 * height),
                new Node2D(5, length, 0.50 * height),
                new Node2D(6, 0.0, 0.75 * height),
                new Node2D(7, length, 0.75 * height),
                new Node2D(8, 0.0, height),
                new Node2D(9, length, height)
            };

            for (int i = 0; i < 10; ++i) model.NodesDictionary.Add(i, nodes[i]);
            #endregion

            #region Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            ContinuumElement2D[] elements =
            {
                factory.CreateQuad4(new Node2D[] { nodes[0], nodes[1], nodes[3], nodes[2]}),
                factory.CreateQuad4(new Node2D[] { nodes[2], nodes[3], nodes[5], nodes[4]}),
                factory.CreateQuad4(new Node2D[] { nodes[4], nodes[5], nodes[7], nodes[6]}),
                factory.CreateQuad4(new Node2D[] { nodes[6], nodes[7], nodes[9], nodes[8]})
            };

            for (int i = 0; i < 4; ++i)
            {
                var elementWrapper = new Element() { ID = i, ElementType = elements[i] };
                foreach (Node node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, elementWrapper);
            }
            #endregion

            #region Constraints
            var constrainedNodes = new int[] { 0, 1 };
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(DOFType.X);
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(DOFType.Y);
            }
            #endregion

            #region Loads
            model.Loads.Add(new Load() { Amount = maxLoad, Node = model.NodesDictionary[8], DOF = DOFType.X });
            #endregion

            model.ConnectDataStructures();
            return model;
        }

        private static void SolveLinearStatic(Model model)
        {
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
            string outputDirectory = @"C:\Users\Serafeim\Desktop\Presentation\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true, LogStrains = true, LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}