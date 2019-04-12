using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class EulerBeam2DLinearTest
    {
        [Fact]
        public void TestEulerBeam2DLinearBendingExample()
        {
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 2000.0;
            int nElems = 2;
            int monitorNode = 3;

            // Create new material
            var material = new ElasticMaterial()
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
            model.SubdomainsDictionary.Add(1, new Subdomain(1));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
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
                model.SubdomainsDictionary[1].Elements.Add(element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = StructuralDof.TranslationY });

            // Solver
            //var solverBuilder = new MSolve.Solvers.Iterative.PcgSolver.Builder(
            //    (new LinearAlgebra.Iterative.ConjugateGradient.PcgAlgorithm.Builder()).Build());
            //var solver = solverBuilder.BuildSolver(model);
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the anlaysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            double solutionNorm = solver.LinearSystems[1].Solution.Norm2();
            double rhsNorm = solver.LinearSystems[1].RhsVector.Norm2();

            Assert.Equal(31.388982074929341, solver.LinearSystems[1].Solution[4], 12);
        }
    }
}
