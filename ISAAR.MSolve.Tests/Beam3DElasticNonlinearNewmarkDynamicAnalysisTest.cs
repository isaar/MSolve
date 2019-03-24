using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class Beam3DElasticNonlinearNewmarkDynamicAnalysisTest
    {
        //[Fact]
        private static void TestBeam3DElasticNonlinearNewmarkDynamicAnalysisExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 20000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            double density = 7.85;
            int nNodes = 3;
            int nElems = 2;
            int monitorNode = 3;

            // Create new 3D material
            ElasticMaterial3D material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 300.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 600.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            int subdomainID = 0;
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = subdomainID });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node> elementNodes = new List<Node>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new Beam3DCorotationalQuaternion(elementNodes, material, density, beamSection);

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
                var b = beam.MassMatrix(element);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.SubdomainsDictionary[subdomainID].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[subdomainID]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[subdomainID] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.SubdomainsDictionary[subdomainID]) };
            int increments = 10;
            int totalDOFs = model.TotalDOFs;
            int maximumIteration = 120;
            int iterationStepsForMatrixRebuild = 500;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static or Dynamic
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(148.936792350562, linearSystems[subdomainID].Solution[7], 12);
        }

        private static void TestBeam3DElasticNonlinearNewmarkDynamicAnalysisExample_v2()
        {
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 20000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            double density = 7.85;
            int nNodes = 3;
            int nElems = 2;
            int monitorNode = 3;

            // Create new 3D material
            var material = new ElasticMaterial3D_v2
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node_v2> nodes = new List<Node_v2>();
            Node_v2 node1 = new Node_v2 { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node_v2 node2 = new Node_v2 { ID = 2, X = 300.0, Y = 0.0, Z = 0.0 };
            Node_v2 node3 = new Node_v2 { ID = 3, X = 600.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            // Model creation
            Model_v2 model = new Model_v2();

            // Add a single subdomain to the model
            int subdomainID = 0;
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < nElems; iElem++)
            {
                // element nodes
                IList<Node_v2> elementNodes = new List<Node_v2>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 1]);

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                var beam = new Beam3DCorotationalQuaternion_v2(elementNodes, material, density, beamSection);

                // Create elements
                var element = new Element_v2()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                var a = beam.StiffnessMatrix(element);
                var b = beam.MassMatrix(element);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[subdomainID].Elements.Add(element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load_v2() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural_v2(model, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer_v2.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.MaxIterationsPerIncrement = 120;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 500;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater_v2(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer_v2 childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer_v2.Builder(model, solver, provider, childAnalyzer, 0.28, 3.36);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            NewmarkDynamicAnalyzer_v2 parentAnalyzer = parentAnalyzerBuilder.Build();

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(148.936792350562, solver.LinearSystems[subdomainID].Solution[7], 12);
        }
    }
}
