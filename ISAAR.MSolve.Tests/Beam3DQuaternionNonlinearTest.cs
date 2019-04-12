using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class Beam3DQuaternionNonlinearTest
    {
        private const string output = @"C:\Users\Serafeim\Desktop\output.txt";

        public void SolveNLBeam()
        {
            var m = new Model();
            m.NodesDictionary.Add(1, new Node() { ID = 1, X = 0, Y = 0, Z = 0 });
            m.NodesDictionary.Add(2, new Node() { ID = 2, X = 5, Y = 0, Z = 0 });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            m.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            m.ElementsDictionary.Add(1, new Element()
            {
                ID = 1,
                ElementType = new Beam3DCorotationalQuaternion(m.Nodes, new ElasticMaterial3D() { YoungModulus = 2.1e6, PoissonRatio = 0.2 }, 1,
                new BeamSection3D(0.06, 0.0002, 0.00045, 0.000818, 0.05, 0.05))
            });
            m.ElementsDictionary[1].AddNodes(m.Nodes);
            m.SubdomainsDictionary.Add(1, new Subdomain(1));
            m.SubdomainsDictionary[1].Elements.Add(m.ElementsDictionary[1]);
            m.Loads.Add(new Load() { Node = m.NodesDictionary[2], Amount = 100, DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(m);

            // Problem type
            var provider = new ProblemStructural(m, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(m, solver, provider, increments);
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(m, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void CantileverYBeam3DQuaternionNonlinearTest()
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
            int nNodes = 3;
            int nElems = 2;
            int monitorNode = 3;

            // Create new 3D material
            var material = new ElasticMaterial3D
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
                var beam = new Beam3DCorotationalQuaternion(elementNodes, material, 7.85, beamSection);

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
                //var writer = new FullMatrixWriter();
                //writer.WriteToFile(a, output);


                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].Elements.Add(element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-3;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 7 });

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
            Assert.Equal(148.936792350562, log.DOFValues[7], 2);
        }

        [Fact]
        public void PlaneFrameTest()
        {
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 500000.0;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            int nNodes = 4;
            int nElems = 3;
            int monitorNode = 2;

            // Create new 3D material
            var material = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 0.0, Y = 100.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 100.0, Y = 100.0, Z = 0.0 };
            Node node4 = new Node { ID = 4, X = 100.0, Y = 0.0, Z = 0.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain(1));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain first and last nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.Z });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotX });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotY });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

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
                var beam = new Beam3DCorotationalQuaternion(elementNodes, material, 7.85, beamSection);

                // Create elements
                var element = new Element()
                {
                    ID = iElem + 1,
                    ElementType = beam
                };

                var a = beam.StiffnessMatrix(element);

                // Add nodes to the created element
                element.AddNode(model.NodesDictionary[iNode]);
                element.AddNode(model.NodesDictionary[iNode + 1]);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].Elements.Add(element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.X });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-3;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 0 });

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            // Check output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
            Assert.Equal(120.1108698752, log.DOFValues[0], 2);
        }
    }
}
