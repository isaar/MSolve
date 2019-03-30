using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;

namespace ISAAR.MSolve.Tests
{
    public class Beam2DNewmarkDynamicAanalysisTest
    {
        //[Fact]
        public void LinearElasticBeam2DNewmarkDynamicAnalysisTest()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000;
            double poissonRatio = 0.3;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double density = 7.85;
            double nodalLoad = 1000.0;
            int totalNodes = 2;

            // Newmark Dynamic parameters
            double alpha = 0.25;
            double delta = 0.50;
            double timeStep = 1.0;
            double totalTime = 100;

            // Material Definition
            ElasticMaterial material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 300.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            // Create a new Beam2D element
            var beam = new EulerBeam2D(youngModulus)
            {
                Density = density,
                SectionArea = area,
                MomentOfInertia = inertiaZ
            };

            var element = new Element()
            {
                ID = 1,
                ElementType = beam
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);

            // Element Stiffness Matrix
            var a = beam.StiffnessMatrix(element);
            var b = beam.MassMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element.ID, element);

            // define loads
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[totalNodes], DOF = DOFType.Y });
            model.ConnectDataStructures();
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            
            // Choose child analyzer -> Child: Linear Analyzer
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            
            // Choose parent analyzer -> Parent: Static or Dynamic
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, alpha, delta, timeStep, totalTime);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            Assert.Equal(2.2840249264795207, linearSystems[0].Solution[1], 8);
        }

        //[Fact]
        public void NonLinearElasticBeam2DNewmarkDynamicAnalysisTest()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000;
            double poissonRatio = 0.3;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double density = 7.85;
            double nodalLoad = 1000.0;
            int totalNodes = 2;

            // Newton-Raphson parameters
            int increments = 10;

            // Newmark Dynamic parameters
            double alpha = 0.25;
            double delta = 0.50;
            double timeStep = 1.0;
            double totalTime = 10.0;

            // Material Definition
            ElasticMaterial material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 300.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);
            
            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            // Create new Beam2D section and element
            var beamSection = new BeamSection2D(area, inertiaZ);
            var beam = new Beam2DCorotational(nodes, material, density, beamSection);

            var element = new Element()
            {
                ID = 1,
                ElementType = beam
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);

            // Element Stiffness Matrix
            var a = beam.StiffnessMatrix(element);
            var b = beam.MassMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(element.ID, element);

            // define loads
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[totalNodes], DOF = DOFType.Y });
            model.ConnectDataStructures();
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[0] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int totalDOFs = model.TotalDOFs;
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
            provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static or Dynamic
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, alpha, delta, timeStep, totalTime);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            Assert.Equal(2.2840249264795207, linearSystems[0].Solution[1], 8);
        }

        public void LinearElasticBeam2DNewmarkDynamicAnalysisTest_v2()
        {
            double youngModulus = 21000;
            double poissonRatio = 0.3;
            double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double density = 7.85;
            double nodalLoad = 1000.0;
            int totalNodes = 2;

            var material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node_v2> nodes = new List<Node_v2>();
            Node_v2 node1 = new Node_v2 { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node_v2 node2 = new Node_v2 { ID = 2, X = 300.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);

            // Model creation
            Model_v2 model = new Model_v2();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain_v2(1));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });

            // Create a new Beam2D element
            var beam = new EulerBeam2D_v2(youngModulus)
            {
                Density = density,
                SectionArea = area,
                MomentOfInertia = inertiaZ
            };

            var element = new Element_v2()
            {
                ID = 1,
                ElementType = beam
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);

            // Element Stiffness Matrix
            var a = beam.StiffnessMatrix(element);
            var b = beam.MassMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[1].Elements.Add(element);

            // define loads
            model.Loads.Add(new Load_v2 { Amount = nodalLoad, Node = model.NodesDictionary[totalNodes], DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural_v2(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer_v2.Builder(model, solver, provider, childAnalyzer, 0.28, 3.36);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            NewmarkDynamicAnalyzer_v2 parentAnalyzer = parentAnalyzerBuilder.Build();
     
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(2.2840249264795207, solver.LinearSystems[1].Solution[1], 8);
        }
    }
}
