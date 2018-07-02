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
    public class Beam3DQuaternionNonlinearTest
    {
        [Fact]
        public void SolveNLBeam()
        {
            var m = new Model();
            m.NodesDictionary.Add(1, new Node() { ID = 1, X = 0, Y = 0, Z = 0 });
            m.NodesDictionary.Add(2, new Node() { ID = 2, X = 5, Y = 0, Z = 0 });
            m.NodesDictionary[1].Constraints.AddRange(new[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ });
            m.ElementsDictionary.Add(1, new Element()
            {
                ID = 1,
                ElementType = new Beam3DCorotationalQuaternion(m.Nodes, new ElasticMaterial3D() { YoungModulus = 2.1e6, PoissonRatio = 0.2 }, 1,
                new BeamSection3D(0.06, 0.0002, 0.00045, 0.000818, 0.05, 0.05))
            });
            m.ElementsDictionary[1].AddNodes(m.Nodes);
            m.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            m.SubdomainsDictionary[1].ElementsDictionary.Add(1, m.ElementsDictionary[1]);
            m.Loads.Add(new Load() { Node = m.NodesDictionary[2], Amount = 100, DOF = DOFType.Y });

            m.ConnectDataStructures();
            VectorExtensions.AssignTotalAffinityCount();

            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, m.Subdomains[0].Forces);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(m.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(m.Subdomains[0]) };
            int increments = 10;
            int totalDOFs = m.TotalDOFs;
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(m, linearSystems);
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
                provider, increments, totalDOFs);

            // Choose parent analyzer -> Parent: Static
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }


        [Fact]
        public void CantileverYBeam3DQuaternionNonlinearTest()
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
            Node node1 = new Node { ID = 1, X = 0.0,   Y = 0.0, Z = 0.0 };
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

                // Add beam element to the element and subdomains dictionary of the model
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

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 10;
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

            Assert.Equal(148.936792350562, linearSystems[1].Solution[7], 2);
        }        

        [Fact]
        public void PlaneFrameTest()
        {
            VectorExtensions.AssignTotalAffinityCount();
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
            ElasticMaterial3D material = new ElasticMaterial3D
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
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain first and last nodes of the model
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[1].Constraints.Add(DOFType.Z);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);
            model.NodesDictionary[4].Constraints.Add(DOFType.X);
            model.NodesDictionary[4].Constraints.Add(DOFType.Y);
            model.NodesDictionary[4].Constraints.Add(DOFType.Z);
            model.NodesDictionary[4].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[4].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[4].Constraints.Add(DOFType.RotZ);

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
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = DOFType.X });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 10;
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

            Assert.Equal(120.1108698752, linearSystems[1].Solution[0], 2);
        }

        [Fact]
        public void SpaceFrameTest()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 21000.0;
            double poissonRatio = 0.3;
            double nodalLoad = 30000.0;
            double area = 91.04;
            double inertiaY = 2843.0; //8091.0; //
            double inertiaZ = 8091.0; //2843.0; //
            double torsionalInertia = 76.57;
            double effectiveAreaY = 91.04;
            double effectiveAreaZ = 91.04;
            int nNodes = 8;
            int nElems = 8;
            int monitorNode = 8;

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
            Node node3 = new Node { ID = 3, X = 0.0, Y = 300.0, Z = 0.0 };
            Node node4 = new Node { ID = 4, X = 300.0, Y = 300.0, Z = 0.0 };
            Node node5 = new Node { ID = 5, X = 0.0, Y = 0.0, Z = 300.0 };
            Node node6 = new Node { ID = 6, X = 300.0, Y = 0.0, Z = 300.0 };
            Node node7 = new Node { ID = 7, X = 0.0, Y = 300.0, Z = 300.0 };
            Node node8 = new Node { ID = 8, X = 300.0, Y = 300.0, Z = 300.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);
            nodes.Add(node5);
            nodes.Add(node6);
            nodes.Add(node7);
            nodes.Add(node8);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain first and last nodes of the model
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[1].Constraints.Add(DOFType.Z);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);

            model.NodesDictionary[2].Constraints.Add(DOFType.X);
            model.NodesDictionary[2].Constraints.Add(DOFType.Y);
            model.NodesDictionary[2].Constraints.Add(DOFType.Z);
            model.NodesDictionary[2].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[2].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[2].Constraints.Add(DOFType.RotZ);

            model.NodesDictionary[3].Constraints.Add(DOFType.X);
            model.NodesDictionary[3].Constraints.Add(DOFType.Y);
            model.NodesDictionary[3].Constraints.Add(DOFType.Z);
            model.NodesDictionary[3].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[3].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[3].Constraints.Add(DOFType.RotZ);

            model.NodesDictionary[4].Constraints.Add(DOFType.X);
            model.NodesDictionary[4].Constraints.Add(DOFType.Y);
            model.NodesDictionary[4].Constraints.Add(DOFType.Z);
            model.NodesDictionary[4].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[4].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[4].Constraints.Add(DOFType.RotZ);

            // Generate elements of the structure
            int iNode = 1;
            for (int iElem = 0; iElem < 4; iElem++)
            {
                // element nodes
                IList<Node> elementNodes = new List<Node>();
                elementNodes.Add(model.NodesDictionary[iNode]);
                elementNodes.Add(model.NodesDictionary[iNode + 4]);

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
                element.AddNode(model.NodesDictionary[iNode + 4]);

                // Add beam element to the element and subdomains dictionary of the model
                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);
                iNode++;
            }

            // ElementID = 5
            // element nodes
            IList<Node> elementNodes1 = new List<Node>();
            elementNodes1.Add(model.NodesDictionary[5]);
            elementNodes1.Add(model.NodesDictionary[6]);

            // Create new Beam3D section and element
            var beamSection1 = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
            var beam1 = new Beam3DCorotationalQuaternion(elementNodes1, material, 7.85, beamSection1);

            // Create elements
            var element1 = new Element()
            {
                ID = 5,
                ElementType = beam1
            };

            var a1 = beam1.StiffnessMatrix(element1);

            // Add nodes to the created element
            element1.AddNode(model.NodesDictionary[5]);
            element1.AddNode(model.NodesDictionary[6]);

            // Add beam element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element1.ID, element1);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element1.ID, element1);
            
            // ElementID = 6
            // element nodes
            IList<Node> elementNodes2 = new List<Node>();
            elementNodes2.Add(model.NodesDictionary[7]);
            elementNodes2.Add(model.NodesDictionary[8]);

            // Create new Beam3D section and element
            var beamSection2 = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
            var beam2 = new Beam3DCorotationalQuaternion(elementNodes2, material, 7.85, beamSection2);

            // Create elements
            var element2 = new Element()
            {
                ID = 6,
                ElementType = beam2
            };

            var a2 = beam2.StiffnessMatrix(element2);

            // Add nodes to the created element
            element2.AddNode(model.NodesDictionary[7]);
            element2.AddNode(model.NodesDictionary[8]);

            // Add beam element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element2.ID, element2);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element2.ID, element2);
            
            // ElementID = 7
            // element nodes
            IList<Node> elementNodes3 = new List<Node>();
            elementNodes3.Add(model.NodesDictionary[5]);
            elementNodes3.Add(model.NodesDictionary[7]);

            // Create new Beam3D section and element
            var beamSection3 = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
            var beam3 = new Beam3DCorotationalQuaternion(elementNodes3, material, 7.85, beamSection3);

            // Create elements
            var element3 = new Element()
            {
                ID = 7,
                ElementType = beam3
            };

            var a3 = beam3.StiffnessMatrix(element3);

            // Add nodes to the created element
            element3.AddNode(model.NodesDictionary[5]);
            element3.AddNode(model.NodesDictionary[7]);

            // Add beam element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element3.ID, element3);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element3.ID, element3);
            
            // ElementID = 8
            // element nodes
            IList<Node> elementNodes4 = new List<Node>();
            elementNodes4.Add(model.NodesDictionary[6]);
            elementNodes4.Add(model.NodesDictionary[8]);

            // Create new Beam3D section and element
            var beamSection4 = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
            var beam4 = new Beam3DCorotationalQuaternion(elementNodes4, material, 7.85, beamSection4);

            // Create elements
            var element4 = new Element()
            {
                ID = 8,
                ElementType = beam4
            };

            var a4 = beam4.StiffnessMatrix(element4);

            // Add nodes to the created element
            element4.AddNode(model.NodesDictionary[6]);
            element4.AddNode(model.NodesDictionary[8]);

            // Add beam element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element4.ID, element4);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element4.ID, element4);
            
            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[7], DOF = DOFType.Y });
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[8], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            int increments = 10;
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

            Assert.Equal(235.1740724707, linearSystems[1].Solution[13], 2);
        }
    }
}
