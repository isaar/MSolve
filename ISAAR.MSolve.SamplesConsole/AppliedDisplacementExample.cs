using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;

namespace ISAAR.MSolve.SamplesConsole
{
    public class AppliedDisplacementExample
    {
        private const int subdomainID = 0;

        private static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node( id: 1, x: 0.0, y:  0.0, z: 0.0 );
            Node node2 = new Node( id: 2, x: 1.0, y:  0.0, z: 0.0 );
            nodes.Add(node1);
            nodes.Add(node2);
            return nodes;
        }

        public static void Run()
        {
            double youngModulus = 200.0e06;
            double poissonRatio = 0.3;
            double nodalLoad = 25.0;

            // Create a new elastic 3D material
            var material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = CreateNodes();

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain(subdomainID));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });

            model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY, Amount = -4.16666666666667E-07 });

            //Create a new Beam2D element
            var beam = new EulerBeam2D(youngModulus)
            {
                SectionArea = 1,
                MomentOfInertia = .1
            };


            var element = new Element()
            {
                ID = 1,
                ElementType = beam
            };

            //// Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);

            //var a = beam.StiffnessMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[subdomainID].Elements.Add(element);

            // Add nodal load values at the top nodes of the model
            //model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.Y });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Output requests
            var logFactory = new TotalDisplacementsLog.Factory(model.SubdomainsDictionary[subdomainID]);
            logFactory.WatchDof(model.NodesDictionary[2], StructuralDof.TranslationX);
            logFactory.WatchDof(model.NodesDictionary[2], StructuralDof.RotationZ);
            childAnalyzer.LogFactories[subdomainID] = logFactory;

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Choose dof types X, Y, rotZ to log for node 5
            var logger = (TotalDisplacementsLog)(childAnalyzer.Logs[subdomainID][0]); //There is a list of logs for each subdomain and we want the first one
            double[] results = {
                logger.GetDisplacementAt(model.NodesDictionary[2], StructuralDof.TranslationX),
                logger.GetDisplacementAt(model.NodesDictionary[2], StructuralDof.RotationZ) };
            

            double[] expected = new double[] { 0, -4.16666666666667E-07, -6.25E-07 };

            for (int i = 0; i < expected.Length - 1; i++)
            {
                //if (Math.Abs(expected[i] - results[i]) > 1e-14)
                //{
                //    throw new SystemException("Failed beam2D test, results don't coincide for dof no: " + i + ", expected displacement: " + expected[i] + ", calculated displacement: " + results[i]);
                //}
                Console.WriteLine(results[i]);

            }
            Console.WriteLine("ran beam2d2 test");
        }
    }
}

