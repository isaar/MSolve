using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;

namespace ISAAR.MSolve.SamplesConsole.Solvers
{
    public static class FloatingSubdomainMatrices
    {
        private const string outputDirectory = @"C:\Users\Serafeim\Desktop\FETI";
        private const int subdomainID = 0;

        public static void WriteStiffnessOfContinuum2DStructure()
        {
            //   ____ ____
            //  |    |    |
            //  |____|____|
            //  |    |    |
            //  |____|____|

            // Model with 1 subdomain
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

            // Material
            double thickness = 1.0;
            var material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = 2.1E7,
                PoissonRatio = 0.3
            };
            var dynamicProperties = new DynamicMaterial(1.0, 0.0, 0.0);

            // Generate mesh
            double domainLength = 2.0;
            double domainHeight = 2.4;
            var meshGenerator = new UniformMeshGenerator2D(0.0, 0.0, domainLength, domainHeight, 10, 10);
            (IReadOnlyList<Node> vertices, IReadOnlyList<CellConnectivity> cells) = meshGenerator.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Add Quad4 elements to the model
            var factory = new ContinuumElement2DFactory(thickness, material, dynamicProperties);
            for (int e = 0; e < cells.Count; ++e)
            {
                ContinuumElement2D element = factory.CreateElement(cells[e].CellType, cells[e].Vertices);
                var elementWrapper = new Element() { ID = e, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(e, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            SkylineSolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis to build the stiffness matrix
            parentAnalyzer.Initialize();
            parentAnalyzer.BuildMatrices();

            // Print the stiffness matrix
            //var writer = new MatlabWriter();
            var writer = new RawArraysWriter();
            writer.WriteToMultipleFiles((SkylineMatrix)solver.LinearSystems[subdomainID].Matrix, 
                outputDirectory + @"\quad4_20x20_stiffness.txt");
        }

        public static void WriteStiffnessOfContinuum3DStructure()
        {
            //    _____ ____
            //   /    /    /|
            //  /____/____/ |
            //  |    |    | |
            //  |____|____|/|
            //  |    |    | /
            //  |____|____|/
            //  

            // Model with 1 subdomain
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

            // Material
            var material = new ElasticMaterial3D()
            {
                YoungModulus = 2.1E7,
                PoissonRatio = 0.3
            };
            var dynamicProperties = new DynamicMaterial(1.0, 0.0, 0.0);

            // Generate mesh
            double lengthX = 2.0;
            double lengthY = 2.4;
            double lengthZ = 2.2;
            var meshGenerator = new UniformMeshGenerator3D(0.0, 0.0, 0.0, lengthX, lengthY, lengthZ, 10, 10, 10);
            (IReadOnlyList<Node> vertices, IReadOnlyList<CellConnectivity> cells) = meshGenerator.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Add Quad4 elements to the model
            var factory = new ContinuumElement3DFactory(material, dynamicProperties);
            for (int e = 0; e < cells.Count; ++e)
            {
                ContinuumElement3D element = factory.CreateElement(cells[e].CellType, cells[e].Vertices);
                var elementWrapper = new Element() { ID = e, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(e, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            SkylineSolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis to build the stiffness matrix
            parentAnalyzer.Initialize();
            parentAnalyzer.BuildMatrices();

            // Print the stiffness matrix
            //var writer = new MatlabWriter();
            var writer = new RawArraysWriter();
            writer.WriteToMultipleFiles((SkylineMatrix)solver.LinearSystems[subdomainID].Matrix,
                outputDirectory + @"\hexa8_10x10x10_stiffness.txt");
        }
    }
}
