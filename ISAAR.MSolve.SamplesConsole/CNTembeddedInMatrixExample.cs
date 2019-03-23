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
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.FEM.Embedding;
using System.Linq;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using System.IO;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Analyzers.NonLinear;

namespace ISAAR.MSolve.SamplesConsole
{
    public class CNTembeddedInMatrixExample
    {
        private const string outputDirectory = @"E:\GEORGE_DATA\DESKTOP\MSolveResults";
        private const int subdomainID = 0;

        public static void EmbeddedCNTinMatrix_NewtonRaphson()
        {
            VectorExtensions.AssignTotalAffinityCount();

            // No. of increments
            int increments = 1000;

            // Model creation
            var model = new Model_v2();

            // Subdomains
            //model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = 1 });
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Variables
            int monitorNode = 10001;
            DOFType monitorDof = DOFType.Z;

            // Choose model
            EmbeddedModelBuilder.EmbeddedExample(model);

            // Boundary Conditions - Left End [End-1]
            for (int iNode = 1; iNode <= 100; iNode++)
            {
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Z });
            }

            // Boundary Conditions - Bottom End [End-3]
            for (int iNode = 1; iNode <= 10001; iNode += 100)
            {
                for (int j = 0; j < 10; j++)
                {
                    model.NodesDictionary[iNode + j].Constraints.Add(new Constraint { DOF = DOFType.Y });
                }
            }

            // Boundary Conditions - Right End [End-5]
            for (int iNode = 1; iNode <= 10091; iNode += 10)
            {
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.X });
            }

            //// Boundary Conditions - Left End [End-6]
            //for (int iNode = 10; iNode <= 10100; iNode += 10)
            //{
            //    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.X });
            //}

            //Compression Loading
            double nodalLoad = -1.0; //0.40;
            for (int iNode = 10001; iNode <= 10100; iNode++) //[End-4]
            {
                model.Loads.Add(new Load_v2() { Amount = nodalLoad, Node = model.NodesDictionary[iNode], DOF = DOFType.Z });
            }

            // Choose linear equation system solver
            //var solverBuilder = new SkylineSolver.Builder();
            //SkylineSolver solver = solverBuilder.BuildSolver(model);
            var solverBuilder = new SuiteSparseSolver.Builder();
            SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural_v2(model, solver);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer            
            var childAnalyzerBuilder = new LoadControlAnalyzer_v2.Builder(model, solver, provider, increments) { ResidualTolerance = 1E-03 };
            LoadControlAnalyzer_v2 childAnalyzer = childAnalyzerBuilder.Build();

            // Choose parent analyzer -> Parent: Static
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Request output
            string outputFile = outputDirectory + "\\CNT-Embedded-3D_Results.txt";
            var logger = new TotalLoadsDisplacementsPerIncrementLog(model.SubdomainsDictionary[subdomainID], increments,
                model.NodesDictionary[monitorNode], monitorDof, outputFile);
            childAnalyzer.IncrementalLogs.Add(subdomainID, logger);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        public static void EmbeddedCNTinMatrix_DisplacementControl()
        {
            VectorExtensions.AssignTotalAffinityCount();

            // No. of increments
            int increments = 10;

            // Model creation
            var model = new Model_v2();

            // Subdomains
            //model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = 1 });
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Variables
            int monitorNode = 10001;
            DOFType monitorDof = DOFType.Z;

            // Choose model
            EmbeddedModelBuilder.EmbeddedExample(model);

            // Boundary Conditions - Left End [End-1]
            for (int iNode = 1; iNode <= 100; iNode++)
            {
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Z });
            }

            // Boundary Conditions - Bottom End [End-3]
            for (int iNode = 1; iNode <= 10001; iNode += 100)
            {
                for (int j = 0; j < 10; j++)
                {
                    model.NodesDictionary[iNode + j].Constraints.Add(new Constraint { DOF = DOFType.Y });
                }
            }

            // Applied Displacements [End-4]
            double nodalDisplacement = -2.0;
            for (int iNode = 10001; iNode <= 10100; iNode++)
            {
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Z, Amount = nodalDisplacement });
            }

            // Choose linear equation system solver
            //var solverBuilder = new SkylineSolver.Builder();
            //SkylineSolver solver = solverBuilder.BuildSolver(model);
            var solverBuilder = new SuiteSparseSolver.Builder();
            SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural_v2(model, solver);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var childAnalyzerBuilder = new DisplacementControlAnalyzer_v2.Builder(model, solver, provider, increments);
            DisplacementControlAnalyzer_v2 childAnalyzer = childAnalyzerBuilder.Build();

            // Choose parent analyzer -> Parent: Static
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Request output
            string outputFile = outputDirectory + "\\CNT-Embedded-3D_Results.txt";
            var logger = new TotalLoadsDisplacementsPerIncrementLog(model.SubdomainsDictionary[subdomainID], increments,
                model.NodesDictionary[monitorNode], monitorDof, outputFile);
            childAnalyzer.IncrementalLogs.Add(subdomainID, logger);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        public static class EmbeddedModelBuilder
        {
            public static void EmbeddedExample(Model_v2 model)
            {
                HostElementsBuilder(model);
                EmbeddedElementsBuilder(model);
                //var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key <= 2500000).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 2500000).Select(kv => kv.Value), true);
                //var embeddedGrouping = new EmbeddedGrouping_v2(model, model.ElementsDictionary.Where(x => x.Key <= 2500).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 2500).Select(kv => kv.Value), true);
                //var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key <= 312500).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 312500).Select(kv => kv.Value), true);
                //var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key <= 1250).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 1250).Select(kv => kv.Value), true);
                //var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key <= 90).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 90).Select(kv => kv.Value), true);
                //var embeddedGrouping = new EmbeddedGrouping_v2(model, model.ElementsDictionary.Where(x => x.Key <= 90).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 90).Select(kv => kv.Value), true);
                var embeddedGrouping = new EmbeddedGrouping_v2(model, model.ElementsDictionary.Where(x => x.Key <= 8100).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 8100).Select(kv => kv.Value), true);
            }

            public static void HostElementsBuilder(Model_v2 model)
            {
                string workingDirectory = @"E:\GEORGE_DATA\DESKTOP\input files"; //"..\..\..\Resources\Beam3DInputFiles";

                string MatrixGeometryFileName = "MATRIX_3D-L_x=10-L_y=10-L_z=100-9x9x100-Geometry_MSolve.inp";
                //"MATRIX_3D-L_x=5-L_y=5-L_z=100-50x50x1000-Geometry_MSolve.inp"; 
                //"MATRIX_3D-L_x=5-L_y=5-L_z=100-5x5x100-Geometry_MSolve.inp"; 
                //"MATRIX_3D-L_x=5-L_y=5-L_z=100-25x25x500-Geometry_MSolve.inp";
                //"MATRIX_3D-L_x=1-L_y=1-L_z=10-5x5x50-Geometry_MSolve.inp";
                //"MATRIX_3D-L_x=30-L_y=30-L_z=100-3x3x10-Geometry_MSolve.inp";
                //"MATRIX_3D-L_x=30-L_y=30-L_z=100-29x29x100-Geometry_MSolve.inp";
                //"MATRIX_3D-L_x=10-L_y=10-L_z=100-19x19x200-Geometry_MSolve.inp";
                //"MATRIX_3D-L_x=10-L_y=10-L_z=100-9x9x100-Geometry_MSolve.inp";

                string MatrixGonnectivityFileName = "MATRIX_3D-L_x=10-L_y=10-L_z=100-9x9x100-ConnMatr_MSolve.inp";
                //"MATRIX_3D-L_x=5-L_y=5-L_z=100-50x50x1000-ConnMatr_MSolve.inp"; 
                //"MATRIX_3D-L_x=5-L_y=5-L_z=100-5x5x100-ConnMatr_MSolve.inp"; 
                //"MATRIX_3D-L_x=5-L_y=5-L_z=100-25x25x500-ConnMatr_MSolve.inp";
                //"MATRIX_3D-L_x=1-L_y=1-L_z=10-5x5x50-ConnMatr_MSolve.inp";
                //"MATRIX_3D-L_x=30-L_y=30-L_z=100-3x3x10-ConnMatr_MSolve.inp";
                //"MATRIX_3D-L_x=30-L_y=30-L_z=100-29x29x100-ConnMatr_MSolve.inp";
                //"MATRIX_3D-L_x=10-L_y=10-L_z=100-19x19x200-ConnMatr_MSolve.inp";
                //"MATRIX_3D-L_x=10-L_y=10-L_z=100-9x9x100-ConnMatr_MSolve.inp";

                int matrixNodes = File.ReadLines(workingDirectory + '\\' + MatrixGeometryFileName).Count();
                int matrixElements = File.ReadLines(workingDirectory + '\\' + MatrixGonnectivityFileName).Count();
                
                // Nodes Geometry                
                using (TextReader reader = File.OpenText(workingDirectory + '\\' + MatrixGeometryFileName))
                {
                    for (int i = 0; i < matrixNodes; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int nodeID = int.Parse(bits[0]);
                        double nodeX = double.Parse(bits[1]);
                        double nodeY = double.Parse(bits[2]);
                        double nodeZ = double.Parse(bits[3]);
                        model.NodesDictionary.Add(nodeID, new Node_v2 { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    }
                }                

                // Create Material
                var solidMaterial = new ElasticMaterial3D_v2()
                {
                    YoungModulus = 1.00,
                    PoissonRatio = 0.30,
                };

                // Generate elements
                using (TextReader reader = File.OpenText(workingDirectory + '\\' + MatrixGonnectivityFileName))
                {
                    for (int i = 0; i < matrixElements; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int elementID = int.Parse(bits[0]);
                        int node1 = int.Parse(bits[1]);
                        int node2 = int.Parse(bits[2]);
                        int node3 = int.Parse(bits[3]);
                        int node4 = int.Parse(bits[4]);
                        int node5 = int.Parse(bits[5]);
                        int node6 = int.Parse(bits[6]);
                        int node7 = int.Parse(bits[7]);
                        int node8 = int.Parse(bits[8]);
                        // Hexa8NL element definition
                        var hexa8NLelement = new Element_v2()
                        {
                            ID = elementID,
                            ElementType = new Hexa8NonLinear_v2(solidMaterial, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
                        };
                        // Add nodes to the created element
                        hexa8NLelement.AddNode(model.NodesDictionary[node1]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node2]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node3]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node4]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node5]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node6]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node7]);
                        hexa8NLelement.AddNode(model.NodesDictionary[node8]);
                        // Add Hexa element to the element and subdomains dictionary of the model
                        model.ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);
                        //model.SubdomainsDictionary[0].ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);
                        model.SubdomainsDictionary[0].Elements.Add(hexa8NLelement);
                    }                    
                }                
            }

            public static void EmbeddedElementsBuilder(Model_v2 model)
            {
                // define mechanical properties
                double youngModulus = 16710.0; // 5490; // 
                double shearModulus = 8080.0; // 871; // 
                double poissonRatio = 2.15; // 0.034; //(youngModulus / (2 * shearModulus)) - 1;
                double area = 5.594673861218848d - 003;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
                double inertiaY = 2.490804749753243D - 006; //1058.55;
                double inertiaZ = 2.490804749753243D - 006; // 1058.55;
                double torsionalInertia = 4.981609499506486D - 006; //496.38;
                double effectiveAreaY = area;
                double effectiveAreaZ = area;
                string workingDirectory = @"E:\GEORGE_DATA\DESKTOP\input files"; //"..\..\..\Resources\Beam3DInputFiles";

                string CNTgeometryFileName = "CNT-8-8-L=100-h=3-Geometry.inp";
                //"CNT-8-8-L=100-Geometry-2.inp"; 
                //"EmbeddedCNT-8-8-L=100-h=0-k=1-EBE-L=1-NumberOfCNTs=1-Geometry_beam.inp"; 
                //"CNT-8-8-L=100-h=3-Geometry.inp";

                string CNTconnectivityFileName = "CNT-8-8-L=100-h=3-ConnMatr.inp";
                //"CNT-8-8-L=100-ConnMatr-2.inp"; 
                //"EmbeddedCNT-8-8-L=100-h=0-k=1-EBE-L=1-NumberOfCNTs=1-ConnMatr_beam.inp"; 
                //"CNT-8-8-L=100-h=3-ConnMatr.inp";

                int CNTNodes = File.ReadLines(workingDirectory + '\\' + CNTgeometryFileName).Count();
                int CNTElems = File.ReadLines(workingDirectory + '\\' + CNTconnectivityFileName).Count();

                // Geometry
                using (TextReader reader = File.OpenText(workingDirectory + '\\' + CNTgeometryFileName))
                {
                    for (int i = 0; i < CNTNodes; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int nodeID = int.Parse(bits[0]) + 10100; // matrixNodes
                        double nodeX = double.Parse(bits[1]);
                        double nodeY = double.Parse(bits[2]);
                        double nodeZ = double.Parse(bits[3]);
                        model.NodesDictionary.Add(nodeID, new Node_v2 { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    }
                }

                // Create new 3D material
                var beamMaterial = new ElasticMaterial3D_v2
                {
                    YoungModulus = youngModulus,
                    PoissonRatio = poissonRatio,
                };

                // Create new Beam3D section and element
                var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
                // element nodes

                using (TextReader reader = File.OpenText(workingDirectory + '\\' + CNTconnectivityFileName))
                {
                    for (int i = 0; i < CNTElems; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int elementID = int.Parse(bits[0]) + 8100; // matrixElements
                        int node1 = int.Parse(bits[1]) + 10100; // matrixNodes
                        int node2 = int.Parse(bits[2]) + 10100; // matrixNodes
                        // element nodes
                        var elementNodes = new List<Node_v2>();
                        elementNodes.Add(model.NodesDictionary[node1]);
                        elementNodes.Add(model.NodesDictionary[node2]);
                        // create element
                        var beam_1 = new Beam3DCorotationalQuaternion_v2(elementNodes, beamMaterial, 7.85, beamSection);
                        var beamElement = new Element_v2 { ID = elementID, ElementType = beam_1 };
                        // Add nodes to the created element
                        beamElement.AddNode(model.NodesDictionary[node1]);
                        beamElement.AddNode(model.NodesDictionary[node2]);
                        // beam stiffness matrix
                        // var a = beam_1.StiffnessMatrix(beamElement);
                        // Add beam element to the element and subdomains dictionary of the model
                        model.ElementsDictionary.Add(beamElement.ID, beamElement);
                        //model.SubdomainsDictionary[0].ElementsDictionary.Add(beamElement.ID, beamElement);
                        model.SubdomainsDictionary[0].Elements.Add(beamElement);
                    }                    
                }
            }            
        }
    }
}
