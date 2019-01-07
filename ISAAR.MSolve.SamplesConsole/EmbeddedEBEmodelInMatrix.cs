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

namespace ISAAR.MSolve.SamplesConsole
{
    public class EmbeddedEBEmodelInMatrix
    {
        private const string outputDirectory = @"E:\GEORGE_DATA\DESKTOP\MSolveResults";
        private const int subdomainID = 0;

        public static void EmbeddedEBEinMatrix_NewtonRaphson()
        {
            VectorExtensions.AssignTotalAffinityCount();

            // Model creation
            var model = new Model_v2();

            // Subdomains
            //model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = 1 });
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Variables
            int monitorNode = 161;
            DOFType monitorDof = DOFType.Y;

            // Choose model
            EmbeddedEBEModelBuilder.EmbeddedExampleBuilder(model);

            //-----------------------------------------------------------------------------------------------------------
            // Model_v2

            // Choose linear equation system solver
            //var solverBuilder = new SkylineSolver.Builder();
            //SkylineSolver solver = solverBuilder.BuildSolver(model);
            var solverBuilder = new SuiteSparseSolver.Builder();
            SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural_v2(model, solver);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            int increments = 100;
            var childAnalyzerBuilder = new LoadControlAnalyzer_v2.Builder(model, solver, provider, increments);
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

        public static class EmbeddedEBEModelBuilder
        {
            public static void EmbeddedExampleBuilder(Model_v2 model)
            {
                MatrixModelBuilder(model);
                EBEModelBuilder(model);
                var embeddedGrouping = new EmbeddedGrouping_v2(model, model.ElementsDictionary.Where(x => x.Key <= 90).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key > 90).Select(kv => kv.Value), true);
            }

            public static void MatrixModelBuilder(Model_v2 model)
            {
                string workingDirectory = @"E:\GEORGE_DATA\DESKTOP\input files"; //"..\..\..\Resources\Beam3DInputFiles";

                string MatrixGeometryFileName = "MATRIX_3D-L_x=30-L_y=30-L_z=100-3x3x10-Geometry-2_MSolve.inp";
                
                string MatrixGonnectivityFileName = "MATRIX_3D-L_x=30-L_y=30-L_z=100-3x3x10-ConnMatr_MSolve.inp";
                
                double nodalLoad = 10.0; //0.40;

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
                        model.NodesDictionary.Add(nodeID, new Node { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    }
                }

                // Create Material
                ElasticMaterial3D solidMaterial = new ElasticMaterial3D()
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
                        Element hexa8NLelement = new Element()
                        {
                            ID = elementID,
                            ElementType = new Hexa8NonLinear(solidMaterial, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
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

                // Boundary Conditions
                for (int iNode = 1; iNode <= 16; iNode++) //(int iNode = 1; iNode <= 676; iNode++) //
                {
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.X });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Y });
                    model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.Z });
                    //model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotX });
                    //model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotY });
                    //model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = DOFType.RotZ });
                }

                // Add nodal load values at the top nodes of the model
                for (int iNode = 161; iNode <= 164; iNode++) //(int iNode = 338001; iNode <= 338026; iNode++) //(int iNode = 3601; iNode <= 3606; iNode++) //(int iNode = 2603551; iNode < 2603601; iNode++)
                {
                    model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[iNode], DOF = DOFType.Y });
                }
            }

            public static void EBEModelBuilder(Model_v2 model)
            {
                // define mechanical properties
                double youngModulus = 1.0;
                double shearModulus = 1.0;
                double poissonRatio = 0.30; // 0.034; //(youngModulus / (2 * shearModulus)) - 1;
                double area = 1218.11; //CNT(8,8)-LinearEBE-EBT-L=10nm 
                             //694.77; //CNT(8,8)-LinearEBE-TBT-L=10nm

                double inertiaY = 177.51; //CNT(8,8)-LinearEBE-EBT-L=10nm
                                //100.18; //CNT(8,8)-LinearEBE-TBT-L=10nm

                double inertiaZ = 177.51; //CNT(8,8)-LinearEBE-EBT-L=10nm
                                //100.18; //CNT(8,8)-LinearEBE-TBT-L=10nm

                double torsionalInertia = 168.25; //CNT(8,8)-LinearEBE-EBT-L=10nm
                                        //68.77 //CNT(8,8)-LinearEBE-TBT-L=10nm

                double effectiveAreaY = area;
                double effectiveAreaZ = area;
                string workingDirectory = @"E:\GEORGE_DATA\DESKTOP\input files"; //"..\..\..\Resources\Beam3DInputFiles";

                string CNTgeometryFileName = "EmbeddedCNT-8-8-L=100-h=3-k=1-EBE-L=10-NumberOfCNTs=1-Geometry_beam.inp";
                
                string CNTconnectivityFileName = "EmbeddedCNT-8-8-L=100-h=3-k=1-EBE-L=10-NumberOfCNTs=1-ConnMatr_beam.inp";
                
                int CNTNodes = File.ReadLines(workingDirectory + '\\' + CNTgeometryFileName).Count();
                int CNTElems = File.ReadLines(workingDirectory + '\\' + CNTconnectivityFileName).Count();

                // Geometry
                using (TextReader reader = File.OpenText(workingDirectory + '\\' + CNTgeometryFileName))
                {
                    for (int i = 0; i < CNTNodes; i++)
                    {
                        string text = reader.ReadLine();
                        string[] bits = text.Split(',');
                        int nodeID = int.Parse(bits[0]) + 176; // matrixNodes
                        double nodeX = double.Parse(bits[1]);
                        double nodeY = double.Parse(bits[2]);
                        double nodeZ = double.Parse(bits[3]);
                        model.NodesDictionary.Add(nodeID, new Node { ID = nodeID, X = nodeX, Y = nodeY, Z = nodeZ });
                    }
                }

                // Create new 3D material
                ElasticMaterial3D beamMaterial = new ElasticMaterial3D
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
                        int elementID = int.Parse(bits[0]) + 90; // matrixElements
                        int node1 = int.Parse(bits[1]) + 176; // matrixNodes
                        int node2 = int.Parse(bits[2]) + 176; // matrixNodes
                                                    
                        // element nodes
                        IList<Node> elementNodes = new List<Node>();
                        elementNodes.Add(model.NodesDictionary[node1]);
                        elementNodes.Add(model.NodesDictionary[node2]);
                        // create element
                        var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, beamMaterial, 7.85, beamSection);
                        var beamElement = new Element { ID = elementID, ElementType = beam_1 };
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
