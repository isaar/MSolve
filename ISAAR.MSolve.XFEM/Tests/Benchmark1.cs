using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tests.Tools;

namespace ISAAR.MSolve.XFEM.Tests
{
    class Benchmark1
    {
        private readonly static int NODE_ROWS = 4;
        private readonly static int NODE_COLUMNS = 5;
        private readonly static int ELEMENT_ROWS = 3;
        private readonly static int ELEMENT_COLUMNS = 4;

        private static readonly string expectedGaussPointsPath = @"../../Resources/GaussPoints.txt";
        private static readonly string expectedElementMatricesPath = @"../../Resources/ElementMatrices.txt";
        private static readonly string expectedGlobalMatrixPath = @"../../Resources/GlobalMatrix.txt";

        private Model2D model;
        private IFiniteElementMaterial2D material;
        private XNode2D[] nodes;
        private XContinuumElement2D[] elements;
        private CrackBody2D crackBody;
        private CrackTip2D crackTip;

        private MockIntegrationStrategy[] mockIntegrations;
        private bool mockIntegration;

        public Benchmark1(bool mockIntegration = false)
        {
            this.mockIntegration = mockIntegration;
        }

        private void HandleIntegration()
        {
            int elementsCount = ELEMENT_ROWS * ELEMENT_COLUMNS;
            GaussPoint2D[][] allPoints = OutputReaders.ReadAllGaussPoints(expectedGaussPointsPath, elementsCount);
            mockIntegrations = new MockIntegrationStrategy[elementsCount];
            for (int el = 0; el < elementsCount; ++el)
            {
                mockIntegrations[el] = new MockIntegrationStrategy(allPoints[el], material);
            }
        }

        private void CreateNodes()
        {
            nodes = new XNode2D[NODE_ROWS * NODE_COLUMNS];
            int id = 0;
            double dx = 0.15;
            double dy = 0.10;
            for (int row = 0; row < NODE_ROWS; ++row)
            {
                for (int col = 0; col < NODE_COLUMNS; ++col)
                {
                    nodes[id] = new XNode2D(id, col * dx, row * dy);
                    ++id;
                }
            }
        }
        
        private void CreateElements()
        {
            elements = new XContinuumElement2D[ELEMENT_ROWS * ELEMENT_COLUMNS];
            int id = 0;
            for (int row = 0; row < ELEMENT_ROWS; ++row)
            {
                for (int col = 0; col < ELEMENT_COLUMNS; ++col)
                {
                    int firstNode = row * NODE_COLUMNS + col;
                    XNode2D[] elementNodes = { nodes[firstNode], nodes[firstNode+1],
                        nodes[firstNode + NODE_COLUMNS + 1], nodes[firstNode + NODE_COLUMNS] };
                    if (mockIntegration)
                    {
                        elements[id] = new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes,
                            mockIntegrations[id]);
                    }
                    else
                    {
                        var enrIntegrationRule = new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order7Points13,
                            new IncrementalTriangulator());
                        elements[id] = new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes,
                            new HomogeneousIntegration2D(enrIntegrationRule, material));
                    }
                    id++;
                }
            }
        }

        private void HandleEnrichments() // Most of these should not be done manually
        {
            // Create enrichments
            var crackStart = new CartesianPoint2D(0.00, 0.18);
            var intersection = new CartesianPoint2D(0.15, 0.1575);
            var crackEnd = new CartesianPoint2D(0.20, 0.15);
            var polyline = new Polyline2D(crackStart, crackEnd);
            crackBody = new CrackBody2D(polyline);
            crackTip = new CrackTip2D(CrackTip2D.TipPosition.CurveEnd, polyline);

            // Mesh geometry interaction
            polyline.ElementIntersections.Add(elements[4], new CartesianPoint2D[] { crackStart, intersection });
            polyline.ElementIntersections.Add(elements[5], new CartesianPoint2D[] { intersection });

            // Enrich nodes
            crackBody.EnrichNode(nodes[5]);
            crackBody.EnrichNode(nodes[10]);
            crackTip.EnrichNode(nodes[6]);
            crackTip.EnrichNode(nodes[7]);
            crackTip.EnrichNode(nodes[11]);
            crackTip.EnrichNode(nodes[12]);

            //Enrich elements. 
            crackBody.EnrichElement(elements[4]);
            crackTip.EnrichElement(elements[5]);
        }

        private void CreateModel()
        {
            material = ElasticMaterial2DPlainStrain.Create(2.0e6, 0.3, 1.0);
            model = new Model2D();

            CreateNodes();
            foreach (var node in nodes) model.AddNode(node);

            HandleIntegration();
            CreateElements();
            for (int el = 0; el < elements.Length; ++el) model.AddElement(new Element2D(el, elements[el]));

            HandleEnrichments();
        }

        private void CheckGaussPoints()
        {
            var checker = new GaussPointChecker(expectedGaussPointsPath, 1.0e-8, true);
            checker.CheckElementGaussPoints(elements);
        }

        private void CheckElementStiffnessMatrices()
        {
            var checker = new ElementMatrixChecker(expectedElementMatricesPath, 1.0e-8, false);
            checker.CheckElementMatrices(elements);
        }

        private Matrix2D<double> ReadGlobalStiffnessMatrix()
        {
            string[] lines = File.ReadAllLines(expectedGlobalMatrixPath);
            int order = lines.Length;
            string lastLine = lines[lines.Length - 1];
            if (string.IsNullOrEmpty(lastLine) || string.IsNullOrWhiteSpace(lastLine)) --order; //Last row might have output an extra newline

            var matrix = new Matrix2D<double>(order, order);
            for (int row = 0; row <order; ++row)
            {
                string[] words = lines[row].Split(' ');
                for (int col = 0; col <order; ++col)
                {
                    matrix[row, col] = double.Parse(words[col]);
                }
            }
           
            return matrix;
        }


        //private void CheckGlobalStiffnessMatrix()
        //{
        //    Console.WriteLine("Checking global stiffness matrix...");
        //    string header = "Global matrix: ";

        //    // Retrieve the matrices
        //    Matrix2D<double> correctMatrix = ReadGlobalStiffnessMatrix();
        //    SkylineMatrix2D<double> globalMatrix = SingleGlobalSkylineAssembler.BuildGlobalMatrix(model);

        //    // Check dimensions first
        //    if (globalMatrix.Rows != correctMatrix.Rows) throw new ArgumentException(header +"Non matching rows.");
        //    if (globalMatrix.Columns != correctMatrix.Columns) throw new ArgumentException(header 
        //        + "Non matching columns.");

        //    // Check each entry
        //    for (int row = 0; row < globalMatrix.Rows; ++row)
        //    {
        //        for (int col = 0; col < globalMatrix.Columns; ++col)
        //        {
        //            if (!AreEqual(globalMatrix[row, col], correctMatrix[row, col]))
        //                throw new ArgumentException(header + "Error at K[" + row + ", " + col + "]");
        //        }
        //    }
        //    Console.WriteLine("Global stiffness matrix is correct!");
        //}

        private void CheckSolutionVector()
        {

        }

        private void Extra()
        {
            var element = elements[5];
            SymmetricMatrix2D<double> kss, kee;
            Matrix2D<double> kes;
            kss = element.BuildStandardStiffnessMatrix();
            element.BuildEnrichedStiffnessMatrices(out kes, out kee);
            MatrixPrinter.PrintElementMatrices(5, kss, kes, kee);
        }

        public static void Main()
        {
            Benchmark1 benchmark = new Benchmark1(false);
            benchmark.CreateModel();

            benchmark.CheckGaussPoints();
            //benchmark.CheckElementStiffnessMatrices();
            //benchmark.CheckGlobalStiffnessMatrix();
        }
    }
}
