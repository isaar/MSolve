using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.Analysis;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip;
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
        private static readonly string expectedDofEnumerationPath = @"../../Resources/DofEnumeration.txt";


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

        private void HandleConstraints(Model2D model)
        {
            XNode2D bottomRight = nodes[NODE_COLUMNS - 1];
            XNode2D topRight = nodes[NODE_ROWS * NODE_COLUMNS -1];
            model.AddConstraint(bottomRight, StandardDOFType.X);
            model.AddConstraint(bottomRight, StandardDOFType.Y);
            model.AddConstraint(topRight, StandardDOFType.X);
            model.AddConstraint(topRight, StandardDOFType.Y);
        }

        private void CreateLoads(Model2D model)
        {
            XNode2D topLeft = nodes[(NODE_ROWS - 1) * NODE_COLUMNS];
            model.AddNodalLoad(topLeft, StandardDOFType.Y, 1000.0);
        }

        private void HandleEnrichments() // Most of these should not be done manually
        {
            // Create enrichments
            var crackStart = new CartesianPoint2D(0.00, 0.18);
            var intersection = new CartesianPoint2D(0.15, 0.1575);
            var crackEnd = new CartesianPoint2D(0.20, 0.15);
            var polyline = new Polyline2D(crackStart, crackEnd);
            crackBody = new CrackBody2D(polyline);
            crackTip = new CrackTip2D(CrackTip2D.TipCurvePosition.CurveEnd, polyline, new SingleElementEnrichment(2.0),
                new HomogeneousMaterialAuxiliaryStates(), new HomogeneousSIFCalculator(material));

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
            for (int el = 0; el < elements.Length; ++el) model.AddElement(elements[el]);

            HandleConstraints(model);
            CreateLoads(model);

            HandleEnrichments();

            model.EnumerateDofs();
        }

        public static void CheckSolution(Model2D model)
        {
            LinearStaticAnalysis analysis = new LinearStaticAnalysis(model);
            analysis.Solve();
            analysis.PrintSolution();
        }

        public static void Main()
        {
            bool mockGaussPoints = true;
            Benchmark1 benchmark = new Benchmark1(mockGaussPoints);
            benchmark.CreateModel();

            // Gauss Points
            var gpChecker = new GaussPointChecker(expectedGaussPointsPath, 1.0e-8, true);
            //gpChecker.CheckElementGaussPoints(benchmark.elements);

            // Element stiffness matrices
            var elementChecker = new ElementMatrixChecker(expectedElementMatricesPath, 1.0e-8, false);
            //elementChecker.CheckElementMatrices(benchmark.elements);

            // Global stiffness matrix
            //DofEnumerationChecker.PrintEnumeration(benchmark.model);
            var globalChecker = new GlobalMatrixChecker(expectedGlobalMatrixPath, expectedDofEnumerationPath, 1.0e-5, true);
            globalChecker.CheckGlobalMatrix(benchmark.model);
            //globalChecker.PrintGlobalMatrix(benchmark.model);

            // Solution
            //CheckSolution(benchmark.model);
            //double[] forces = benchmark.model.CalculateForces();
            //Console.WriteLine("Forces:");
            //for (int i = 0; i < forces.Length; ++i) Console.WriteLine(i + ": " + forces[i]);
        }
    }
}
