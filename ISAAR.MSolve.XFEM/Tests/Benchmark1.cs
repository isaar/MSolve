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
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests
{
    class Benchmark1
    {
        private static int NODE_ROWS = 4;
        private static int NODE_COLUMNS = 5;
        private static int ELEMENT_ROWS = 3;
        private static int ELEMENT_COLUMNS = 4;
        private static double TOLERANCE = 1e-4;

        private Model2D model;
        private IFiniteElementMaterial2D material;
        private XNode2D[] nodes;
        private XElement2D[] elements;
        private CrackBody2D crackBody;
        private CrackTip2D crackTip;

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
            IIntegrationRule2D integrationRule = new RectangularSubgridIntegration2D(2, GaussLegendre2D.Order2x2);
            var integrationFactory = new HomogeneousIntegration2D.Factory(integrationRule, material);

            elements = new XElement2D[ELEMENT_ROWS * ELEMENT_COLUMNS];
            int id = 0;
            for (int row = 0; row < ELEMENT_ROWS; ++row)
            {
                for (int col = 0; col < ELEMENT_COLUMNS; ++col)
                {
                    int firstNode = row * NODE_COLUMNS + col;
                    XNode2D[] elementNodes = { nodes[firstNode], nodes[firstNode+1],
                        nodes[firstNode + NODE_COLUMNS + 1], nodes[firstNode + NODE_COLUMNS] };
                    elements[id++] = new XElement2D(new IsoparametricQuad4(elementNodes, integrationFactory));
                }
            }
        }

        private void HandleEnrichments() // Most of these should not be done manually
        {
            // Create enrichments
            var crackStart = new CartesianPoint2D(0.00, 0.15);
            var crackEnd = new CartesianPoint2D(0.20, 0.15);
            crackBody = new CrackBody2D(new Line2D(crackStart, crackEnd));
            crackTip = new CrackTip2D(crackEnd, 0.0);

            // Enrich nodes
            crackBody.EnrichNode(nodes[5]);
            crackBody.EnrichNode(nodes[10]);
            crackTip.EnrichNode(nodes[6]);
            crackTip.EnrichNode(nodes[7]);
            crackTip.EnrichNode(nodes[11]);
            crackTip.EnrichNode(nodes[12]);

            //Enrich elements
            crackBody.EnrichElement(elements[4]);
            crackTip.EnrichElement(elements[5]);
        }

        private void CreateModel()
        {
            material = ElasticMaterial2DPlainStrain.Create(2.0e6, 0.3, 1.0);
            model = new Model2D();

            CreateNodes();
            foreach (var node in nodes) model.AddNode(node);

            CreateElements();
            for (int el = 0; el < elements.Length; ++el) model.AddElement(new Element2D(el, elements[el]));

            HandleEnrichments();
        }

        private static void BlendingElement0() // Delete this
        {
            Benchmark1 benchmark = new Benchmark1();
            benchmark.CreateModel();

            SymmetricMatrix2D<double> kss, kee;
            Matrix2D<double> kes;
            kss = benchmark.elements[0].BuildStdStiffnessMatrix();
            benchmark.elements[0].BuildEnrichedStiffnessMatrices(out kes, out kee);

            Console.WriteLine("k0ss = ");
            Console.WriteLine(kss);
            Console.WriteLine("k0es = ");
            Console.WriteLine(kes);
            Console.WriteLine("k0ee =  ");
            Console.WriteLine(kee);
        }

        private bool AreEqual(double actual, double expected)
        {
            if (Math.Abs(expected) < TOLERANCE && Math.Abs(actual) < TOLERANCE) return true;
            else if (Math.Abs(expected) < TOLERANCE && Math.Abs(actual) >= TOLERANCE) return false;
            else return (Math.Abs(1.0 - actual / expected) < TOLERANCE) ? true : false;
        }

        private Matrix2D<double>[] ReadElementStiffnessMatrices()
        {
            string path = @"../../Resources/ElementMatrices.txt";
            string[] lines = File.ReadAllLines(path);
            var elementMatrices = new Matrix2D<double>[elements.Length];

            const int START = -1;
            int currentElement = START;
            List<List<double>> currentMatrix = null;
            foreach (string line in lines)
            {
                string[] words = line.Split(' ');
                if (line[0] == '*') // New matrix
                {
                    // Build and store the finished matrix
                    if (currentElement != START)
                    {
                        int rows = currentMatrix.Count;
                        int cols = currentMatrix[0].Count;
                        Matrix2D<double> finishedMatrix = new Matrix2D<double>(rows, cols);
                        for (int r = 0; r < rows; ++r)
                        {
                            for (int c = 0; c < cols; ++c)
                            {
                                finishedMatrix[r, c] = currentMatrix[r][c];
                            }
                        }
                        elementMatrices[currentElement] = finishedMatrix;
                    }
                    if (words.Length > 1) // The last line only has a *.
                    {
                        currentElement = int.Parse(words[1]);
                        currentMatrix = new List<List<double>>();
                    }
                }
                else
                {
                    // Read a new row of doubles
                    var newRow = new List<double>();
                    foreach (string number in words)
                    {
                        newRow.Add(double.Parse(number));
                    }
                    currentMatrix.Add(newRow);
                }
            }
            return elementMatrices;
        }

        private void CheckElementStiffnessMatrices()
        {
            Console.WriteLine("Checking element stiffness matrices...");

            Matrix2D<double>[] correctMatrices = ReadElementStiffnessMatrices();
            for (int el = 0; el < elements.Length; ++el)
            {
                string header = "Element " + el + ": ";

                // Retrieve the matrices
                Matrix2D<double> correctK = correctMatrices[el];
                XElement2D element = elements[el];
                SymmetricMatrix2D<double> kss, kee;
                Matrix2D<double> kes;
                kss = element.BuildStdStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out kes, out kee);

                // Check dimensions first
                //int stdDofsCount = element.StandardFiniteElement.DofsCount;
                //int enrDofsCount = element.CountArtificialDofs();
                if (kss.Rows + kes.Rows != correctK.Rows) throw new ArgumentException(header +"Non matching rows.");
                if (kss.Columns + kee.Columns != correctK.Columns) throw new ArgumentException(header 
                    + "Non matching columns.");

                // Check Kss entrywise
                for (int row = 0; row < kss.Rows; ++row)
                {
                    for (int col = 0; col < kss.Columns; ++col)
                    {
                        if (!AreEqual(kss[row, col], correctK[row, col]))
                            throw new ArgumentException(header + "Error at Kss[" + row + ", " + col + "]");
                    }
                }

                //    // Check Kes entrywise
                //    for (int row = 0; row < kes.Rows; ++row)
                //    {
                //        for (int col = 0; col < kes.Columns; ++col)
                //        {
                //            if (!AreEqual(kes[row, col], correctK[kss.Rows + row, col]))
                //                throw new ArgumentException("Error at Kes[" + row + ", " + col + "]");
                //        }
                //    }

                //    // Check Kee entrywise
                //    for (int row = 0; row < kee.Rows; ++row)
                //    {
                //        for (int col = 0; col < kee.Columns; ++col)
                //        {
                //            if (!AreEqual(kee[row, col], correctK[kss.Rows + row, kss.Columns + col]))
                //                throw new ArgumentException("Error at Kee[" + row + ", " + col + "]");
                //        }
                //    }
            }
            Console.WriteLine("Element stiffness matrices are correct!");
        }


        private Matrix2D<double> ReadGlobalStiffnessMatrix()
        {
            string path = @"../../Resources/GlobalMatrix.txt";
            string[] lines = File.ReadAllLines(path);

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


        private void CheckGlobalStiffnessMatrix()
        {
            Console.WriteLine("Checking global stiffness matrix...");
            string header = "Global matrix: ";

            // Retrieve the matrices
            Matrix2D<double> correctMatrix = ReadGlobalStiffnessMatrix();
            SkylineMatrix2D<double> globalMatrix = SingleGlobalSkylineAssembler.BuildGlobalMatrix(model);

            // Check dimensions first
            if (globalMatrix.Rows != correctMatrix.Rows) throw new ArgumentException(header +"Non matching rows.");
            if (globalMatrix.Columns != correctMatrix.Columns) throw new ArgumentException(header 
                + "Non matching columns.");

            // Check each entry
            for (int row = 0; row < globalMatrix.Rows; ++row)
            {
                for (int col = 0; col < globalMatrix.Columns; ++col)
                {
                    if (!AreEqual(globalMatrix[row, col], correctMatrix[row, col]))
                        throw new ArgumentException(header + "Error at K[" + row + ", " + col + "]");
                }
            }
            Console.WriteLine("Global stiffness matrix is correct!");
        }

        private void CheckSolutionVector()
        {

        }

        public static void Main()
        {
            Benchmark1 benchmark = new Benchmark1();
            benchmark.CreateModel();

            benchmark.CheckElementStiffnessMatrices();
            benchmark.CheckGlobalStiffnessMatrix();
        }
    }
}
