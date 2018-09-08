using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Tests.Tools;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    /// <summary>
    /// See Khoei (2015), 7.6.1.
    /// Main units: cm, kN
    /// </summary>
    class DCB3x1
    {
        private static readonly bool integrationWithTriangles = false;
        private static Matrix expectedK_Node6;
        private static Matrix expectedK_Node7_El1;
        private static Matrix expectedK_Node7_El2;
        private static Matrix expectedK_Node7_Global;
        private static double[] expectedSolutionNodes5_6 =
            { -8.17e-3, -50e-3, -8.17e-3, 50e-3, 15.69e-3, 49.88e-3, -15.69e-3, 49.88e-3 };

        private static void StiffnessNode6()
        {
            // Expected entries of element 2 stiffness submatrix, corresponding to node 6
            var array = new double[,]
            {
                { 1.154, 0.481, -0.481, -0.240 },
                { 0.481, 1.154, -0.240, -0.962 },
                { -0.481, -0.240, 0.962, 0.481 },
                { -0.240, -0.962, 0.481, 1.923 }
            };
            expectedK_Node6 = Matrix.CreateFromArray(array);
            expectedK_Node6.ScaleIntoThis(1e6);
        }

        private static void StiffnessNode7Element1()
        {
            // Expected entries of element 1 stiffness submatrix, corresponding to node 7
            var array = new double[,]
            {
                { 1.154, 0.481, 1.568, 0.544, -0.444, 0.12, -0.847, 0.016, 0.378, -0.337 },
                { 0.481, 1.154, 0.575, 2.668, -0.114, 0.175, -0.165, -0.271, -0.055, -0.358 },
                { 1.568, 0.575, 12.432, 4.896, -0.824, -1.69, -3.114, -3.125, 0.134, 0.537 },
                { 0.544, 2.668, 4.896, 17.018, -1.359, -3.322, -2.444, -6.366, 0.465, 3.07 },
                { -0.444, -0.114, -0.824, -1.359, 2.639, -0.459, 1.648, -0.253, -1.869, 0.939 },
                { 0.12, 0.175, -1.69, -3.322, -0.459, 3.921, -0.214, 4.699, 0.909, -2.729 },
                { -0.847, -0.165, -3.114, -2.444, 1.648, -0.214, 4.063, -0.01, -1.386, 0.645 },
                { 0.016, -0.271, -3.125, -6.366, -0.253, 4.699, -0.01, 7.896, 0.685, -2.804 },
                { 0.378, -0.055, 0.134, 0.465, -1.869, 0.909, -1.386, 0.685, 3.081, -0.859 },
                { -0.337, -0.358, 0.537, 3.07, 0.939, -2.729, 0.645, -2.804, -0.859, 4.694 }
            };
            expectedK_Node7_El1 = Matrix.CreateFromArray(array);
            expectedK_Node7_El1.ScaleIntoThis(1e6);
        }

        private static void StiffnessNode7Element2()
        {
            // Expected entries of element 2 stiffness submatrix, corresponding to node 7
            var array = new double[,]
            {
               { 1.154, -0.481, 1.715, -0.752, -0.586, 0.006, -0.886, -0.052, 0.600, -0.087 },
               { -0.481, 1.154, -0.856, 3.184, 0.092, -0.325, 0.106, -0.417, -0.156, 0.417 },
               { 1.715, -0.856, 13.871, -5.198, -2.050, 0.718, -3.537, 1.216, 1.561, -0.610 },
               { -0.752, 3.184, -5.198, 26.465, 0.679, -5.280, 1.135, -9.091, -0.753, 3.713 },
               { -0.586, 0.092, -2.050, 0.679, 1.098, -0.169, 1.822, -0.235, -0.843, 0.223 },
               { 0.006, -0.325, 0.718, -5.280, -0.169, 2.052, -0.246, 3.583, 0.234, -1.222 },
               { -0.886, 0.106, -3.537, 1.135, 1.822, -0.246, 3.075, -0.350, -1.317, 0.322 },
               { -0.052, -0.417, 1.216, -9.091, -0.235, 3.583, -0.350, 6.348, 0.330, -1.973 },
               { 0.600, -0.156, 1.561, -0.753, -0.843, 0.234, -1.317, 0.330, 0.869, -0.282 },
               { -0.087, 0.417, -0.610, 3.713, 0.223, -1.222, 0.322, -1.973, -0.282, 1.180 }
            };
            expectedK_Node7_El2 = Matrix.CreateFromArray(array);
            expectedK_Node7_El2.ScaleIntoThis(1e6);
        }

        private static void StiffnessNode7Global()
        {
            // Expected entries of global stiffness submatrix, corresponding to node 7
            var array = new double[,]
            {
                { 2.308, 0.000, 3.283, -0.208, -1.030, 0.126, -1.733, -0.036, 0.979, -0.424 },
                { 0.000, 2.308, -0.282, 5.852, -0.022, -0.150, -0.060, -0.687, -0.211, 0.059 },
                { 3.283, -0.282, 26.303, -0.302, -2.874, -0.972, -6.651, -1.910, 1.695, -0.073 },
                { -0.208, 5.852, -0.302, 43.483, -0.680, -8.601, -1.310, -15.456, -0.289, 6.783 },
                { -1.030, -0.022, -2.874, -0.680, 3.736, -0.628, 3.470, -0.488, -2.712, 1.162 },
                { 0.126, -0.150, -0.972, -8.601, -0.628, 5.973, -0.460, 8.282, 1.142, -3.951 },
                { -1.733, -0.060, -6.651, -1.310, 3.470, -0.460, 7.139, -0.360, -2.703, 0.966 },
                { -0.036, -0.687, -1.910, -15.456, -0.488, 8.282, -0.360, 14.244, 1.015, -4.777 },
                { 0.979, -0.211, 1.695, -0.289, -2.712, 1.142, -2.703, 1.015, 3.950, -1.141 },
                { -0.424, 0.059, -0.073, 6.783, 1.162, -3.951, 0.966, -4.777, -1.141, 5.874 }
            };
            expectedK_Node7_Global = Matrix.CreateFromArray(array);
            expectedK_Node7_Global.ScaleIntoThis(1e6);
        }

        static DCB3x1()
        {
            StiffnessNode6();
            StiffnessNode7Element1();
            StiffnessNode7Element2();
            StiffnessNode7Global();
        }

        public static void Run()
        {
            DCB3x1 benchmark = new DCB3x1(20);
            benchmark.CreateMaterial();
            Model2D model = benchmark.CreateModel();
            benchmark.HandleEnrichment(model);
            //ISolver solver = new DenseSolver(model);
            //ISolver solver = new SkylineSolverOLD(model);
            ISolver solver = new SkylineSolver(model);
            //ISolver solver = new CholeskySuiteSparseSolver(model);
            //ISolver solver = new PCGSolver(0.8, 1e-6);
            solver.Initialize();
            solver.Solve();

            benchmark.CheckStiffnessNode6(model);
            benchmark.CheckStiffnessNode7Element1(model);
            benchmark.CheckStiffnessNode7Element2(model);
            benchmark.CheckGlobalStiffnessNode7(model, solver);
            benchmark.CheckSolution(model, solver);

            benchmark.PrintAllStiffnesses(model);
            benchmark.PrintDisplacements(model, solver);
        }

        private readonly double h; // element length
        private readonly double E = 2e6, v = 0.3;
        private readonly SubmatrixChecker checker;
        private HomogeneousElasticMaterial2D globalHomogeneousMaterial;
        private BasicExplicitCrack2D crack;


        public DCB3x1(double elementLength)
        {
            this.h = elementLength;
            this.checker = new SubmatrixChecker(1e-6); // The heaviside dofs are slightly off
        }

        private void CreateMaterial()
        {
            globalHomogeneousMaterial = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
        }

        public Model2D CreateModel()
        {
            var model = new Model2D();

            //Nodes
            model.AddNode(new XNode2D(0, 0.0, 0.0));
            model.AddNode(new XNode2D(1, h, 0.0));
            model.AddNode(new XNode2D(2, h, h));
            model.AddNode(new XNode2D(3, 0.0, h));
            model.AddNode(new XNode2D(4, 2 * h, 0.0));
            model.AddNode(new XNode2D(5, 3 * h, 0.0));
            model.AddNode(new XNode2D(6, 3 * h, h));
            model.AddNode(new XNode2D(7, 2 * h, h));

            // Elements
            IReadOnlyList<XNode2D> nodes = model.Nodes;
            XNode2D[][] connectivity = new XNode2D[3][];
            connectivity[0] = new XNode2D[] { nodes[0], nodes[1], nodes[2], nodes[3] };
            connectivity[1] = new XNode2D[] { nodes[1], nodes[4], nodes[7], nodes[2] };
            connectivity[2] = new XNode2D[] { nodes[4], nodes[5], nodes[6], nodes[7] };

            IIntegrationStrategy2D<XContinuumElement2D> integration, jIntegration;
            if (integrationWithTriangles)
            {
                ITriangulator2D triangulator = new IncrementalTriangulator();
                integration = new IntegrationForCrackPropagation2D(
                    new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order2Points3, crack, triangulator),
                    new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order2Points3, crack, triangulator));
                jIntegration = new IntegrationWithSubtriangles(GaussQuadratureForTriangle.Order3Points4, crack,
                    triangulator);
            }
            else
            {
                integration = new IntegrationForCrackPropagation2D(
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
                jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);
            }

            for (int e = 0; e < 3; ++e)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    connectivity[e], materialField, integration, jIntegration));
            }

            // Boundary conditions
            model.AddConstraint(model.Nodes[0], DisplacementDof.X, 0.0);
            model.AddConstraint(model.Nodes[0], DisplacementDof.Y, 0.0);
            model.AddConstraint(model.Nodes[3], DisplacementDof.X, 0.0);
            model.AddConstraint(model.Nodes[3], DisplacementDof.Y, 0.0);
            model.AddConstraint(model.Nodes[5], DisplacementDof.Y, -0.05);
            model.AddConstraint(model.Nodes[6], DisplacementDof.Y, 0.05);

            return model;
        }

        public void HandleEnrichment(Model2D model)
        {
            crack = new BasicExplicitCrack2D();
            var boundary = new RectangularBoundary(0, 3 * h, 0, h);
            crack.Mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements, boundary);

            // Create enrichments          
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack, new SignFunctionOpposite2D());
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);
            //crackTip = new CrackTip2D(CrackTip2D.TipCurvePosition.CurveStart, polyline, new SingleElementEnrichment(),
            //    2.0, new HomogeneousMaterialAuxiliaryStates(globalHomogeneousMaterial),
            //    new HomogeneousSIFCalculator(globalHomogeneousMaterial));

            // Mesh geometry interaction
            var crackMouth = new CartesianPoint2D(3 * h, 0.5 * h);
            var crackTip = new CartesianPoint2D(1.5 * h, 0.5 * h);
            crack.InitializeGeometry(crackMouth, crackTip);
            crack.UpdateEnrichments();
        }

        private void CheckStiffnessNode6(Model2D model)
        {
            Matrix Kss = model.Elements[2].BuildStandardStiffnessMatrix();
            (Matrix Kee, Matrix Kes) = model.Elements[2].BuildEnrichedStiffnessMatricesLower();
            //PrintElementMatrices(Kss, Kes, Kee, 1e-6);

            int[] stdDofsNode6 = new int[] { 4, 5 }; // All nodes: 0-1, 2-3, 4-5, 6-7 (2 dofs/node)
            int[] enrDofsNode6 = new int[] { 10, 11 }; // All nodes: 0-7, 8-9, 10-11, 12-19 (2 or 8 dofs/node)
            Matrix node6Submatrix = checker.ExtractRelevant(Kss, Kes, Kee, stdDofsNode6, enrDofsNode6);
            RoundMatrix(node6Submatrix);
            Console.Write("Checking stiffness submatrix of node 6 - element 2: ");
            checker.Check(expectedK_Node6, node6Submatrix);
        }

        private void CheckStiffnessNode7Element1(Model2D model)
        {
            Matrix Kss = model.Elements[1].BuildStandardStiffnessMatrix();
            (Matrix Kee, Matrix Kes) = model.Elements[1].BuildEnrichedStiffnessMatricesLower();
            //PrintElementMatrices(Kss, Kes, Kee, 1e-6);

            int[] stdDofsNode7 = new int[] { 4, 5 }; // All nodes: 0-1, 2-3, 4-5, 6-7 (2 dofs/node)
            int[] enrDofsNode7 = new int[] { 16, 17, 18, 19, 20, 21, 22, 23 }; // All nodes: 0-7, 8-15, 16-23, 24-31 (8 dofs/node)
            Matrix node7Submatrix = checker.ExtractRelevant(Kss, Kes, Kee, stdDofsNode7, enrDofsNode7);
            RoundMatrix(node7Submatrix);
            Console.Write("Checking stiffness submatrix of node 7 - element 1: ");
            checker.Check(expectedK_Node7_El1, node7Submatrix);
        }

        private void CheckStiffnessNode7Element2(Model2D model)
        {
            Matrix Kss = model.Elements[2].BuildStandardStiffnessMatrix();
            (Matrix Kee, Matrix Kes) = model.Elements[2].BuildEnrichedStiffnessMatricesLower();
            //PrintElementMatrices(Kss, Kes, Kee, 1e-6);

            int[] stdDofsNode7 = new int[] { 6, 7 }; // All nodes: 0-1, 2-3, 4-5, 6-7 (2 dofs/node)
            int[] enrDofsNode7 = new int[] { 12, 13, 14, 15, 16, 17, 18, 19 }; // All nodes: 0-7, 8-9, 10-11, 12-19 (2 or 8 dofs/node)
            Matrix node7Submatrix = checker.ExtractRelevant(Kss, Kes, Kee, stdDofsNode7, enrDofsNode7);
            RoundMatrix(node7Submatrix);
            Console.Write("Checking stiffness submatrix of node 7 - element 2: ");
            checker.Check(expectedK_Node7_El2, node7Submatrix);
        }

        private void CheckGlobalStiffnessNode7(Model2D model, ISolver solver)
        {
            (DokSymmetric Kuu, DokRowMajor Kuc) = (new GlobalDOKAssembler()).BuildGlobalMatrix(model, solver.DofOrderer);

            var dofsOfNode7 = new List<int>();
            foreach (int dof in solver.DofOrderer.GetStandardDofsOf(model.Nodes[7])) dofsOfNode7.Add(dof);
            foreach (int dof in solver.DofOrderer.GetEnrichedDofsOf(model.Nodes[7])) dofsOfNode7.Add(dof);
            dofsOfNode7.Sort();

            Matrix node7Submatrix = checker.ExtractRelevant(Kuu, dofsOfNode7);
            RoundMatrix(node7Submatrix);
            Console.Write("Checking stiffness submatrix of node 7 - global: ");
            //Console.WriteLine("Knode7 = ");
            //MatrixUtilities.PrintDense(node7Submatrix);
            checker.Check(expectedK_Node7_Global, node7Submatrix);
        }

        // WARNING: The heaviside enriched dofs are slightly off!
        private void CheckSolution(Model2D model, ISolver solver)
        {
            var displacements = new List<double>();

            // Standard dofs
            int dofUx5 = solver.DofOrderer.GetStandardDofOf(model.Nodes[5], DisplacementDof.X);
            displacements.Add(solver.Solution[dofUx5]);
            displacements.Add(-0.05);
            int dofUx6 = solver.DofOrderer.GetStandardDofOf(model.Nodes[6], DisplacementDof.X);
            displacements.Add(solver.Solution[dofUx6]);
            displacements.Add(0.05);

            // Enriched dofs Concatenate IEnumerables and then toArray()
            var enrichedDofs = solver.DofOrderer.GetEnrichedDofsOf(model.Nodes[5]).
                Concat(solver.DofOrderer.GetEnrichedDofsOf(model.Nodes[6]));
            foreach (int dof in enrichedDofs) displacements.Add(solver.Solution[dof]);
            RoundList(displacements);
            Console.Write("Checking the solution vector for nodes 5, 6: ");
            checker.Check(expectedSolutionNodes5_6, displacements);
        }

        private void PrintAllStiffnesses(Model2D model)
        {
            Console.WriteLine("\n-------------------- Element Matrices ------------------");
            for (int i = 0; i < model.Elements.Count; ++i)
            {
                Matrix Kss = model.Elements[i].BuildStandardStiffnessMatrix();
                (Matrix Kee, Matrix Kes) = model.Elements[i].BuildEnrichedStiffnessMatricesLower();

                Console.WriteLine("Element " + i);
                PrintElementMatrices(Kss, Kes, Kee, 1e-6);
            }
            Console.WriteLine();
        }

        private void PrintDisplacements(Model2D model, ISolver solver)
        {
            double[,] nodalDisplacements = solver.DofOrderer.GatherNodalDisplacements(model, solver.Solution);
           
            Console.WriteLine("\n-------------------- Standard displacements ------------------");
            for (int i = 0; i < nodalDisplacements.GetLength(0); ++i)
            {
                Console.WriteLine("Node " + i + ", dof X :\t Ux = " + nodalDisplacements[i, 0]);
                Console.WriteLine("Node " + i + ", dof Y :\t Uy = " + nodalDisplacements[i, 1]);
            }
            Console.WriteLine();

            Console.WriteLine("\n-------------------- Artificial displacements ------------------");
            Console.WriteLine(solver.DofOrderer.GatherEnrichedNodalDisplacements(model, solver.Solution));
        }

        private void PrintElementMatrices(Matrix kss, Matrix kes, Matrix kee, double scale)
        {
            var writer = new FullMatrixWriter();

            Console.WriteLine("Kss = " + 1/scale + " * ");
            writer.WriteToConsole(kss.Scale(scale));
            Console.WriteLine();

            Console.WriteLine("Kes = " + 1 / scale + " * ");
            writer.WriteToConsole(kes.Scale(scale));
            Console.WriteLine();

            Console.WriteLine("Kee = " + 1 / scale + " * ");
            writer.WriteToConsole(kee.Scale(scale));
            Console.WriteLine();
        }

        private void RoundMatrix(Matrix matrix) //TODO: add this as a Matrix method and in an interface returning a copy ofc
        {
            for (int r = 0; r < matrix.NumRows; ++r)
            {
                for (int c = 0; c < matrix.NumColumns; ++c)
                {
                    matrix[r, c] = 1e6 * Math.Round(matrix[r, c] * 1e-6, 3);
                }
            }
        }

        private void RoundList(IList<double> list)
        {
            for (int i = 0; i < list.Count; ++i) list[i] = Math.Round(list[i], 5);
        }
    }
}
