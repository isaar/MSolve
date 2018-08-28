using System;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Solvers.MenkBordas;
using ISAAR.MSolve.XFEM.Tests.Subdomains;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class TestMenkBordasSolver
    {
        private static readonly FullVectorWriter writer = new FullVectorWriter(false)
        {
            ArrayFormat = Array1DFormat.PlainVertical,
            NumericFormat = new GeneralNumericFormat()
        };

        public static void Run()
        {
            // Build the model
            (Model2D model, ISingleCrack crack) = SubdomainTest2.CreateModel();
            IDomainDecomposer decomposer = SubdomainTest2.DefinePartition(model);

            // Print matrices
            //(MenkBordasSystem sys, MenkBordasVector uExpected) = BuildCustomSystem(model, cluster);
            //PrintMatrices(sys);
            //PrintRhsVectors(sys, uExpected);

            // Solve with CG
            //SolveWithCG(model, cluster);

            // Solve with direct
            //SolveWithLU(model, cluster);

            // Solve actual system with the dedicated solver
            UseMenkBordasSolver(model, crack, decomposer);
            //UseSkylineSolver(model);
        }

        //public static (MenkBordasSystem sys, MenkBordasVector uExpected) BuildCustomSystem(Model2D model, XCluster2D cluster)
        //{
        //    MenkBordasSystem sys = BuildMatrices(model, cluster);
        //    MenkBordasVector uExpected = BuildLhs(sys);
        //    (MenkBordasMatrix K, Vector rhs) = sys.BuildSystem();
        //    MenkBordasVector f = K.MultiplyRight(uExpected);
        //    sys.bs = f.Vs;
        //    foreach (var subvector in f.Ve) sys.be.Add(subvector);
        //    sys.CheckDimensions();
        //    return (sys, uExpected);
        //}

        //public static MenkBordasVector BuildLhs(MenkBordasSystem sys)
        //{
        //    MenkBordasSystem.Dimensions dim = sys.CountDimensions();

        //    int seed = 2;
        //    double scale = 1000.0;
        //    var rand = new Random(seed);

        //    var xs = Vector.CreateZero(dim.NumDofsStd);
        //    for (int i = 0; i < xs.Length; ++i) xs[i] = scale * rand.NextDouble();

        //    var xe = new Vector[dim.NumSubdomains];
        //    for (int i = 0; i < dim.NumSubdomains; ++i)
        //    {
        //        xe[i] = Vector.CreateZero(dim.SubdomainEnds[i] - dim.SubdomainStarts[i]);
        //        for (int j = 0; j < xe[j].Length; ++j) xe[i][j] = scale * rand.NextDouble();
        //    }

        //    return new MenkBordasVector(sys.numSubdomains, 0, xs, xe, null);
        //}

        //public static MenkBordasSystem BuildMatrices(Model2D model, XCluster2D cluster)
        //{
        //    var sys = new MenkBordasSystem(cluster.Subdomains.Count);
        //    var assembler = new XClusterMatrixAssembler();
        //    (DOKSymmetricColMajor globalKss, DOKRowMajor globalKsc) = assembler.BuildStandardMatrices(model, cluster.DofOrderer);
        //    sys.Kss = globalKss;
        //    Dictionary<XSubdomain2D, SignedBooleanMatrix> booleanMatrices =
        //        assembler.BuildSubdomainSignedBooleanMatrices(cluster);
        //    foreach (var subdomain in cluster.Subdomains)
        //    {
        //        (DOKSymmetricColMajor Kee, DOKRowMajor Kes, DOKRowMajor Kec) =
        //            assembler.BuildSubdomainMatrices(subdomain, cluster.DofOrderer);
        //        sys.Kee.Add(Kee);
        //        var KesCSR = Kes.BuildCSRMatrix(true);
        //        sys.Kes.Add(KesCSR);
        //        sys.Kse.Add(KesCSR.TransposeToCSR());
        //        sys.B.Add(booleanMatrices[subdomain]);
        //    }
        //    return sys;
        //}

        //public static void PrintMatrices(MenkBordasSystem sys)
        //{
        //    Console.WriteLine("************************* Subdomain Matrices *************************");
        //    (MenkBordasMatrix K, Vector f) = sys.BuildSystem();
        //    K.WriteToConsole();
        //    Console.WriteLine();

        //    Console.WriteLine("************************* Global Matrix *************************");
        //    Matrix denseK = K.CopyToDense();
        //    FullMatrixWriter.NumericFormat = new GeneralNumericFormat();
        //    (new FullMatrixWriter(denseK)).WriteToConsole();
        //    Console.WriteLine();
        //}

        //public static void PrintRhsVectors(MenkBordasSystem sys, MenkBordasVector uExpected)
        //{
        //    (MenkBordasMatrix K, Vector denseF) = sys.BuildSystem();

        //    Console.WriteLine("************************* Global RHS Vector *************************");
        //    var formatting = Array1DFormat.PlainVertical;
        //    FullVectorWriter.NumericFormat = new GeneralNumericFormat();
        //    (new FullVectorWriter(denseF, false, formatting)).WriteToConsole();
        //    Console.WriteLine();

        //    Console.WriteLine("************************* Global RHS Vector computed only with dense *************************");
        //    Vector denseU = uExpected.CopyToDense();
        //    FullVectorWriter.NumericFormat = new GeneralNumericFormat();
        //    (new FullVectorWriter(K.CopyToDense() * denseU, false, formatting)).WriteToConsole();
        //    Console.WriteLine();
        //}

        public static void SolveWithCG(Model2D model, XCluster2D cluster)
        {
            //(MenkBordasSystem sys, MenkBordasVector uExpected) = BuildCustomSystem(model, cluster);
            //var cg = new MenkBordasCG(1000, 1e-10);
            //(MenkBordasVector u, IterativeStatistics stats) = cg.Solve(sys);
            //Console.WriteLine(stats);
            //MenkBordasVector diff = u.Axpy(-1, uExpected);
            //double error = Math.Sqrt(diff.DotProduct(diff)) / Math.Sqrt(uExpected.DotProduct(uExpected));
            //Console.WriteLine("Normalized error = " + error);
            //Console.WriteLine();
            //Console.WriteLine("Solution expected = ");
            //uExpected.WriteToConsole();
            //Console.WriteLine();
            //Console.WriteLine("Solution computed = ");
            //u.WriteToConsole();
        }

        //public static void SolveWithLU(Model2D model, XCluster2D cluster)
        //{
        //    (MenkBordasSystem sys, MenkBordasVector uExpected) = BuildCustomSystem(model, cluster);
        //    (MenkBordasMatrix K, Vector denseF) = sys.BuildSystem();
        //    Matrix denseK = K.CopyToDense();
        //    Vector denseUExpected = uExpected.CopyToDense();
        //    Vector denseU = denseK.FactorLU().SolveLinearSystem(denseF);

        //    var formatting = Array1DFormat.PlainVertical;
        //    FullVectorWriter.NumericFormat = new GeneralNumericFormat();
        //    Console.WriteLine("U expected (all dofs): ");
        //    (new FullVectorWriter(denseUExpected, false, formatting)).WriteToConsole();
        //    Console.WriteLine();
        //    Console.WriteLine("U computed (all dofs): ");
        //    (new FullVectorWriter(denseU, false, formatting)).WriteToConsole();
        //}

        public static void UseMenkBordasSolver(Model2D model, ISingleCrack crack, IDomainDecomposer decomposer)
        {
            int maxIterations = 1700;
            double tolerance = double.Epsilon;
            var solver = new MenkBordasSolver(model, crack, decomposer, maxIterations, tolerance, 
                new StandardPreconditionerCholesky.Builder(model), new EnrichedPreconditioningNaive());
            solver.Initialize();
            solver.Solve();

            //Console.WriteLine("U computed (all dofs): ");
            //(new FullVectorWriter(solver.Solution, false, formatting)).WriteToConsole();

            int elementID = 0;
            //TODO: isn't this computed in the solver as well?
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(solver.DofOrderer);
            foreach (var element in model.Elements)
            {
                Console.WriteLine($"Element {elementID++}: ");

                Vector stdU = solver.DofOrderer.ExtractDisplacementVectorOfElementFromGlobal(
                    element, solver.Solution, constrainedDisplacements);
                Console.WriteLine("Standard displacements: ");
                writer.WriteToConsole(stdU);
                Console.WriteLine();

                if (element.CountEnrichedDofs() > 0)
                {
                    Vector enrU = solver.DofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, solver.Solution);
                    Console.WriteLine("Enriched displacements: ");
                    writer.WriteToConsole(enrU);
                    Console.WriteLine();
                }
            }
        }

        public static void UseSkylineSolver(Model2D model)
        {
            var solver = new SkylineSolver(model);
            solver.Initialize();
            solver.Solve();

            //Console.WriteLine("U computed (all dofs): ");
            //(new FullVectorWriter(solver.Solution, false, formatting)).WriteToConsole();

            int elementID = 0;
            //TODO: isn't this computed in the solver as well?
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(solver.DofOrderer);
            foreach (var element in model.Elements)
            {
                Console.WriteLine($"Element {elementID++}: ");

                Vector stdU = solver.DofOrderer.ExtractDisplacementVectorOfElementFromGlobal(
                    element, solver.Solution, constrainedDisplacements);
                Console.WriteLine("Standard displacements: ");
                writer.WriteToConsole(stdU);
                Console.WriteLine();

                if (element.CountEnrichedDofs() > 0)
                {
                    Vector enrU = solver.DofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, solver.Solution);
                    Console.WriteLine("Enriched displacements: ");
                    writer.WriteToConsole(enrU);
                    Console.WriteLine();
                }
            }
        }
    }
}
