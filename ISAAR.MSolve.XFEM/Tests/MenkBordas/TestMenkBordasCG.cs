using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Solvers.Algorithms;
using ISAAR.MSolve.XFEM.Solvers.Algorithms.MenkBordas;

namespace ISAAR.MSolve.XFEM.Tests.MenkBordas
{
    class TestMenkBordasCG
    {
        public static void Run()
        {
            // Build the model
            Model2D model = SubdomainTest2.CreateModel();
            XCluster2D cluster = SubdomainTest2.CreateSubdomains(model);
            cluster.OrderDofs(model);

            // Print matrices
            //(MenkBordasSystem sys, MenkBordasVector uExpected) = BuildCustomSystem(model, cluster);
            //PrintMatrices(sys);
            //PrintRhsVectors(sys, uExpected);

            // Solve with CG
            SolveWithCG(model, cluster);

            // Solve with direct
            //SolveWithLU(model, cluster);
        }

        public static (MenkBordasSystem sys, MenkBordasVector uExpected) BuildCustomSystem(Model2D model, XCluster2D cluster)
        {
            MenkBordasSystem sys = BuildMatrices(model, cluster);
            MenkBordasVector uExpected = BuildLhs(sys);
            MenkBordasMatrix K = sys.BuildMatrix();
            MenkBordasVector f = K.MultiplyRight(uExpected);
            sys.bs = f.Vs;
            foreach (var subvector in f.Ve) sys.be.Add(subvector);
            sys.bc = f.Vc;
            sys.CheckDimensions();
            return (sys, uExpected);
        }

        public static MenkBordasVector BuildLhs(MenkBordasSystem sys)
        {
            int seed = 2;
            double scale = 1000.0;
            var rand = new Random(seed);

            var xs = Vector.CreateZero(sys.Kss.NumColumns);
            for (int i = 0; i < xs.Length; ++i) xs[i] = scale * rand.NextDouble();

            var xe = new Vector[sys.Kee.Count];
            for (int sub = 0; sub < sys.Kee.Count; ++sub)
            {
                xe[sub] = Vector.CreateZero(sys.Kee[sub].NumColumns);
                for (int i = 0; i < xe[sub].Length; ++i) xe[sub][i] = scale * rand.NextDouble();
            }

            var xc = Vector.CreateZero(sys.numContinuityEquations); // This should be 0 right?

            return new MenkBordasVector(sys.numSubdomains, sys.numContinuityEquations, xs, xe, xc);
        }

        public static MenkBordasSystem BuildMatrices(Model2D model, XCluster2D cluster)
        {
            var sys = new MenkBordasSystem(cluster.Subdomains.Count);
            var assembler = new XClusterMatrixAssembler();
            (DOKRowMajor globalKss, DOKRowMajor globalKsc) = assembler.BuildStandardMatrices(model, cluster.DofOrderer);
            sys.Kss = globalKss.BuildCSRMatrix(true);
            Dictionary<XSubdomain2D, SignedBooleanMatrix> booleanMatrices =
                assembler.BuildSubdomainSignedBooleanMatrices(cluster);
            sys.numContinuityEquations = booleanMatrices[cluster.Subdomains[0]].NumRows;
            foreach (var subdomain in cluster.Subdomains)
            {
                (DOKSymmetricColMajor Kee, DOKRowMajor Kes, DOKRowMajor Kec) =
                    assembler.BuildSubdomainMatrices(subdomain, cluster.DofOrderer);
                sys.Kee.Add(DOKRowMajor.CreateFromSparseMatrix(Kee).BuildCSRMatrix(true)); // Not the most efficient, but ok for testing
                var KesCSR = Kes.BuildCSRMatrix(true);
                sys.Kes.Add(KesCSR);
                sys.Kse.Add(KesCSR.TransposeToCSR());
                sys.B.Add(booleanMatrices[subdomain]);
            }
            return sys;
        }

        public static void PrintMatrices(MenkBordasSystem sys)
        {
            Console.WriteLine("************************* Subdomain Matrices *************************");
            MenkBordasMatrix K = sys.BuildMatrix();
            K.WriteToConsole();
            Console.WriteLine();

            Console.WriteLine("************************* Global Matrix *************************");
            Matrix denseK = K.CopyToDense();
            FullMatrixWriter.NumericFormat = new GeneralNumericFormat();
            (new FullMatrixWriter(denseK)).WriteToConsole();
            Console.WriteLine();
        }

        public static void PrintRhsVectors(MenkBordasSystem sys, MenkBordasVector uExpected)
        {
            Console.WriteLine("************************* Subdomain RHS Vectors *************************");
            MenkBordasVector f = sys.BuildRhsVector();
            f.WriteToConsole();
            Console.WriteLine();

            Console.WriteLine("************************* Global RHS Vector *************************");
            var formatting = new Array1DFormatting("", "", "\n");
            Vector denseF = f.CopyToDense();
            FullVectorWriter.NumericFormat = new GeneralNumericFormat();
            (new FullVectorWriter(denseF, false, formatting)).WriteToConsole();
            Console.WriteLine();

            Console.WriteLine("************************* Global RHS Vector computed only with dense *************************");
            Vector denseU = uExpected.CopyToDense();
            FullVectorWriter.NumericFormat = new GeneralNumericFormat();
            (new FullVectorWriter(sys.BuildMatrix().CopyToDense() * denseU, false, formatting)).WriteToConsole();
            Console.WriteLine();
        }

        public static void SolveWithCG(Model2D model, XCluster2D cluster)
        {
            (MenkBordasSystem sys, MenkBordasVector uExpected) = BuildCustomSystem(model, cluster);
            var cg = new MenkBordasCG(1000, 1e-10);
            (MenkBordasVector u, IterativeStatistics stats) = cg.Solve(sys);
            Console.WriteLine(stats);
            MenkBordasVector diff = u.Axpy(-1, uExpected);
            double error = Math.Sqrt(diff.DotProduct(diff)) / Math.Sqrt(uExpected.DotProduct(uExpected));
            Console.WriteLine("Normalized error = " + error);
            Console.WriteLine();
            Console.WriteLine("Solution expected = ");
            uExpected.WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("Solution computed = ");
            u.WriteToConsole();
        }

        public static void SolveWithLU(Model2D model, XCluster2D cluster)
        {
            (MenkBordasSystem sys, MenkBordasVector uExpected) = BuildCustomSystem(model, cluster);
            MenkBordasMatrix K = sys.BuildMatrix();
            MenkBordasVector f = sys.BuildRhsVector();
            Matrix denseK = K.CopyToDense();
            Vector denseF = f.CopyToDense();
            Vector denseUExpected = uExpected.CopyToDense();
            Vector denseU = denseK.FactorLU().SolveLinearSystem(denseF);

            var formatting = new Array1DFormatting("", "", "\n");
            FullVectorWriter.NumericFormat = new GeneralNumericFormat();
            Console.WriteLine("U expected (all dofs): ");
            (new FullVectorWriter(denseUExpected, false, formatting)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("U computed (all dofs): ");
            (new FullVectorWriter(denseU, false, formatting)).WriteToConsole();
        }
    }
}
