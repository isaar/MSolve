using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: remove checks and prints
//TODO: allow various preconditioners for Kss
//TODO: use AMD first. 
//TODO: Alternatively consider grouping all boundary dofs together. Perhaps this could reduce the QR effort to certain 
//      submatrices and even allow decoupling of Q and R. Ideally this ordering could coexist with AMD's by using a permutation
//      vector when multiplying with Q. Perhaps the permutation vector or the multiplication could also take advantage of empty
//      submatrices in Q.
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasSolver: ISolver
    {
        private readonly Model2D model;
        private readonly XCluster2D cluster;
        private readonly int maxIterations;
        private readonly double tolerance;

        private MenkBordasSystem system;

        /// <summary>
        /// Only entries related to constrained standard dofs are stored. These will not be modifed as the crack grows.
        /// </summary>
        private Vector Uc;

        public MenkBordasSolver(Model2D model, XCluster2D cluster, int maxIterations, double tolerance)
        {
            this.model = model;
            this.cluster = cluster;
            this.maxIterations = maxIterations;
            this.tolerance = tolerance;
            Logger = new SolverLogger();
        }

        public IDofOrderer DofOrderer { get { return cluster.DofOrderer; } }

        public SolverLogger Logger { get; }

        public Vector Solution { get; private set; }

        /// <summary>
        /// Create and store anything that pertains to the standard dofs and will not change as the crack propagates. 
        /// </summary>
        public void Initialize() //TODO: I should also set up the domain decomposition, irregardless of current enrichments.
        {
            var watch = new Stopwatch();
            watch.Start();

            // Standard dofs are not divided into subdomains and will not change over time.
            cluster.OrderDofs(model);
            int numStdDofs = cluster.DofOrderer.NumStandardDofs;
            var assembler = new XClusterMatrixAssembler();
            (DOKSymmetricColMajor globalKss, DOKRowMajor globalKsc) = assembler.BuildStandardMatrices(model, cluster.DofOrderer);
            //if (!Kss.IsSymmetric(1e-10)) throw new AsymmetricMatrixException(
            //    "Stiffness matrix corresponding to std-std dofs is not symmetric");

            /* 
             * The extended linear system is:
             * [Kcc Kcu; Kuc Kuu] * [uc; uu] = [Fc; Fu]
             * where c are the standard constrained dofs, f are the standard free dofs, e are the enriched dofs and 
             * u = Union(f,c) are both the dofs with unknown left hand side vectors: uu = [uf; ue].
             * To solve the system (for the unknowns ul):
             * i) Kuu * uu = Fu - Kuc * uc = Feff
             * ii) uu = Kuu \ Feff 
             */
            Vector globalFu = model.CalculateFreeForces(DofOrderer); //TODO: wasted space on enriched dofs. They will always be 0.
            this.Uc = model.CalculateConstrainedDisplacements(DofOrderer); 
            Vector Fu = globalFu.Slice(0, numStdDofs);
            Vector Fs = Fu - globalKsc.MultiplyRight(Uc);

            this.system = new MenkBordasSystem(globalKss, Fs);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        public void Solve() //TODO: only update the subdomains with at least 1 modified element
        {
            var watch = new Stopwatch();
            watch.Start();

            //int numSubdomains = cluster.Subdomains.Count;
            // TODO: should I order the dofs (cluster.OrderDofs(model)) again?

            /// Create the linear system
            var assembler = new XClusterMatrixAssembler();

            /// Signed boolean matrices and continuity equations
            system.SetBooleanMatrices(assembler.BuildSubdomainSignedBooleanMatrices(cluster));

            /// Subdomain enriched matrices and rhs
            foreach (var subdomain in cluster.Subdomains)
            {
                (DOKSymmetricColMajor Kee, DOKRowMajor Kes, DOKRowMajor Kec) =
                    assembler.BuildSubdomainMatrices(subdomain, cluster.DofOrderer);
                var KesCSR = Kes.BuildCSRMatrix(true);
                var be = Kec.MultiplyRight(Uc.Scale(-1.0), true); //Fe = 0 - Kec * Uc.
                system.SetSubdomainMatrices(subdomain, Kee, KesCSR, KesCSR.TransposeToCSR(), be);
            }
            //sys.CheckDimensions();

            /// Use an iterative algorithm to solve the symmetric indefinite system
            var minres = new MinimumResidual(maxIterations, tolerance, 0, false, false);

            /// With preconditioning
            (MenkBordasPrecondMatrix precMatrix, Vector rhs) = system.BuildPreconditionedSystem();
            Vector precRhs = precMatrix.PreconditionerTimesVector(rhs, true);
            (Vector precSolution, MinresStatistics stats) = minres.Solve(precMatrix, precRhs);
            Solution = precMatrix.PreconditionerTimesVector(precSolution, false);
            Console.WriteLine(stats);

            #region Debug
            //CheckMultiplication(matrix, rhs);
            //CheckPrecondMultiplication(sys);
            //Vector xExpected = SolveLUandPrintMatrices(matrix, rhs);
            //CompareSolutionWithLU(sys, Solution);
            #endregion

            /// Find the solution vector without multiple dofs
            // TODO:Should I do that? Isn't the ClusterDofOrderer responsible for extracting the correct dofs, when needed?

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        private static void CheckMultiplication(MenkBordasMatrix K, Vector b)
        {
            var rand = new Random();
            var x = Vector.CreateZero(b.Length);
            for (int i = 0; i < x.Length; ++i) x[i] = rand.NextDouble();
            Matrix denseK = K.CopyToDense();

            Vector yExpected = denseK * x;
            Vector yComputed = K.Multiply(x);
            double error = (yComputed - yExpected).Norm2() / yExpected.Norm2();
        }

        //private static void CheckPrecondMultiplication(MenkBordasSystem sys)
        //{
        //    // Build the matrices
        //    (MenkBordasMatrix K, Vector fCopy) = sys.BuildSystem();
        //    (MenkBordasPrecondMatrix Kbar, Vector f) = sys.BuildPreconditionedSystem();
        //    Matrix denseK = K.CopyToDense();

        //    // Build the lhs
        //    var rand = new Random();
        //    var xBar = Vector.CreateZero(f.Length);
        //    for (int i = 0; i < f.Length; ++i)
        //    {
        //        xBar[i] = rand.NextDouble();
        //        //xBar[i] = 1.0;
        //    }

        //    // Solve rigged system
        //    Vector x = Kbar.PreconditionerTimesVector(xBar, false);
        //    Vector yExpected = denseK * x;
        //    Vector yBarExpected = Kbar.PreconditionerTimesVector(yExpected, true);
        //    Vector yBar = Kbar.Multiply(xBar);

        //    double error = (yBar - yBarExpected).Norm2() / yBarExpected.Norm2();
        //    Console.WriteLine("Normalized error in preconditioned multiplication = " + error);
        //}

        //private static void CompareSolutionWithLU(MenkBordasSystem sys, Vector solution)
        //{
        //    (MenkBordasMatrix K, Vector f) = sys.BuildSystem();
        //    Matrix denseK = K.CopyToDense();
        //    Vector xExpected = denseK.FactorLU().SolveLinearSystem(f);
        //    double error = (solution - xExpected).Norm2() / xExpected.Norm2();
        //    Console.WriteLine("Normalized difference |xMINRES - xLU| / |xLU| = " + error);
        //}

        private static Vector SolveLUandPrintMatrices(MenkBordasMatrix K, Vector f)
        {
            Matrix denseK = K.CopyToDense();
            Vector denseX = denseK.FactorLU().SolveLinearSystem(f);

            // Sparse global matrix
            int order = denseK.NumRows;
            var sparseK = DOKColMajor.CreateEmpty(order, order);
            for (int j = 0; j < order; ++j)
            {
                for (int i = 0; i < order; ++i)
                {
                    double val = denseK[i, j];
                    if (val != 0.0) sparseK[i, j] = val;
                }
            }

            // Export the global system to Matlab
            var writer = new MatlabWriter();
            string directory = @"C:\Users\Serafeim\Desktop\GRACM\MenkBordas_debugging\";
            writer.WriteSparseMatrix(sparseK, directory + "global_matrix.txt");
            writer.WriteFullVector(f, directory + "global_rhs.txt");
            writer.WriteFullVector(denseX, directory + "solution.txt");

            // Export the submatrices to Matlab
            K.WriteToFiles(directory);

            return denseX;
        }
    }
}
