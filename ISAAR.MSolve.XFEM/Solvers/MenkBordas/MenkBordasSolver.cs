using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: remove checks and prints
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasSolver: ISolver
    {
        private readonly Model2D model;
        private readonly XCluster2D cluster;
        private readonly int maxIterations;
        private readonly double tolerance;

        /// <summary>
        /// Only entries related to free standard dofs are stored. These will not be modifed as the crack grows.
        /// </summary>
        private CSRMatrix Kss;

        /// <summary>
        /// Only entries related to free standard dofs are stored. These will not be modifed as the crack grows.
        /// </summary>
        private Vector Fs;

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
        public void Initialize()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Standard dofs are not divided into subdomains and will not change over time.
            cluster.OrderDofs(model);
            int numStdDofs = cluster.DofOrderer.NumStandardDofs;
            var assembler = new XClusterMatrixAssembler();
            (DOKRowMajor globalKss, DOKRowMajor globalKsc) = assembler.BuildStandardMatrices(model, cluster.DofOrderer);
            this.Kss = globalKss.BuildCSRMatrix(true);
            if (!Kss.IsSymmetric(1e-10)) throw new AsymmetricMatrixException(
                "Stiffness matrix corresponding to std-std dofs is not symmetric");

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
            this.Fs = Fu - globalKsc.MultiplyRight(Uc);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        public void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            int numSubdomains = cluster.Subdomains.Count;
            // TODO: should I order the dofs (cluster.OrderDofs(model)) again?

            // Create the linear system
            var assembler = new XClusterMatrixAssembler();
            var sys = new MenkBordasSystem(cluster.Subdomains.Count);
            sys.Kss = Kss;
            sys.bs = Fs;

            // Signed boolean matrices and continuity equations
            Dictionary<XSubdomain2D, SignedBooleanMatrix> booleanMatrices =
                assembler.BuildSubdomainSignedBooleanMatrices(cluster);
            sys.numContinuityEquations = booleanMatrices[cluster.Subdomains[0]].NumRows;
            sys.bc = Vector.CreateZero(sys.numContinuityEquations);

            // Subdomain enriched matrices and rhs
            foreach (var subdomain in cluster.Subdomains)
            {
                (DOKSymmetricColMajor Kee, DOKRowMajor Kes, DOKRowMajor Kec) =
                    assembler.BuildSubdomainMatrices(subdomain, cluster.DofOrderer);

                sys.Kee.Add(DOKRowMajor.CreateFromSparseMatrix(Kee).BuildCSRMatrix(true)); // Not the most efficient, but ok for testing
                if (!Kee.IsSymmetric()) throw new AsymmetricMatrixException(
                    "Stiffness matrix corresponding to enr-enr dofs of a subdomain is not symmetric");
                var KesCSR = Kes.BuildCSRMatrix(true);
                sys.Kes.Add(KesCSR);
                sys.Kse.Add(KesCSR.TransposeToCSR());
                sys.B.Add(booleanMatrices[subdomain]);

                //Fe = 0 - Kec * Uc. TODO: do this without building the CSR matrix
                sys.be.Add(Kec.MultiplyRight(Uc.Scale(-1.0), true)); 
            }
            sys.CheckDimensions();

            // Use an iterative algorithm to solve the system
            //var cg = new MenkBordasCG(maxIterations, tolerance); // CG is not suited for this indefinite matrix
            //(MenkBordasVector u, IterativeStatistics stats) = cg.Solve(sys);
            var minres = new MenkBordasMinres(maxIterations, tolerance);
            (MenkBordasVector u, MinresStatistics stats) = minres.Solve(sys);
            Solution = u.CopyToDense();
            Console.WriteLine(stats);

            #region Debug
            Vector xExpected = SolveLUandPrintMatrices(sys);
            #endregion

            // Find the solution vector without multiple dofs
            // TODO:Should I do that? Isn't the ClusterDofOrderer responsible for extracting the correct dofs, when needed?

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        private static Vector SolveLUandPrintMatrices(MenkBordasSystem sys)
        {
            // Use a dense direct algorithm to solve the system (for testing)
            MenkBordasMatrix K = sys.BuildMatrix();
            MenkBordasVector F = sys.BuildRhsVector();

            Matrix denseK = K.CopyToDense();
            Vector denseF = F.CopyToDense();
            Vector denseX = denseK.FactorLU().SolveLinearSystem(denseF);

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
            writer.WriteFullVector(denseF, directory + "global_rhs.txt");
            writer.WriteFullVector(denseX, directory + "solution.txt");

            // Export the submatrices to Matlab
            writer.WriteSparseMatrix((ISparseMatrix)(K.Kss), directory + "Kss.txt");
            writer.WriteFullVector(F.Vs, directory + "fs.txt");
            for (int i = 0; i < sys.numSubdomains; ++i)
            {
                writer.WriteSparseMatrix((ISparseMatrix)(K.Kee[i]), directory + $"Kee{i+1}.txt");
                writer.WriteSparseMatrix((ISparseMatrix)(K.Kes[i]), directory + $"Kes{i+1}.txt");
                writer.WriteSparseMatrix((ISparseMatrix)(K.B[i]), directory + $"B{i+1}.txt");
                writer.WriteFullVector(F.Ve[i], directory + $"fe{i+1}.txt");
            }
            writer.WriteFullVector(F.Vc, directory + "fc.txt"); // This vector is 0.

            return denseX;
        }
    }
}
