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
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Output.VTK;

//TODO: remove checks and prints
//TODO: allow various preconditioners for Kss
//TODO: use AMD first. 
//TODO: Alternatively consider grouping all boundary dofs together. Perhaps this could reduce the QR effort to certain 
//      submatrices and even allow decoupling of Q and R. Ideally this ordering could coexist with AMD's by using a permutation
//      vector when multiplying with Q. Perhaps the permutation vector or the multiplication could also take advantage of empty
//      submatrices in Q.
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasSolver: ISolver //TODO: dispose of MenkBordasSystem
    {
        private readonly ICrackDescription crack;
        private readonly Model2D model;
        private readonly IDomainDecomposer decomposer;
        private readonly int maxIterations;
        private readonly double tolerance;
        private readonly string subdomainsDirectory;

        private XCluster2D cluster;
        private ISet<XSubdomain2D> enrichedSubdomains;
        private HashSet<XSubdomain2D> tipEnrichedSubdomains;
        private int iteration;
        private MenkBordasSystem system;

        /// <summary>
        /// Only entries related to constrained standard dofs are stored. These will not be modifed as the crack grows.
        /// </summary>
        private Vector Uc;

        public MenkBordasSolver(Model2D model, ICrackDescription crack, IDomainDecomposer decomposer, int maxIterations,
            double tolerance, string subdomainsDirectory = null)
        {
            this.model = model;
            this.crack = crack;
            this.decomposer = decomposer;
            this.maxIterations = maxIterations;
            this.tolerance = tolerance;
            this.subdomainsDirectory = subdomainsDirectory;
            Logger = new SolverLogger("MenkBordasSolver");
        }

        public IDofOrderer DofOrderer { get { return cluster.DofOrderer; } }

        public SolverLogger Logger { get; }

        public Vector Solution { get; private set; }

        /// <summary>
        /// Create and store anything that pertains to the standard dofs and will not change as the crack propagates. 
        /// </summary>
        public void Initialize() //TODO: I should also set up the domain decomposition, irregardless of current enrichments.
        {
            iteration = 0;
            tipEnrichedSubdomains = new HashSet<XSubdomain2D>();
            var watch = new Stopwatch();

            // Partion the domain into subdomains
            watch.Start();
            cluster = decomposer.CreateSubdomains();
            watch.Stop();
            Logger.LogDuration(iteration, "domain decomposition", watch.ElapsedMilliseconds);

            // Standard dofs are not divided into subdomains and will not change over time.
            watch.Restart();
            cluster.OrderStandardDofs(model);
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
            Vector Fs = model.CalculateStandardForces(DofOrderer);
            this.Uc = model.CalculateConstrainedDisplacements(DofOrderer); 
            Vector bs = Fs - globalKsc.MultiplyRight(Uc);
            watch.Stop();
            Logger.LogDuration(iteration, "linear system assembly", watch.ElapsedMilliseconds);

            // Preconditioner for Kss
            watch.Restart();
            this.system = new MenkBordasSystem(globalKss, bs);
            watch.Stop();
            Logger.LogDuration(iteration, "standard preconditioner", watch.ElapsedMilliseconds);
        }

        public void Solve() //TODO: only update the subdomains with at least 1 modified element
        {
            ++iteration;
            var watch = new Stopwatch();

            // Possibly update the domain decomposition
            watch.Start();
            decomposer.UpdateSubdomains(cluster);
            watch.Stop();
            Logger.LogDuration(iteration, "domain decomposition", watch.ElapsedMilliseconds);
            if (subdomainsDirectory != null) WriteDecomposition(subdomainsDirectory, cluster);

            watch.Restart();
            
            /// Order the dofs of subdomains that are modified.
            if (enrichedSubdomains == null) enrichedSubdomains = cluster.FindEnrichedSubdomains();
            HashSet<XSubdomain2D> previousTipSubdomains = tipEnrichedSubdomains;
            tipEnrichedSubdomains = FindTipEnrichedSubdomains();
            var modifiedSubdomains = new SortedSet<XSubdomain2D>();
            foreach (var subdomain in cluster.Subdomains)
            {
                if (tipEnrichedSubdomains.Contains(subdomain) || previousTipSubdomains.Contains(subdomain))
                {
                    modifiedSubdomains.Add(subdomain);
                    enrichedSubdomains.Add(subdomain);
                }
            }
            cluster.DofOrderer.OrderSubdomainDofs(enrichedSubdomains, modifiedSubdomains, crack);
            #region debug
            // These are inefficient, but useful for debugging
            //system.ClearSubdomains();
            //SortedSet<XSubdomain2D> modifiedSubdomains = cluster.FindEnrichedSubdomains();
            //cluster.DofOrderer.OrderSubdomainDofs(modifiedSubdomains, modifiedSubdomains, crack);

            //Console.Write("Enriched subdomains: ");
            //foreach (var subdomain in enrichedSubdomains) Console.Write(subdomain.ID + " ");
            //Console.WriteLine();
            //if (enrichedSubdomains.Count > 2)
            //{
            //    Console.WriteLine();
            //}
            #endregion

            /// Subdomain enriched matrices and rhs
            var assembler = new XClusterMatrixAssembler();
            foreach (var subdomain in modifiedSubdomains)
            {
                (DOKSymmetricColMajor Kee, DOKRowMajor Kes, DOKRowMajor Kec) =
                    assembler.BuildSubdomainMatrices(subdomain, cluster.DofOrderer);
                var KesCSR = Kes.BuildCSRMatrix(true);
                var be = Kec.MultiplyRight(Uc.Scale(-1.0), true); //Fe = 0 - Kec * Uc.
                system.SetSubdomainMatrices(subdomain, Kee, KesCSR, KesCSR.TransposeToCSR(), be);
            }

            /// Signed boolean matrices and continuity equations
            var booleanMatrices = assembler.BuildSubdomainSignedBooleanMatrices(cluster);
            if (booleanMatrices.Count > 1) system.SetBooleanMatrices(booleanMatrices);
            #region debug
            //foreach (var subdomainB in booleanMatrices)
            //{
            //    Console.WriteLine("Subdomain " + subdomainB.Key.ID);
            //    (new FullMatrixWriter(subdomainB.Value)).WriteToConsole();
            //    Console.WriteLine();
            //}
            #endregion
            watch.Stop();
            Logger.LogDuration(iteration, "linear system assembly", watch.ElapsedMilliseconds);

            /// Enriched preconditioners 
            watch.Restart();
            (MenkBordasPrecondMatrix precMatrix, Vector rhs) = system.BuildPreconditionedSystem();
            watch.Stop();
            Logger.LogDuration(iteration, "enriched preconditioner", watch.ElapsedMilliseconds);

            /// Use an iterative algorithm to solve the symmetric indefinite system
            watch.Restart();
            var minres = new MinimumResidual(maxIterations, tolerance, 0, false, false);
            Vector precRhs = precMatrix.PreconditionerTimesVector(rhs, true);
            (Vector precSolution, MinresStatistics stats) = minres.Solve(precMatrix, precRhs);
            Solution = precMatrix.PreconditionerTimesVector(precSolution, false);
            watch.Stop();
            Logger.LogDuration(iteration, "iterative algorithm", watch.ElapsedMilliseconds);
            Console.WriteLine(stats);

            Logger.LogDofs(iteration, DofOrderer.NumStandardDofs + DofOrderer.NumEnrichedDofs);
            #region Debug
            //CheckMultiplication(matrix, rhs);
            //CheckPrecondMultiplication(sys);
            //Vector xExpected = SolveLUandPrintMatrices(matrix, rhs);
            //CompareSolutionWithLU(sys, Solution);
            #endregion
        }

        private HashSet<XSubdomain2D> FindTipEnrichedSubdomains()
        {
            var tipSubdomains = new HashSet<XSubdomain2D>();
            IEnumerable<ISet<XNode2D>> allTipNodes = crack.CrackTipNodesNew.Values;
            foreach (var subdomain in cluster.Subdomains)
            {
                if (IsEnrichedSubdomain(subdomain, allTipNodes)) tipSubdomains.Add(subdomain);
            }
            return tipSubdomains;
        }

        private bool IsEnrichedSubdomain(XSubdomain2D subdomain, IEnumerable<ISet<XNode2D>> allTipNodes)
        {
            foreach (var tipNodes in allTipNodes)
            {
                foreach (var node in tipNodes)
                {
                    if (subdomain.AllNodes.Contains(node)) return true;
                }
            }
            return false;
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

        private void WriteDecomposition(string subdomainsDirectory, XCluster2D cluster)
        {
            var writer = new DomainDecompositionWriter();

            //writer.WriteRegions(directory + "regions.vtk", regions);
            writer.WriteSubdomainElements(subdomainsDirectory + $"\\subdomains_{iteration-1}.vtk", cluster.Subdomains);
            writer.WriteBoundaryNodes(subdomainsDirectory + $"\\boundaryNodes_{iteration-1}.vtk", cluster.Subdomains);
        }
    }
}
