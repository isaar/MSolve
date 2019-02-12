using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Optimization.Logging;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP
{
    /// <summary>
    /// C# translation of the Matlab code published in "A 99 line topology optimization code written in Matlab, O. Sigmund, 1991".
    /// I tried to keep the code structure and variable names as close to the original as possible, with the exception of the
    /// extensions presented in the paper, which required some structural changes.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class TopologySimp99Lines
    {
        public enum BoundaryConditions
        {
            MbbBeam, ShortCantilever, Cantilever2LoadCases
        }

        public enum PassiveElements
        {
            No, HoleInCantilever
        }

        public enum OptimAlgorithm
        {
            OC, MMA
        }

        // Number of elements along the x, y axes
        private readonly int nelx, nely;

        /// <summary>The prescribed volume fraction.</summary>
        private readonly double volfrac;

        /// <summary>The exponent of the penalty function.</summary>
        private readonly double penal;

        /// <summary>The radius outside which the convolution operator of the filer is zero.</summary>
        private readonly double rmin;

        private readonly Func<Matrix, IList<Vector>> feAnalysis;
        private readonly Func<Matrix, bool[,]> applyPassiveElements;
        private readonly Func<Matrix, Matrix, bool[,], Matrix> optimAlgorithm;

        public TopologySimp99Lines(int nelx, int nely, double volfrac, double penal, double rmin, 
            BoundaryConditions bc, PassiveElements passive, OptimAlgorithm alg)
        {
            this.nelx = nelx;
            this.nely = nely;
            this.volfrac = volfrac;
            this.penal = penal;
            this.rmin = rmin;

            if (bc == BoundaryConditions.MbbBeam) feAnalysis = FEAnalysisHalfMbbBeam;
            else if (bc == BoundaryConditions.ShortCantilever) feAnalysis = FEAnalysisCantilever;
            else if (bc == BoundaryConditions.Cantilever2LoadCases) feAnalysis = FEAnalysisCantilever2LoadCases;

            if (passive == PassiveElements.No) applyPassiveElements = NoPassiveElements;
            else if (passive == PassiveElements.HoleInCantilever) applyPassiveElements = PassiveElementsForHoleInCantilever;

            if (alg == OptimAlgorithm.OC) optimAlgorithm = OptimalityCriteriaUpdate;
            else if (alg == OptimAlgorithm.MMA) optimAlgorithm = MethodMovingAsymptotesUpdate;
        }

        public (double compliance, Matrix densities, ObjectiveFunctionLogger logger) TopologyOptimization()
        {
            // Initialize
            var x = Matrix.CreateWithValue(nely, nelx, volfrac); // design variables: element densities
            bool[,] passive = applyPassiveElements(x); // extension for passive elements
            var dc = Matrix.CreateZero(nely, nelx); // compliance sensitivities
            int loop = 0;
            double change = 1.0;
            double c = double.MaxValue;
            var logger = new ObjectiveFunctionLogger();

            // Start iterations
            while (change > 0.01)
            {
                ++loop;
                Matrix xold = x.Copy();

                // FE analysis
                IList<Vector> U = feAnalysis(x);

                // Objective function and sensitivity analysis
                Matrix Ke = ElementStiffness(); // common for all elements
                c = 0.0;
                for (int ely = 1; ely <= nely; ++ely)
                {
                    for (int elx = 1; elx <= nelx; ++elx)
                    {
                        int n1 = (nely + 1) * (elx - 1) + ely;
                        int n2 = (nely + 1) * elx + ely;
                        dc[ely - 1, elx - 1] = 0.0;
                        int[] elementDofs = GetElementDofs(elx, ely);
                        foreach (Vector Ui in U)
                        {
                            Vector Ue = Ui.GetSubvector(elementDofs);
                            double unitCompliance = Ue * (Ke * Ue);
                            double xe = x[ely - 1, elx - 1];
                            c += Math.Pow(xe, penal) * unitCompliance;
                            dc[ely - 1, elx - 1] -= Math.Pow(xe, penal - 1) * unitCompliance * penal;
                        }
                    }
                }

                // Filtering of sensitivities
                dc = MeshIndependencyFilter(x, dc);

                // Design update by the optimality criteria method
                x = OptimalityCriteriaUpdate(x, dc, passive);
                change = (x - xold).MaxAbsolute();

                // Print results
                Console.WriteLine($"Iteration = {loop}, Compliance = {c}, Volume = {x.Sum() / (nelx*nely)}, Change = {change}");
                //Console.WriteLine("Densities:");
                //var writer = new FullMatrixWriter();
                //writer.WriteToConsole(x);
                //Console.WriteLine("-----------------------------------");
                logger.Log(c);
            }

            return (c, x, logger);
        }

        private Matrix OptimalityCriteriaUpdate(Matrix x, Matrix dc, bool[,] passive)
        {
            double l1 = 0.0;
            double l2 = 1E5;
            double move = 0.2;
            double xmin = 0.001;
            double xmax = 1.0;
            Matrix xnew = null;
            while (l2-l1 > 1E-4)
            {
                double lmid = 0.5 * (l1 + l2);

                // Update densities
                xnew = x.DoEntrywise(dc, (xe, dce) =>
                {
                    double xeBe = xe * Math.Sqrt(-dce / lmid); // η = 1/2 => Be^η = Sqrt(Be)
                    double xMinus = Math.Max(xmin, xe - move);
                    double xPlus = Math.Min(xmax, xe + move);
                    if (xeBe <= xMinus) return xMinus;
                    else if (xeBe < xPlus) return xeBe;
                    else return xPlus;
                });

                // Extension for passive elements
                for (int ely = 0; ely < nely; ++ely)
                {
                    for (int elx = 0; elx < nelx; ++elx)
                    {
                        if (passive[ely, elx]) xnew[ely, elx] = 0.001;
                    }
                }

                // Bi-sectioning
                double constraint = xnew.Sum() - volfrac * nelx * nely ;
                if (constraint > 0) l1 = lmid;
                else l2 = lmid;
            }
            return xnew;
        }

        private Matrix MeshIndependencyFilter(Matrix x, Matrix dc)
        {
            var dcn = Matrix.CreateZero(nely, nelx);
            int rminRounded = (int)Math.Round(rmin);
            for (int i = 1; i <= nelx; ++i)
            {
                for (int j = 1; j <= nely; ++j)
                {
                    double sum = 0.0;
                    for (int k = Math.Max(i - rminRounded, 1); k <= Math.Min(i + rminRounded, nelx); ++k)
                    {
                        for (int m = Math.Max(j - rminRounded, 1); m <= Math.Min(j + rminRounded, nely); ++m)
                        {
                            double fac = rmin - Math.Sqrt((i - k) * (i - k) + (j - m) * (j - m));
                            if (fac > 0)
                            {
                                sum += fac;
                                dcn[j - 1, i - 1] += fac * x[m - 1, k - 1] * dc[m - 1, k - 1];
                            }
                        }
                    }
                    dcn[j - 1, i - 1] /= x[j - 1, i - 1] * sum;
                }
            }
            return dcn;
        }

        /// <summary>
        /// Element stiffness matrix
        /// </summary>
        private static Matrix ElementStiffness()
        {
            double E = 1.0; // Young's modulus
            double nu = 0.3; // Poisson ratio
            double[] k = { 0.5 - nu / 6.0, 0.125 + nu / 8.0, -0.25 - nu / 12.0, -0.125 + 3 * nu / 8.0,
                -0.25 + nu / 12.0, -0.125 - nu / 8.0, nu / 6.0, 0.125 - 3 * nu / 8.0 }; // unique stiffness matrix entries
            var Ke = Matrix.CreateFromArray(new double[,]
            {
                 { k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7] },
                 { k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2] },
                 { k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1] },
                 { k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4] },
                 { k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3] },
                 { k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6] },
                 { k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5] },
                 { k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0] }
            });
            Ke.ScaleIntoThis(E / (1 - nu * nu));
            return Ke;
        }

        private IList<Vector> FEAnalysisHalfMbbBeam(Matrix x)
        {
            int numAllDofs = 2 * (nelx + 1) * (nely + 1);
            var K = DokSymmetric.CreateEmpty(numAllDofs);
            var F = Vector.CreateZero(numAllDofs);
            var U = Vector.CreateZero(numAllDofs);

            Matrix Ke = ElementStiffness();
            int numElementDofs = Ke.NumRows;
            int[] elementDofsLocal = Enumerable.Range(0, numElementDofs).ToArray();

            // Global stiffness matrix assembly
            for (int ely = 1; ely <= nely; ++ely)
            {
                for (int elx = 1; elx <= nelx; ++elx)
                {
                    int[] elementDofsGlobal = GetElementDofs(elx, ely);
                    double density = Math.Pow(x[ely - 1, elx - 1], penal);
                    K.AddSubmatrixSymmetric(density * Ke, elementDofsLocal, elementDofsGlobal);
                }
            }

            // Define loads and supports for half MBB beam
            F[1] = -1.0;
            var fixedDofs = new HashSet<int>(); //TODO: Use LINQ to simplify this
            for (int i = 0; i < 2 * (nely + 1); i += 2) fixedDofs.Add(i);
            fixedDofs.Add(numAllDofs - 1);
            int[] freeDofs = Enumerable.Range(0, numAllDofs).Except(fixedDofs).ToArray();

            // Solve linear system
            DokSymmetric Kf = K.GetSubmatrix(freeDofs);
            Vector Ff = F.GetSubvector(freeDofs);
            Vector Uf = CholeskyCSparseNet.Factorize(Kf.BuildSymmetricCscMatrix(true)).SolveLinearSystem(Ff);
            for (int i = 0; i < freeDofs.Length; ++i) U[freeDofs[i]] = Uf[i];
            //U.CopyNonContiguouslyFrom(freeDofs, Uf, Enumerable.Range(0, freeDofs.Length).ToArray()); // alternative way, but probably slower.
            //foreach (int i in fixedDofs) U[i] = 0.0; // They are 0 by default

            return new Vector[] { U };
        }

        private IList<Vector> FEAnalysisCantilever(Matrix x)
        {
            int numAllDofs = 2 * (nelx + 1) * (nely + 1);
            var K = DokSymmetric.CreateEmpty(numAllDofs);
            var F = Vector.CreateZero(numAllDofs);
            var U = Vector.CreateZero(numAllDofs);

            Matrix Ke = ElementStiffness();
            int numElementDofs = Ke.NumRows;
            int[] elementDofsLocal = Enumerable.Range(0, numElementDofs).ToArray();

            // Global stiffness matrix assembly
            for (int ely = 1; ely <= nely; ++ely)
            {
                for (int elx = 1; elx <= nelx; ++elx)
                {
                    int[] elementDofsGlobal = GetElementDofs(elx, ely);
                    double density = Math.Pow(x[ely - 1, elx - 1], penal);
                    K.AddSubmatrixSymmetric(density * Ke, elementDofsLocal, elementDofsGlobal);
                }
            }

            // Define loads and supports for cantilever beam
            F[numAllDofs - 1] = -1.0;
            IEnumerable<int> fixedDofs = Enumerable.Range(0, 2 * (nely + 1));
            int[] freeDofs = Enumerable.Range(0, numAllDofs).Except(fixedDofs).ToArray();

            // Solve linear system
            DokSymmetric Kf = K.GetSubmatrix(freeDofs);
            Vector Ff = F.GetSubvector(freeDofs);
            Vector Uf = CholeskyCSparseNet.Factorize(Kf.BuildSymmetricCscMatrix(true)).SolveLinearSystem(Ff);
            for (int i = 0; i < freeDofs.Length; ++i) U[freeDofs[i]] = Uf[i];
            //U.CopyNonContiguouslyFrom(freeDofs, Uf, Enumerable.Range(0, freeDofs.Length).ToArray()); // alternative way, but probably slower.
            //foreach (int i in fixedDofs) U[i] = 0.0; // They are 0 by default

            return new Vector[] { U };
        }

        private IList<Vector> FEAnalysisCantilever2LoadCases(Matrix x)
        {
            int numLoadCases = 2;
            int numAllDofs = 2 * (nelx + 1) * (nely + 1);
            var K = DokSymmetric.CreateEmpty(numAllDofs);
            var F = new Vector[] { Vector.CreateZero(numAllDofs), Vector.CreateZero(numAllDofs) };
            var U = new Vector[] { Vector.CreateZero(numAllDofs), Vector.CreateZero(numAllDofs) };

            Matrix Ke = ElementStiffness();
            int numElementDofs = Ke.NumRows;
            int[] elementDofsLocal = Enumerable.Range(0, numElementDofs).ToArray();

            // Global stiffness matrix assembly
            for (int ely = 1; ely <= nely; ++ely)
            {
                for (int elx = 1; elx <= nelx; ++elx)
                {
                    int[] elementDofsGlobal = GetElementDofs(elx, ely);
                    double density = Math.Pow(x[ely - 1, elx - 1], penal);
                    K.AddSubmatrixSymmetric(density * Ke, elementDofsLocal, elementDofsGlobal);
                }
            }

            // Define supports for cantilever beam
            IEnumerable<int> fixedDofs = Enumerable.Range(0, 2 * (nely + 1));
            int[] freeDofs = Enumerable.Range(0, numAllDofs).Except(fixedDofs).ToArray();

            // Load case 1: unit load towards -y at bottom right corner
            F[0][2 * (nelx + 1) * (nely + 1) - 1] = -1.0;

            // Load case 2: unit load towards +y at top right corner
            F[1][2 * (nelx) * (nely + 1) + 1] = 1.0;

            // Solve linear system
            CholeskyCSparseNet factorizedKf = CholeskyCSparseNet.Factorize(
                K.GetSubmatrix(freeDofs).BuildSymmetricCscMatrix(true)); // only factorize the matrix once
            for (int c = 0; c < numLoadCases; ++c)
            {
                Vector Ff = F[c].GetSubvector(freeDofs);
                Vector Uf = factorizedKf.SolveLinearSystem(Ff);
                for (int i = 0; i < freeDofs.Length; ++i) U[c][freeDofs[i]] = Uf[i];
                //U[c].CopyNonContiguouslyFrom(freeDofs, Uf, Enumerable.Range(0, freeDofs.Length).ToArray()); // alternative way, but probably slower.
                //foreach (int i in fixedDofs) U[c][i] = 0.0; // They are 0 by default
            }
            return U;
        }

        private bool[,] NoPassiveElements(Matrix x) => new bool[nely, nelx];

        private bool[,] PassiveElementsForHoleInCantilever(Matrix x)
        {
            var passive = new bool[nely, nelx];
            for (int ely = 1; ely <= nely; ++ely)
            {
                for (int elx = 1; elx <= nelx; ++elx)
                {
                    double distance = Math.Sqrt(Math.Pow(ely - nely / 2.0, 2) + Math.Pow(elx - nelx / 3.0, 2));
                    if (distance < nely / 3.0)
                    {
                        passive[ely - 1, elx - 1] = true;
                        x[ely - 1, elx - 1] = 0.001;
                    }
                }
            }
            return passive;
        }

        /// <summary>
        /// Simplified MMA by Bendsoe. Unfortunately it converges slower than Optimality Criteria and oscillates.
        /// </summary>
        private Matrix MethodMovingAsymptotesUpdate(Matrix x, Matrix dc, bool[,] passive)
        {
            double xlow = 0.001, xhigh = 1.0;
            Matrix L = x - 0.1 * (xhigh - xlow) * Matrix.CreateWithValue(nely, nelx, 1.0);
            Matrix xMinusLSquared = (x - L).Square();
            Matrix high = xMinusLSquared.MultiplyEntrywise(
                dc.DoEntrywise(L, (dci, Li) => -dci / Math.Pow(xlow - Li, 2)));
            Matrix low = xMinusLSquared.MultiplyEntrywise(
                dc.DoEntrywise(L, (dci, Li) => -dci / Math.Pow(xhigh - Li, 2)));
            double l2 = high.Min();
            double l1 = low.Max();

            Matrix xnew = null;
            for (int i = 0; i < 50; ++i) // 50 is totally arbitrary
            {
                double lmid = 0.5 * (l1 + l2);
                Matrix dcRoot = dc.DoToAllEntries(dci => Math.Sqrt(-dci / lmid));
                xnew = L + x.DoEntrywise(L, (xi, Li) => Math.Abs(xi - Li)).MultiplyEntrywise(dcRoot);
                xnew.DoToAllEntriesIntoThis(xi => Math.Max(xlow, xi));

                // Extension for passive elements
                for (int ely = 0; ely < nely; ++ely)
                {
                    for (int elx = 0; elx < nelx; ++elx)
                    {
                        if (passive[ely, elx]) xnew[ely, elx] = 0.001;
                    }
                }

                if (xnew.Sum() - volfrac * nelx * nely > 0) l2 = lmid;
                else l1 = lmid;
            }

            return xnew;
        }

        /// <summary>
        /// The global indices of the element's dofs, in 0-based numbering.
        /// </summary>
        /// <param name="elx">1-based numbering</param>
        /// <param name="ely">1-based numbering</param>
        private int[] GetElementDofs(int elx, int ely)
        {
            int n1 = (nely + 1) * (elx - 1) + ely;
            int n2 = (nely + 1) * elx + ely;
            return new int[]
            {
                2 * n1 - 2,     2 * n1 - 1,     2 * n2 - 2,     2 * n2 - 1,
                2 * n2,         2 * n2 + 1,     2 * n1,         2 * n1 + 1
            };
        }
    }
}
