using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP
{
    public static class TopologySimp99Lines
    {
        public static void TopologyOptimization(int nelx, int nely, double volfrac, double penal, double rmin)
        {
            // Initialize
            var x = Matrix.CreateWithValue(nely, nelx, volfrac); // design variables: element densities
            var dc = Matrix.CreateZero(nely, nelx); // compliance sensitivities
            int loop = 0;
            double change = 1.0;
            var writer = new FullMatrixWriter();
            double c = double.MaxValue;

            // Start iterations
            while (change > 0.01)
            {
                ++loop;
                Matrix xold = x.Copy();

                // FE analysis
                Vector U = FEAnalysis(nelx, nely, x, penal);

                // Objective function and sensitivity analysis
                Matrix Ke = ElementStiffness(); // common for all elements
                c = 0.0;
                for (int ely = 1; ely <= nely; ++ely)
                {
                    for (int elx = 1; elx <= nelx; ++elx)
                    {
                        int n1 = (nely + 1) * (elx - 1) + ely;
                        int n2 = (nely + 1) * elx + ely;
                        int[] elementDofs = GetElementDofs(nelx, nely, elx, ely);
                        Vector Ue = U.GetSubvector(elementDofs);
                        double unitCompliance = Ue * (Ke * Ue);
                        double xe = x[ely - 1, elx - 1];
                        c += Math.Pow(xe, penal) * unitCompliance;
                        dc[ely - 1, elx - 1] = - Math.Pow(xe, penal - 1) * unitCompliance * penal;
                    }
                }

                // Filtering of sensitivities
                dc = MeshIndependencyFilter(nelx, nely, rmin, x, dc);

                // Design update by the optimality criteria method
                x = OptimalityCriteriaUpdate(nelx, nely, x, volfrac, dc);

                // Print results
                change = (x - xold).MaxAbsolute();
                Console.WriteLine($"Iteration = {loop}, Compliance = {c}, Volume = {x.Sum() / (nelx*nely)}, Change = {change}");

                // Print densities
                Console.WriteLine("Densities:");
                writer.WriteToConsole(x);
                Console.WriteLine("-----------------------------------");
            }
        }

        private static Matrix OptimalityCriteriaUpdate(int nelx, int nely, Matrix x, double volfrac, Matrix dc)
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
                
                // Bi-sectioning
                if (xnew.Sum() - volfrac * nelx * nely > 0) l1 = lmid;
                else l2 = lmid;
            }
            return xnew;
        }

        private static Vector FEAnalysis(int nelx, int nely, Matrix x, double penal)
        {
            //TODO: Use sparse matrix DOK -> Skyline or DOK -> CSC and use CSparse.NET

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
                    int[] elementDofsGlobal = GetElementDofs(nelx, nely, elx, ely);
                    double density = Math.Pow(x[ely - 1, elx - 1], penal);
                    K.AddSubmatrixSymmetric(density * Ke, elementDofsLocal, elementDofsGlobal);
                    //for (int i = 0; i < numElementDofs; ++i)
                    //{
                    //    for (int j = 0; j < numElementDofs; ++j)
                    //    {
                    //        K.AddToEntry(elementDofsGlobal[i], elementDofsGlobal[j], density * Ke[i, j]); 
                    //    }
                    //}
                }
            }

            // Define loads and supports (half MBB-beam)
            F[1] = -1;
            var fixedDofs = new HashSet<int>(); //TODO: Use LINQ to simplify this madness
            for (int i = 0; i < 2 * (nely + 1); i += 2) fixedDofs.Add(i);
            fixedDofs.Add(numAllDofs - 1);
            int[] freeDofs = Enumerable.Range(0, numAllDofs).Except(fixedDofs).ToArray();
            DokSymmetric Kf = K.GetSubmatrix(freeDofs);
            Vector Ff = F.GetSubvector(freeDofs);
            Vector Uf = CholeskyCSparseNet.Factorize(Kf.BuildSymmetricCscMatrix(true)).SolveLinearSystem(Ff);
            for (int i = 0; i < freeDofs.Length; ++i) U[freeDofs[i]] = Uf[i];
            U.CopyNonContiguouslyFrom(freeDofs, Uf, Enumerable.Range(0, freeDofs.Length).ToArray()); // alternative way, but probably slower.
            foreach (int i in fixedDofs) U[i] = 0.0; // They are 0 by default

            return U;
        }

        private static Matrix MeshIndependencyFilter(int nelx, int nely, double rmin, Matrix x, Matrix dc)
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
                                dcn[j-1, i-1] += fac * x[m-1, k-1] * dc[m-1, k-1];
                            }
                        }
                    }
                    dcn[j-1, i-1] /= x[j-1, i-1] * sum; 
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

        /// <summary>
        /// The global indices of the element's dofs, in 0-based numbering.
        /// </summary>
        /// <param name="nelx"></param>
        /// <param name="nely"></param>
        /// <param name="elx">1-based numbering</param>
        /// <param name="ely">1-based numbering</param>
        private static int[] GetElementDofs(int nelx, int nely, int elx, int ely)
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
