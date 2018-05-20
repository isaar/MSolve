using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Allow it to update some submatrices, while keeping others the same. 
//TODO: Break up Q and L into submatrices or store the submatrices that assemble the intermediate matrix before QR factorization.
//TODO: I need a new DOK class:
//      1) The same DOK should produce the full CSR (for multiplications) and the upper CSC (for factorization).
//          Fortunately the arrays for the upper CSC are identical to the ones for lower CSR.
//      2) Sorting the DOK should only be done once. It actually costs substantially. Also a lot of other operations can only 
//          once, when creating the full and lower CSR.
//      3) Ok having the full CSR/CSC matrix is good for multiplications. However can't I factorize using just that? Do I need
//          a separate upper copy? I think SuiteSparse can ignore the lower/upper entries
//TODO: Q and L are dense now. Try to take advantage of any sparsity. Perhaps create some by appropriate DD or permutations.
//      Also try to take advantage of the fact that B are boolean matrices, with at most 2 entries per column
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasPreconditioner: IDisposable
    {
        public readonly int numSubdomains;
        public readonly int numEquations;
        public readonly int numDofsAll;
        public readonly int numDofsStd;
        public readonly int[] subdomainStarts;
        public readonly int[] subdomainEnds;
        public readonly int equationsStart;

        public readonly CholeskySuiteSparse Ps; //TODO: abstract this so other preconditioners can be used.
        public readonly CholeskySuiteSparse[] Pe;
        public readonly Matrix Q;
        public readonly Matrix L;

        private MenkBordasPreconditioner(int numSubdomains, int numEquations, int numDofsStd, int numDofsAll,
            int[] subdomainStarts, int[] subdomainEnds, int equationsStart,
            CholeskySuiteSparse Ps, CholeskySuiteSparse[] Pe, Matrix Q, Matrix L)
        {
            this.numSubdomains = numSubdomains;
            this.numEquations = numEquations;
            this.numDofsStd = numDofsStd;
            this.numDofsAll = numDofsAll;
            this.subdomainStarts = subdomainStarts;
            this.subdomainEnds = subdomainEnds;
            this.equationsStart = equationsStart;

            this.Ps = Ps;
            this.Pe = Pe;
            this.Q = Q;
            this.L = L;
        }

        public static MenkBordasPreconditioner Create(int numSubdomains, int numEquations, int numDofsStd, int numDofsAll,
            int[] subdomainStarts, int[] subdomainEnds, int equationsStart,
            DOKSymmetricColMajor Kss, IReadOnlyList<DOKSymmetricColMajor> Kee, IReadOnlyList<SignedBooleanMatrix> B)
        {
            // Do not use "using(...){...}". We want the unmanaged memory to persist.
            var (valuesStd, rowIndicesStd, colOffsetsStd) = Kss.BuildSymmetricCSCArrays(true);

            // Standard preconditioner
            var Ps = CholeskySuiteSparse.Factorize(Kss.NumRows, valuesStd.Length, valuesStd, rowIndicesStd, colOffsetsStd,
                true, SuiteSparseOrdering.Natural);

            // Matrix that will undergo QR:
            // (B*Pe)^T = [B1*inv(U1) B2*inv(U2) ...]^T = [inv(L1)*B1^T inv(L2)*B2^T ...]
            // Dims: B = numEquations -by- numDofsEnr, Pe = numDofsEnr -by- numDofsEnr, (B*Pe)^T = numDofsEnr -by- numEquations
            int numDofsEnr = numDofsAll - numDofsStd - numEquations;
            var BPeTransp = Matrix.CreateZero(numDofsEnr, numEquations);

            var Pe = new CholeskySuiteSparse[numSubdomains];
            for (int i = 0; i < numSubdomains; ++i)
            {
                // Enriched preconditioner = cholesky factor U
                var (valuesEnr, rowIndicesEnr, colOffsetsEnr) = Kee[i].BuildSymmetricCSCArrays(true);
                Pe[i] = CholeskySuiteSparse.Factorize(Kee[i].NumRows, valuesEnr.Length, valuesEnr, rowIndicesEnr, colOffsetsEnr,
                    true, SuiteSparseOrdering.Natural);

                // Contribution to the matrix that will undergo QR
                Matrix contribution = Pe[i].ForwardSubstitution(B[i].CopyToFullMatrix(true));
                BPeTransp.SetSubmatrix(subdomainStarts[i] - numDofsStd, 0, contribution);
            }

            // LQ factorization 
            //TODO: various optimizations might be possible here
            var qr = BPeTransp.FactorQR();
            Matrix L = qr.GetFactorR().Slice(0, numEquations, 0, numEquations).Transpose(); //TODO: should probably use a packed UpperTriangular R
            Matrix Q = qr.GetFactorQ().Slice(0, numDofsEnr, 0, numEquations).Transpose(); //TODO: MKL has routines that only build some columns of Q!!!

            #region Debug
            //The problem is with Q or with its multiplications (L is ok though). Print Q and Q1 and check it against matlab

            //FullMatrixWriter.NumericFormat = new ExponentialFormat { NumDecimalDigits = 4 };
            //Console.WriteLine("Before QR: Pe^T * B^T = ");
            //(new FullMatrixWriter(BPeTransp)).WriteToConsole();
            //Console.WriteLine();

            //Console.WriteLine("Q = ");
            //(new FullMatrixWriter(qr.GetFactorQ().Slice(0, numDofsEnr, 0, numEquations))).WriteToConsole();
            //Console.WriteLine();

            //Console.WriteLine("R = ");
            //(new FullMatrixWriter(qr.GetFactorR().Slice(0, numEquations, 0, numEquations))).WriteToConsole();
            //Console.WriteLine();
            #endregion

            return new MenkBordasPreconditioner(numSubdomains, numEquations, numDofsStd, numDofsAll,
                subdomainStarts, subdomainEnds, equationsStart,
                Ps, Pe, Q, L);
        }

        public void Dispose()
        {
            if (Ps != null) Ps.Dispose();
            for (int i = 0; i < numSubdomains; ++i)
            {
                if (Pe[i] != null) Pe[i].Dispose();
            }
        }

        public Vector Multiply(Vector x, bool transpose)
        {
            if (transpose)
            {
                var y = Vector.CreateZero(numDofsAll);

                // ys = Ps * xs = inv(Us^T) * xs
                Vector ys = Ps.ForwardSubstitution(x.Slice(0, numDofsStd));
                y.SetSubvector(ys, 0);

                for (int i = 0; i < numSubdomains; ++i)
                {
                    // ye = Pe * xe = inv(Ue^T) * xe
                    Vector ye = Pe[i].ForwardSubstitution(x.Slice(subdomainStarts[i], subdomainEnds[i]));
                    y.SetSubvector(ye, subdomainStarts[i]);
                }

                // yc = inv(L^T) * xc
                Vector yc = L.Invert() * x.Slice(equationsStart, numDofsAll); //TODO: I MUST do optimizations here
                y.SetSubvector(yc, equationsStart);
                return y;
            }
            else
            {
                var y = Vector.CreateZero(numDofsAll);

                // ys = Ps * xs = inv(Us) * xs
                Vector ys = Ps.BackSubstitution(x.Slice(0, numDofsStd));
                y.SetSubvector(ys, 0);

                for (int i = 0; i < numSubdomains; ++i)
                {
                    // ye = Pe * xe = inv(Ue) * xe
                    Vector ye = Pe[i].BackSubstitution(x.Slice(subdomainStarts[i], subdomainEnds[i]));
                    y.SetSubvector(ye, subdomainStarts[i]);
                }

                // yc = inv(L^T) * xc
                Vector yc = L.Transpose().Invert() * x.Slice(equationsStart, numDofsAll); //TODO: I MUST do optimizations here
                y.SetSubvector(yc, equationsStart);
                return y;
            }
        }
    }
}
