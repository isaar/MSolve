using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasPrecondMatrix : ILinearTransformation<Vector>
    {
        public readonly int numSubdomains;
        public readonly int numEquations;
        public readonly int numDofsAll;
        public readonly int numDofsStd;
        public readonly int[] subdomainStarts;
        public readonly int[] subdomainEnds;
        public readonly int equationsStart;

        public readonly CSRMatrix Kss;
        public readonly CSRMatrix[] Kee;
        public readonly CSRMatrix[] Kes;
        public readonly CSRMatrix[] Kse;
        public readonly SignedBooleanMatrix[] B;

        public readonly MenkBordasPreconditioner prec;

        public MenkBordasPrecondMatrix(int numSubdomains, int numEquations, int numDofsStd, int numDofsAll,
            int[] subdomainStarts, int[] subdomainEnds, int equationsStart,
            CSRMatrix Kss, CSRMatrix[] Kee, CSRMatrix[] Kes, CSRMatrix[] Kse, SignedBooleanMatrix[] B,
            MenkBordasPreconditioner preconditioner)
        {
            this.numSubdomains = numSubdomains;
            this.numEquations = numEquations;
            this.numDofsStd = numDofsStd;
            this.numDofsAll = numDofsAll;
            this.subdomainStarts = subdomainStarts;
            this.subdomainEnds = subdomainEnds;
            this.equationsStart = equationsStart;

            this.Kss = Kss;
            this.Kee = Kee;
            this.Kes = Kes;
            this.Kse = Kse;
            this.B = B;

            this.prec = preconditioner;
        }

        public Vector Multiply(Vector x)
        {
            var y = Vector.CreateZero(numDofsAll);
            var xs = x.Slice(0, numDofsStd);
            var xc = x.Slice(equationsStart, numDofsAll);

            // ys = inv(Us^T) * Kss * inv(Us) * xs
            Vector ys = xs.Copy(); // For cholesky preconditioner
            //Vector ys = prec.Ps.ForwardSubstitution(Kss.MultiplyRight(prec.Ps.BackSubstitution(xs))); // For other preconditioners

            // ye_all = Q^T * xc
            Vector yeAll = prec.Q.MultiplyRight(xc, true);
            y.SetSubvector(yeAll, numDofsStd);

            // yc = Q * xe_all
            Vector yc = prec.Q * x.Slice(numDofsStd, equationsStart); // Rows correspond to the continuity equations. TODO: these entries are sliced twice.
            y.SetSubvector(yc, equationsStart);

            for (int i = 0; i < numSubdomains; ++i)
            {
                var xe = x.Slice(subdomainStarts[i], subdomainEnds[i]);

                // ys += inv(Us^T) * Kse * inv(Ue) * xe
                ys.AddIntoThis(prec.Ps.ForwardSubstitution(Kse[i].MultiplyRight(prec.Pe[i].BackSubstitution(xe))));

                // ye = inv(Ue^T) * Kes * inv(Us) * xs + I * xe
                Vector ye = prec.Pe[i].ForwardSubstitution(Kes[i].MultiplyRight(prec.Ps.BackSubstitution(xs)));
                ye.AddIntoThis(xe);
                y.AddSubvector(ye, subdomainStarts[i]);
            }
            y.SetSubvector(ys, 0);
            return y;
        }
    }
}
