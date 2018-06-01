using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    /// <summary>
    /// An object of this class is responsible for the necessary multiplications used in the Menk-Bordas preconditioning method.
    /// It is not responsible for managing its fields, which are injected during construction.
    /// </summary>
    class MenkBordasPrecondMatrix : ILinearTransformation<Vector>
    {
        private readonly MenkBordasSystem.Dimensions dim;
        private readonly IReadOnlyDictionary<XSubdomain2D, CSRMatrix> Kes;
        private readonly IReadOnlyDictionary<XSubdomain2D, CSRMatrix> Kse;
        private readonly IStandardPreconditioner Ps;
        private readonly IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> Pe;
        private readonly Matrix L;
        private readonly Matrix Q;

        public MenkBordasPrecondMatrix(MenkBordasSystem.Dimensions dimensions,
            IReadOnlyDictionary<XSubdomain2D, CSRMatrix> Kes, IReadOnlyDictionary<XSubdomain2D, CSRMatrix> Kse,
            IStandardPreconditioner Ps, IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> Pe, Matrix L, Matrix Q)
        {
            this.dim = dimensions;
            this.Kes = Kes;
            this.Kse = Kse;
            this.Ps = Ps;
            this.Pe = Pe;
            this.L = L;
            this.Q = Q;
        }

        public Vector Multiply(Vector x)
        {
            var y = Vector.CreateZero(dim.NumDofsAll);
            var xs = x.Slice(0, dim.NumDofsStd);

            // ys = Ps^T * Kss * Ps * xs
            Vector ys = Ps.PreconditionedMatrixTimesVector(xs);
            //Vector ys = prec.Ps.ForwardSubstitution(Kss.MultiplyRight(prec.Ps.BackSubstitution(xs))); // For other preconditioners

            if (Q != null) //TODO: something more explicit is needed
            {
                var xc = x.Slice(dim.EquationsStart, dim.NumDofsAll);

                // ye_all = Q^T * xc
                Vector yeAll = Q.MultiplyRight(xc, true);
                y.SetSubvector(yeAll, dim.NumDofsStd);

                // yc = Q * xe_all
                Vector yc = Q * x.Slice(dim.NumDofsStd, dim.EquationsStart); // Rows correspond to the continuity equations. TODO: these entries are sliced twice.
                y.SetSubvector(yc, dim.EquationsStart);
            }

            foreach (var sub in dim.Subdomains)
            {
                var xe = x.Slice(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]);

                // ys += Ps^T * Kse * inv(Ue) * xe
                ys.AddIntoThis(Ps.PreconditionerTimesVector(Kse[sub].MultiplyRight(Pe[sub].BackSubstitution(xe)), true));

                // ye = inv(Ue^T) * Kes * Ps * xs + I * xe
                Vector ye = Pe[sub].ForwardSubstitution(Kes[sub].MultiplyRight(Ps.PreconditionerTimesVector(xs, false)));
                ye.AddIntoThis(xe);
                y.AddSubvector(ye, dim.SubdomainStarts[sub]);
            }
            y.SetSubvector(ys, 0);
            return y;
        }

        public Vector PreconditionerTimesVector(Vector x, bool transposePreconditioner)
        {
            if (transposePreconditioner)
            {
                var y = Vector.CreateZero(dim.NumDofsAll);

                // ys = Ps^T * xs
                Vector ys = Ps.PreconditionerTimesVector(x.Slice(0, dim.NumDofsStd), true);
                y.SetSubvector(ys, 0);

                foreach (var sub in dim.Subdomains)
                {
                    // ye = Pe^T * xe = inv(Ue^T) * xe
                    Vector ye = Pe[sub].ForwardSubstitution(x.Slice(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]));
                    y.SetSubvector(ye, dim.SubdomainStarts[sub]);
                }

                if (L != null) //TODO: something more explicit is needed
                {
                    // yc = inv(L) * xc
                    Vector yc = L.Invert() * x.Slice(dim.EquationsStart, dim.NumDofsAll); //TODO: I MUST do optimizations here
                    y.SetSubvector(yc, dim.EquationsStart);
                }
                                           
                return y;
            }
            else
            {
                var y = Vector.CreateZero(dim.NumDofsAll);

                // ys = Ps * xs
                Vector ys = Ps.PreconditionerTimesVector(x.Slice(0, dim.NumDofsStd), false);
                y.SetSubvector(ys, 0);

                foreach (var sub in dim.Subdomains)
                {
                    // ye = Pe * xe = inv(Ue) * xe
                    Vector ye = Pe[sub].BackSubstitution(x.Slice(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]));
                    y.SetSubvector(ye, dim.SubdomainStarts[sub]);
                }

                if (L != null) //TODO: something more explicit is needed
                {
                    // yc = inv(L^T) * xc
                    Vector yc = L.Transpose().Invert() * x.Slice(dim.EquationsStart, dim.NumDofsAll); //TODO: I MUST do optimizations here
                    y.SetSubvector(yc, dim.EquationsStart);
                }

                return y;
            }
        }
    }
}
