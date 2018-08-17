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
        private readonly IReadOnlyDictionary<XSubdomain2D, CsrMatrix> Kes;
        private readonly IReadOnlyDictionary<XSubdomain2D, CsrMatrix> Kse;
        private readonly IStandardPreconditioner Ps;
        private readonly IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> Pe;
        private readonly IFactorizationLQ LQ;

        public MenkBordasPrecondMatrix(MenkBordasSystem.Dimensions dimensions,
            IReadOnlyDictionary<XSubdomain2D, CsrMatrix> Kes, IReadOnlyDictionary<XSubdomain2D, CsrMatrix> Kse,
            IStandardPreconditioner Ps, IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> Pe, IFactorizationLQ LQ)
        {
            this.dim = dimensions;
            this.Kes = Kes;
            this.Kse = Kse;
            this.Ps = Ps;
            this.Pe = Pe;
            this.LQ = LQ;
        }

        public Vector Multiply(Vector x)
        {
            var y = Vector.CreateZero(dim.NumDofsAll);
            var xs = x.GetSubvector(0, dim.NumDofsStd);

            // ys = Ps^T * Kss * Ps * xs
            Vector ys = Ps.PreconditionedMatrixTimesVector(xs);

            if (LQ != null) //TODO: something more explicit is needed
            {
                var xc = x.GetSubvector(dim.EquationsStart, dim.NumDofsAll);

                // ye_all = Q^T * xc
                Vector yeAll = LQ.QTimesVector(xc, true);
                y.SetSubvector(yeAll, dim.NumDofsStd);

                // yc = Q * xe_all. Rows correspond to the continuity equations. 
                //TODO: these entries are sliced twice.
                Vector yc = LQ.QTimesVector(x.GetSubvector(dim.NumDofsStd, dim.EquationsStart), false); 
                y.SetSubvector(yc, dim.EquationsStart);
            }

            foreach (var sub in dim.Subdomains)
            {
                var xe = x.GetSubvector(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]);

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
            var y = Vector.CreateZero(dim.NumDofsAll);

            // ys = Ps^T * xs
            Vector ys = Ps.PreconditionerTimesVector(x.GetSubvector(0, dim.NumDofsStd), transposePreconditioner);
            y.SetSubvector(ys, 0);

            if (LQ != null) //TODO: something more explicit is needed
            {
                // yc = inv(L) * xc
                Vector yc = LQ.InverseLTimesVector(x.GetSubvector(dim.EquationsStart, dim.NumDofsAll), transposePreconditioner);
                y.SetSubvector(yc, dim.EquationsStart);
            }

            if (transposePreconditioner)
            {
                foreach (var sub in dim.Subdomains)
                {
                    // ye = Pe^T * xe = inv(Ue^T) * xe
                    Vector ye = Pe[sub].ForwardSubstitution(x.GetSubvector(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]));
                    y.SetSubvector(ye, dim.SubdomainStarts[sub]);
                }
            }
            else
            {
                foreach (var sub in dim.Subdomains)
                {
                    // ye = Pe * xe = inv(Ue) * xe
                    Vector ye = Pe[sub].BackSubstitution(x.GetSubvector(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]));
                    y.SetSubvector(ye, dim.SubdomainStarts[sub]);
                }
            }

            return y;
        }
    }
}
