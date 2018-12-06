using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    /// <summary>
    /// An object of this class is responsible for the necessary multiplications used in the Menk-Bordas preconditioning method.
    /// It is not responsible for managing its fields, which are injected during construction.
    /// </summary>
    class MenkBordasPrecondMatrix : ILinearTransformation
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

        public int NumColumns => dim.NumDofsAll;

        public int NumRows => dim.NumDofsAll;

        public void Multiply(IVectorView x, IVector y)
        {
            y = Vector.CreateZero(dim.NumDofsAll);
            var xs = ((Vector)x).GetSubvector(0, dim.NumDofsStd);

            // ys = Ps^T * Kss * Ps * xs
            Vector ys = Ps.PreconditionedMatrixTimesVector(xs);

            if (LQ != null) //TODO: something more explicit is needed
            {
                var xc = ((Vector)x).GetSubvector(dim.EquationsStart, dim.NumDofsAll);

                // ye_all = Q^T * xc
                Vector yeAll = LQ.QTimesVector(xc, true);
                y.CopySubvectorFrom(dim.NumDofsStd, yeAll, 0, yeAll.Length);

                // yc = Q * xe_all. Rows correspond to the continuity equations. 
                //TODO: these entries are sliced twice.
                Vector yc = LQ.QTimesVector(((Vector)x).GetSubvector(dim.NumDofsStd, dim.EquationsStart), false); 
                y.CopySubvectorFrom(dim.EquationsStart, yc, 0, yc.Length);
            }

            foreach (var sub in dim.Subdomains)
            {
                var xe = ((Vector)x).GetSubvector(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]);

                // ys += Ps^T * Kse * inv(Ue) * xe
                ys.AddIntoThis(Ps.PreconditionerTimesVector(Kse[sub].Multiply(Pe[sub].BackSubstitution(xe)), true));

                // ye = inv(Ue^T) * Kes * Ps * xs + I * xe
                Vector ye = Pe[sub].ForwardSubstitution(Kes[sub].Multiply(Ps.PreconditionerTimesVector(xs, false)));
                ye.AddIntoThis(xe);
                //y.AddSubvector(ye, dim.SubdomainStarts[sub]);
                y.AddSubvectorIntoThis(dim.SubdomainStarts[sub], ye, 0, ye.Length);
            }
            y.CopySubvectorFrom(0, ys, 0, ys.Length);
        }

        public IVector PreconditionerTimesVector(IVector x, bool transposePreconditioner)
        {
            var y = Vector.CreateZero(dim.NumDofsAll);

            // ys = Ps^T * xs
            Vector ys = Ps.PreconditionerTimesVector(((Vector)x).GetSubvector(0, dim.NumDofsStd), transposePreconditioner);
            y.CopySubvectorFrom(0, ys, 0, ys.Length);

            if (LQ != null) //TODO: something more explicit is needed
            {
                // yc = inv(L) * xc
                Vector yc = LQ.InverseLTimesVector(
                    ((Vector)x).GetSubvector(dim.EquationsStart, dim.NumDofsAll), transposePreconditioner);
                y.CopySubvectorFrom(dim.EquationsStart, yc, 0, yc.Length);
            }

            if (transposePreconditioner)
            {
                foreach (var sub in dim.Subdomains)
                {
                    // ye = Pe^T * xe = inv(Ue^T) * xe
                    Vector ye = Pe[sub].ForwardSubstitution(
                        ((Vector)x).GetSubvector(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]));
                    y.CopySubvectorFrom(dim.SubdomainStarts[sub], ye, 0, ye.Length);
                }
            }
            else
            {
                foreach (var sub in dim.Subdomains)
                {
                    // ye = Pe * xe = inv(Ue) * xe
                    Vector ye = Pe[sub].BackSubstitution(
                        ((Vector)x).GetSubvector(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]));
                    y.CopySubvectorFrom(dim.SubdomainStarts[sub], ye, 0, ye.Length);
                }
            }

            return y;
        }
    }
}
