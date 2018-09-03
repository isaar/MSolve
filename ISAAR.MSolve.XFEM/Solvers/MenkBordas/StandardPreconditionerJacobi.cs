using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class StandardPreconditionerJacobi: IStandardPreconditioner
    {
        private readonly int order;
        private readonly CsrMatrix Kss;
        //private readonly Vector invDiagonal;
        private readonly Vector invDiagonalRoot;

        private StandardPreconditionerJacobi(DokRowMajor Kss)
        {
            CsrMatrix.UseMKL = true; //TODO: this should be done elsewhere
            this.Kss = Kss.BuildCsrMatrix(true);
            this.order = Kss.NumColumns;

            // Diagonal preconditioner
            (double[] diagonal, int firstZeroIdx) = Kss.GetDiagonalAsArray();
            if (firstZeroIdx != -1) throw new SingularMatrixException($"Zero diagonal entry at index {firstZeroIdx}");
            var invDiagonal = new double[order];
            var invDiagonalRoot = new double[order];
            for (int i = 0; i < order; ++i)
            {
                //if (Math.Abs(diagonal[i]) < tolerance) throw new SingularMatrixException($"Zero diagonal entry at index {i}");
                //invDiagonal[i] = 1.0 / diagonal[i];
                invDiagonal[i] = 1.0 / Math.Sqrt(diagonal[i]);
            }
            //this.invDiagonal = Vector.CreateFromArray(invDiagonal, false);
            this.invDiagonalRoot = Vector.CreateFromArray(invDiagonal, false);

            Kss.Clear();  // No longer needed.
        }

        public void Dispose()
        {
            // do nothing
        }

        public Vector PreconditionedMatrixTimesVector(Vector x)
        {
            // y = Ps^T * Kss * Ps * x = D * Kss * D * x
            Vector y = invDiagonalRoot.MultiplyEntrywise(x);
            y = Kss.MultiplyRight(y, false);
            return y = invDiagonalRoot.MultiplyEntrywise(y);

            //TODO: perhaps I can do D^2 * Kss * x
            //TODO: Use MKL's vector math for entrywise vector multiplications
        }

        public Vector PreconditionerTimesVector(Vector x, bool transposeThis)
        {
            var result = new double[order];
            return invDiagonalRoot.MultiplyEntrywise(x);
        }

        public class Builder : IStandardPreconditionerBuilder
        {
            private readonly StandardMatrixAssemblerCsr assembler;

            public Builder(Model2D model)
            {
                assembler = new StandardMatrixAssemblerCsr();
                Ordering = new StandardOrderingNatural();
                //Ordering = new StandardOrderingAmd(model); // Pretty much the same, at least for small problems.
            }

            public IStandardMatrixAssembler Assembler { get { return assembler; } }

            public IStandardPreconditioner Build()
            {
                return new StandardPreconditionerJacobi(assembler.Kss);
            }

            public string Name { get; } = "Jacobi";

            public IStandardOrdering Ordering { get; }
        }
    }
}
