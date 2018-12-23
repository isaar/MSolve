using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class StandardPreconditionerCholesky : IStandardPreconditioner, IDisposable
    {
        private readonly CholeskySuiteSparse factorU;

        private StandardPreconditionerCholesky(DokSymmetric Kss)
        {
            int order = Kss.NumColumns;
            var (valuesStd, rowIndicesStd, colOffsetsStd) = Kss.BuildSymmetricCscArrays(true);
            Kss.Clear();  // No longer needed.
            factorU = CholeskySuiteSparse.Factorize(order, valuesStd.Length, valuesStd, rowIndicesStd, colOffsetsStd,
                true);
        }

        ~StandardPreconditionerCholesky()
        {
            ReleaseResources();
        }

        public void Dispose()
        {
            ReleaseResources();
            GC.SuppressFinalize(this);
        }

        public Vector PreconditionedMatrixTimesVector(Vector x)
        {
            // y = inv(Us^T) * Kss * inv(Us) * x = I * x
            return x.Copy();
        }

        public Vector PreconditionerTimesVector(Vector x, bool transposeThis)
        {
            if (transposeThis) // Ps^T * x = inv(Us^T) * x
            {
                return factorU.ForwardSubstitution(x);
            }
            else // Ps * x = inv(Us) * x
            {
                return factorU.BackSubstitution(x);
            }
        }

        private void ReleaseResources()
        {
            factorU.Dispose();
        }

        public class Builder: IStandardPreconditionerBuilder
        {
            private readonly StandardMatrixAssemblerCsc assembler;

            public Builder(Model2D model)
            {
                assembler = new StandardMatrixAssemblerCsc();
                //Ordering = new StandardOrderingNatural(); 
                Ordering = new StandardOrderingAmd(model);

                //TODO: find why the orderings give different results. Is it error built-up by so many back/forward solutions?
            }

            public IStandardMatrixAssembler Assembler { get { return assembler; } }

            public IStandardPreconditioner Build()
            {
                return new StandardPreconditionerCholesky(assembler.Kss);
            }

            public string Name { get; } = "Cholesky";

            public IStandardOrdering Ordering { get; }
        }
    }
}
