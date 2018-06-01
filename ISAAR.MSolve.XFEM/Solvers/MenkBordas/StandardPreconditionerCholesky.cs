using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class StandardPreconditionerCholesky : IStandardPreconditioner, IDisposable
    {
        private readonly CholeskySuiteSparse factorU;

        private StandardPreconditionerCholesky(DOKSymmetricColMajor Kss)
        {
            int order = Kss.NumColumns;
            var (valuesStd, rowIndicesStd, colOffsetsStd) = Kss.BuildSymmetricCSCArrays(true);
            Kss.Clear();  // No longer needed.
            factorU = CholeskySuiteSparse.Factorize(order, valuesStd.Length, valuesStd, rowIndicesStd, colOffsetsStd,
                true, SuiteSparseOrdering.Natural);
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
            public IStandardPreconditioner Build(DOKSymmetricColMajor Kss)
            {
                return new StandardPreconditionerCholesky(Kss);
            }
        }
    }
}
