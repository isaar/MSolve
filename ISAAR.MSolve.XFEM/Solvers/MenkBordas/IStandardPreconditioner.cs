using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    interface IStandardPreconditioner
    {
        Vector PreconditionedMatrixTimesVector(Vector x);
        Vector PreconditionerTimesVector(Vector x, bool transposeThis);
    }

    interface IStandardPreconditionerBuilder
    {
        IStandardOrdering Ordering { get; }
        IStandardPreconditioner Build(DOKSymmetricColMajor Kss);
    }
}
