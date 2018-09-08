using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    interface IEnrichedPreconditioning
    {
        IEnrichedOrdering Ordering { get; }

        IFactorizationLQ CreateContinuityEquationsPreconditioner(MenkBordasSystem.Dimensions dimensions,
             IReadOnlyDictionary<XSubdomain2D, SignedBooleanMatrix> B,
             IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> Pe);

        CholeskySuiteSparse CreateEnrichedPreconditioner(DokSymmetric Kee);
    }

    interface IFactorizationLQ
    {
        Vector InverseLTimesVector(Vector x, bool transposePreconditioner);
        Vector QTimesVector(Vector x, bool transposeQ);
    }
}
