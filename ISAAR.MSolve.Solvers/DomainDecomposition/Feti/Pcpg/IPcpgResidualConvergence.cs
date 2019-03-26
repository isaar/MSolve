using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Pcpg
{
    internal interface IPcpgResidualConvergence
    {
        /// <summary>
        /// Calculates the ratio norm2(f(r(x))) / norm2(g(r(x0))), where f and g are vector functions of the residual 
        /// vector r(x).
        /// </summary>
        /// <param name="lagrangeMultipliers">
        /// The lagrange multipliers which is the unknown vector in the linear system solved by PCPG.
        /// </param>
        /// <param name="projectedPrecondResidual">
        /// The vector resulting from projecting and preconditioning the residual: z = inv(M) * P *  r, where r is the residual, 
        /// M is the preconditioner and P is the projection.
        ///</param>
        double EstimateResidualNormRatio(IVectorView lagrangeMultipliers, IVectorView projectedPrecondResidual);
    }
}
