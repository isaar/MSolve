using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Pcpg
{
    public interface IFeti1PcgResidualCalculator
    {
        /// <summary>
        /// Calculates norm2(f(r(λ))), where f is a vector function of the residual vector r(x) and λ is the vector of lagrange 
        /// multipliers.
        /// </summary>
        /// <param name="lagrangeMultipliers">
        /// The lagrange multipliers which is the unknown vector in the linear system solved by PCPG.
        /// </param>
        /// <param name="projectedPrecondResidual">
        /// The vector resulting from projecting and preconditioning the residual: z = inv(M) * P *  r, where r is the residual, 
        /// M is the preconditioner and P is the projection.
        ///</param>
        double CalcResidualNormRatio(IVectorView lagrangeMultipliers, IVectorView projectedPrecondResidual);

        /// <summary>
        /// Initializes internal data.
        /// </summary>
        void Initialize();
    }
}
