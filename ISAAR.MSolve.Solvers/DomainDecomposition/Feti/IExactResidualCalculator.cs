using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    /// <summary>
    /// Called by PCPG.
    /// </summary>
    internal delegate double CalculateExactResidualNorm(Vector lagrangeMultipliers);

    internal interface IExactResidualCalculator
    {
        /// <summary>
        /// Called by FetiSolver after calculating the global displacements.
        /// </summary>
        /// <param name="globalDisplacements"></param>
        double CalculateExactResidualNorm(IVectorView globalDisplacements);
    }
}
