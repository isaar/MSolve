using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition
{
    public interface IDofSeparator
    {
        /// <summary>
        /// Dofs where Lagrange multipliers will be applied. Depending on the domain decomposition method, these could be all
        /// boundary dofs (e.g. FETI-1) or a subset of them (e.g. boundary dofs minus corner dofs in FETI-DP).
        /// </summary>
        Dictionary<INode, IDofType[]> DualDofs { get; }
    }
}
