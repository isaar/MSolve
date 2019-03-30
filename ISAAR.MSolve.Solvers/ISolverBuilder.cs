using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers
{
    /// <summary>
    /// Builder for classes implementing <see cref="ISolver_v2"/>.
    /// Authors: Gerasimos Sotiropoulos, Serafeim Bakalakos
    /// </summary>
    public interface ISolverBuilder
    {
        /// <summary>
        /// Creates a new solver for the provided model.
        /// </summary>
        /// <param name="model">The model that will be analyzed.</param>
        ISolver_v2 BuildSolver(IStructuralModel_v2 model);
    }
}
