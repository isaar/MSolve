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

        //TODO: This is a bad idea if we need to create multiple solvers. It will either pass default settings to all settings 
        //      but the first or it will use the same settings, which might be mutable objects. If a method needs to create 
        //      many solvers, it would be best to pass a lambda that instantiates one ISolverBuilder, which then will be used to 
        //      create one solver. 
        ISolverBuilder Clone(); 
    }
}
