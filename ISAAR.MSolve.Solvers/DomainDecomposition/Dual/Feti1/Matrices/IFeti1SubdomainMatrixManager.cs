using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Perhaps I should split this into subinterfaces: 
//      1) The matrices Kbi, Kii, diag(Kii) are not always needed.
//      2) The factorization iof Kff depends on if the user chooses analytical vs algebraic computation of the rbms.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices
{
    /// <summary>
    /// Implements the linear algebra operations needed by <see cref="Feti1Solver"/> depending on the underlying matrix storage
    /// format. All the matrices represented by this interface belong to a single subdomain.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IFeti1SubdomainMatrixManager : IFetiSubdomainMatrixManager
    {
        // The FETI-1 solver should access them from here, so that they can be created/disposed together with inv(Kff).
        List<Vector> RigidBodyModes { get; } 

        void InvertKff(double factorizationTolerance, bool inPlace);

        Vector MultiplyInverseKffTimes(Vector vector);
    }
}
