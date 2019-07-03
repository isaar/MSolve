using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: During initialization, the solver and its verious strategies should inform IFetiDPSubdomainMatrixManager what matrices
//      will be necessary. IFetiDPSubdomainMatrixManager should then determine the correct order they must be created in and
//      notify the solver and each strategy when they are ready for consumption. Also once a matrix has been fully used, 
//      it should be cleared to conserve memory. This also applies for Kff.
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

        //TODO: Perhaps this should be done by the dofSeparator since its state will be modified. What about the cases where no 
        //      reordering is applied?
        //TODO: Bad design. All this time the matrix manager had access to only 1 subdomain and now I pass it an object that
        //      stores dof data for all subdomains and which subdomain to use.
        void ReorderInternalDofs(Feti1DofSeparator dofSeparator, ISubdomain subdomain);
    }
}
