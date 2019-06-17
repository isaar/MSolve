using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: During initialization, the solver and its verious strategies should inform IFetiDPSubdomainMatrixManager what matrices
//      will be necessary. IFetiDPSubdomainMatrixManager should then determine the correct order they must be created in and
//      notify the solver and each strategy when they are ready for consumption. Also once a matrix has been fully used, 
//      it should be cleared to conserve memory. This also applies for Kff.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices
{
    /// <summary>
    /// Implements the linear algebra operations needed by <see cref="Feti1Solver"/> depending on the underlying matrix storage
    /// format. All the matrices represented by this interface belong to a single subdomain.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IFetiDPSubdomainMatrixManager : IFetiSubdomainMatrixManager
    {
        Matrix SchurComplementOfRemainderDofs { get; } //TODO: Perhaps static condensations should be handled by a different interface

        void CalcSchurComplementOfRemainderDofs();

        //TODO: E.g. Once Kcc* is calculated Kcc and Krc can be cleared. There are 2 options:
        //      a) Each matrix must be able to be cleared independently, if the FETI-DP solver and its strategies decide when.
        //      b) Otherwise this matrix decides when to clear what and these methods are optional/risky.
        //void ClearKcc();
        //void ClearKcrKrc();

        void ExtractKcc(int[] cornerDofs);  //TODO: perhaps SymmetricMatrix
        void ExtractKcrKrc(int[] cornerDofs, int[] remainderDofs); //TODO: perhaps CSR or CSC
        void ExtractKrr(int[] remainderDofs); //TODO: perhaps SkylineMatrix or SymmetricCSC

        void InvertKrr(bool inPlace);

        Vector MultiplyInverseKrrTimes(Vector vector);
        Vector MultiplyKccTimes(Vector vector);
        Vector MultiplyKcrTimes(Vector vector);
        Vector MultiplyKrcTimes(Vector vector);

        //TODO: Perhaps these reorderings should be done by the dofSeparator since its state will be modified. What about the  
        //      cases where no reordering is applied (e.g. lumped preconditioner)?
        //TODO: Bad design. All this time the matrix manager had access to only 1 subdomain and now I pass it an object that
        //      stores dof data for all subdomains and which subdomain to use.
        void ReorderInternalDofs(FetiDPDofSeparator dofSeparator, ISubdomain subdomain);

        //TODO: This must be called before creating mapping matrices (Br, Bc) or even processing boundary/internal dofs
        void ReorderRemainderDofs(FetiDPDofSeparator dofSeparator, ISubdomain subdomain);
    }
}
