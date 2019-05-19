using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

//TODO: This could be split into an interface with the same name and an IFtiDPCoarseProblemSolver.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public interface IFetiDPInterfaceProblemSolver
    {
        void ClearCoarseProblemMatrix();

        void CreateCoarseProblemMatrix(FetiDPDofSeparator dofSeparator, Dictionary<int, CholeskyFull> factorizedKrr,
            Dictionary<int, Matrix> Krc, Dictionary<int, Matrix> Kcc);

        Vector CreateCoarseProblemRhs(FetiDPDofSeparator dofSeparator, Dictionary<int, CholeskyFull> factorizedKrr,
            Dictionary<int, Matrix> Krc, Dictionary<int, Vector> fr, Dictionary<int, Vector> fbc);

        Vector SolveCoarseProblem(Vector rhs);

        (Vector lagrangeMultipliers, Vector cornerDisplacements) SolveInterfaceProblem(FetiDPFlexibilityMatrix flexibility, 
            IFetiPreconditioner preconditioner, Vector globalFcStar, Vector dr,
            double globalForcesNorm, DualSolverLogger logger);
    }
}
