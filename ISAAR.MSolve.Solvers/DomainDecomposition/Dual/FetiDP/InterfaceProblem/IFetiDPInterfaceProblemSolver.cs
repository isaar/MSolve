using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public interface IFetiDPInterfaceProblemSolver
    {
        (Vector lagrangeMultipliers, Vector cornerDisplacements) Solve(FetiDPFlexibilityMatrix flexibility,
            IFetiPreconditioner preconditioner, Vector disconnectedDisplacements, double globalForcesNorm, 
            DualSolverLogger logger);
    }
}
