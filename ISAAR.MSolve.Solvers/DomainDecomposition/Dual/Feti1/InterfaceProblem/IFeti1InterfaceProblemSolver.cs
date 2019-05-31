using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.InterfaceProblem
{
    public interface IFeti1InterfaceProblemSolver
    {
        Vector CalcLagrangeMultipliers(Feti1FlexibilityMatrix flexibility, IFetiPreconditioner preconditioner,
            Feti1Projection projection, Vector disconnectedDisplacements, Vector rigidBodyModesWork, double globalForcesNorm,
            SolverLogger logger);
    }
}
