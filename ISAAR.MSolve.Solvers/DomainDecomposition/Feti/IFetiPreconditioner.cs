using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public interface IFetiPreconditioner
    {
        void SolveLinearSystem(Vector rhs, Vector lhs);
    }

    public interface IFetiPreconditionerFactory
    {
        IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution, DofSeparator dofSeparator, 
            ContinuityEquationsCalculator continuityEquations, Dictionary<int, IMatrixView> stiffnessMatrices);
    }
}
