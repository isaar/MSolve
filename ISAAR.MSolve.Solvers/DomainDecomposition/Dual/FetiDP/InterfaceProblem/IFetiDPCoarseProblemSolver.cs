using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public interface IFetiDPCoarseProblemSolver
    {
        void ClearCoarseProblemMatrix();

        Vector CreateCoarseProblemRhs(FetiDPDofSeparator dofSeparator,
            Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers,
            Dictionary<int, Vector> fr, Dictionary<int, Vector> fbc);

        void CreateAndInvertCoarseProblemMatrix(FetiDPDofSeparator dofSeparator,
            Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers);

        Vector MultiplyInverseCoarseProblemMatrixTimes(Vector vector);
    }
}
