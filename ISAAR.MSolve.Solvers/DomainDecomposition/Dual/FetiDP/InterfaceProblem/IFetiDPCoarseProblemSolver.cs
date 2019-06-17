using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
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

        //TODO: Perhaps corner nodes of each subdomain should be stored in FetiDPDofSeparator.
        void CreateAndInvertCoarseProblemMatrix(Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains, 
            FetiDPDofSeparator dofSeparator, Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers);

        Vector MultiplyInverseCoarseProblemMatrixTimes(Vector vector);

        void ReorderCornerDofs(FetiDPDofSeparator dofSeparator);
    }
}
