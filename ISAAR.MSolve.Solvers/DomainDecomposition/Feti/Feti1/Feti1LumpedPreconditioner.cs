using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1
{
    public class Feti1LumpedPreconditioner : IFetiPreconditioner
    {
        private readonly Dictionary<int, Matrix> preconditioningBoundarySignedBooleanMatrices;
        private readonly Dictionary<int, Matrix> stiffnessesBoundaryBoundary;
        private readonly int[] subdomainIDs;

        private Feti1LumpedPreconditioner(int[] subdomainIDs, Dictionary<int, Matrix> stiffnessesBoundaryBoundary,
            Dictionary<int, Matrix> preconditioningBoundarySignedBooleanMatrices)
        {
            this.subdomainIDs = subdomainIDs;
            this.preconditioningBoundarySignedBooleanMatrices = preconditioningBoundarySignedBooleanMatrices;
            this.stiffnessesBoundaryBoundary = stiffnessesBoundaryBoundary;
        }

        //TODO: This can be moved to a base class. Only the S matrix is different for these preconditioners. 
        //      Other ones might be different though.
        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            // inv(F) * y = Bpb * Kbb * Bpb^T * y
            lhs.Clear(); //TODO: this should be avoided
            foreach (int id in subdomainIDs)
            {
                Matrix Bpb = preconditioningBoundarySignedBooleanMatrices[id];
                Matrix Kbb = stiffnessesBoundaryBoundary[id];
                Vector subdomainContribution = Bpb.Multiply(Kbb.Multiply(Bpb.Multiply(rhs, true)));
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public class Factory : FetiPreconditionerFactoryBase
        {
            public override IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution,
                DofSeparator dofSeparator, LagrangeMultipliersEnumerator lagrangeEnumerator,
                Dictionary<int, IMatrixView> stiffnessMatrices)
            {
                int[] subdomainIDs = dofSeparator.BoundaryDofs.Keys.ToArray();
                Dictionary<int, Matrix> boundaryBooleans = CalcBoundaryPreconditioningBooleanMatrices(stiffnessDistribution, 
                    dofSeparator, lagrangeEnumerator);
                Dictionary<int, Matrix> stiffnessesBoundaryBoundary = 
                    ExtractStiffnessesBoundaryBoundary(dofSeparator, stiffnessMatrices);
                return new Feti1LumpedPreconditioner(subdomainIDs, stiffnessesBoundaryBoundary, boundaryBooleans);
            }
        }
    }
}
