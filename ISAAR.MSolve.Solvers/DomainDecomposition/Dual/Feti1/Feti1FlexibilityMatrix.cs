using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    public class Feti1FlexibilityMatrix : IInterfaceFlexibilityMatrix
    {
        private readonly Feti1LagrangeMultipliersEnumerator lagrangeEnumerator;
        private readonly Dictionary<int, IFeti1SubdomainMatrixManager> matrixManagers;

        internal Feti1FlexibilityMatrix(Dictionary<int, IFeti1SubdomainMatrixManager> matrixManagers,
            Feti1LagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            this.lagrangeEnumerator = lagrangeEnumerator;
            this.matrixManagers = matrixManagers;
            this.Order = lagrangeEnumerator.NumLagrangeMultipliers;
        }

        public int Order { get; }

        public void Multiply(Vector lhs, Vector rhs)
        {
            rhs.Clear(); //TODO: perhaps this should be done outside.
            foreach (var keyFactor in matrixManagers)
            {
                int id = keyFactor.Key;
                IFeti1SubdomainMatrixManager matrixManager = keyFactor.Value;
                SignedBooleanMatrixColMajor B = lagrangeEnumerator.BooleanMatrices[id];
                Vector FBx = matrixManager.MultiplyInverseKffTimes(B.Multiply(lhs, true)); 
                Vector BFBx = B.Multiply(FBx, false);
                rhs.AddIntoThis(BFBx);
            }
        }

        public Vector Multiply(Vector lhs)
        {
            var rhs = Vector.CreateZero(Order);
            foreach (var keyFactor in matrixManagers)
            {
                int id = keyFactor.Key;
                IFeti1SubdomainMatrixManager matrixManager = keyFactor.Value;
                SignedBooleanMatrixColMajor B = lagrangeEnumerator.BooleanMatrices[id];
                Vector FBx = matrixManager.MultiplyInverseKffTimes(B.Multiply(lhs, true));
                Vector BFBx = B.Multiply(FBx, false);
                rhs.AddIntoThis(BFBx);
            }
            return rhs;
        }
    }
}
