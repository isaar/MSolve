using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1
{
    public class Feti1FlexibilityMatrix : IInterfaceFlexibilityMatrix
    {
        private readonly LagrangeMultipliersEnumerator lagrangeEnumerator;
        private readonly Dictionary<int, SemidefiniteCholeskySkyline> factorizations;

        internal Feti1FlexibilityMatrix(Dictionary<int, SemidefiniteCholeskySkyline> factorizations,
            LagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            this.lagrangeEnumerator = lagrangeEnumerator;
            this.factorizations = factorizations;
            this.Order = lagrangeEnumerator.NumLagrangeMultipliers;
        }

        public int Order { get; }

        public void Multiply(Vector lhs, Vector rhs)
        {
            rhs.Clear(); //TODO: perhaps this should be done outside.
            foreach (var keyFactor in factorizations)
            {
                int id = keyFactor.Key;
                SemidefiniteCholeskySkyline factor = keyFactor.Value;
                SignedBooleanMatrix boolean = lagrangeEnumerator.BooleanMatrices[id];
                Vector FBx = factor.MultiplyGeneralizedInverseMatrixTimesVector(boolean.Multiply(lhs, true)); 
                Vector BFBx = boolean.Multiply(FBx, false);
                rhs.AddIntoThis(BFBx);
            }
        }

        public Vector Multiply(Vector lhs)
        {
            var rhs = Vector.CreateZero(Order);
            foreach (var keyFactor in factorizations)
            {
                int id = keyFactor.Key;
                SemidefiniteCholeskySkyline factor = keyFactor.Value;
                SignedBooleanMatrix boolean = lagrangeEnumerator.BooleanMatrices[id];
                Vector FBx = factor.MultiplyGeneralizedInverseMatrixTimesVector(boolean.Multiply(lhs, true));
                Vector BFBx = boolean.Multiply(FBx, false);
                rhs.AddIntoThis(BFBx);
            }
            return rhs;
        }
    }
}
