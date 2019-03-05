using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    internal class InterfaceFlexibilityMatrix 
    {
        private readonly ContinuityEquationsCalculator continuityEquations;
        private readonly Dictionary<int, SemidefiniteCholeskySkyline> factorizations;

        internal InterfaceFlexibilityMatrix(Dictionary<int, SemidefiniteCholeskySkyline> factorizations,
            ContinuityEquationsCalculator continuityEquations)
        {
            this.continuityEquations = continuityEquations;
            this.factorizations = factorizations;
        }

        internal int Order { get; }

        internal void Multiply(Vector lhs, Vector rhs)
        {
            rhs.Clear(); //TODO: perhaps this should be done outside.
            foreach (var keyFactor in factorizations)
            {
                int id = keyFactor.Key;
                SemidefiniteCholeskySkyline factor = keyFactor.Value;
                SignedBooleanMatrix boolean = continuityEquations.BooleanMatrices[id];
                Vector FBx = factor.MultiplyGeneralizedInverseMatrixTimesVector(boolean.Multiply(lhs, true)); 
                Vector BFBx = boolean.Multiply(FBx, false);
                rhs.AddIntoThis(BFBx);
            }
        }
    }
}
