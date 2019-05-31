using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPFlexibilityMatrix
    {
        private readonly Dictionary<int, SignedBooleanMatrixRowMajor> Br;
        private readonly Dictionary<int, CholeskyFull> factorizedKrr;
        private readonly Dictionary<int, Matrix> Krc;
        private readonly Dictionary<int, Matrix> Lc;
        private readonly int numCornerDofs;

        public FetiDPFlexibilityMatrix(Dictionary<int, CholeskyFull> factorizedKrr, Dictionary<int, Matrix> Krc, 
            FetiDPLagrangeMultipliersEnumerator lagrangeEnumerator, FetiDPDofSeparator dofSeparator)
        {
            this.Br = lagrangeEnumerator.BooleanMatrices;
            this.Lc = dofSeparator.CornerBooleanMatrices;
            this.factorizedKrr = factorizedKrr;
            this.Krc = Krc;
            this.numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            this.Order = lagrangeEnumerator.NumLagrangeMultipliers;
        }

        //TODO: This only matches FIrr, FIrc has a different dimension. Remove this property.
        public int Order { get; }

        public void MultiplyFIrr(Vector lhs, Vector rhs)
        {
            Preconditions.CheckMultiplicationDimensions(Order, lhs.Length);
            Preconditions.CheckSystemSolutionDimensions(Order, rhs.Length);

            // FIrr[s] * x = sum_over_s( Br[s] * (inv(Krr[s]) * (Br[s]^T * x)) )
            rhs.Clear();
            foreach (int s in Br.Keys)
            {
                Vector temp = Br[s].Multiply(lhs, true);
                temp = factorizedKrr[s].SolveLinearSystem(temp);
                temp = Br[s].Multiply(temp);
                rhs.AddIntoThis(temp);
            }
        }

        public Vector MultiplyFIrc(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(numCornerDofs, vector.Length);

            // FIrc[s] * x = sum_over_s( Br[s] * (inv(Krr[s]) * (Krc[s] * (Lc[s] * x))) )
            var result = Vector.CreateZero(Order);
            foreach (int s in Br.Keys)
            {
                Vector temp = Lc[s].Multiply(vector);
                temp = Krc[s].Multiply(temp);
                temp = factorizedKrr[s].SolveLinearSystem(temp);
                temp = Br[s].Multiply(temp);
                result.AddIntoThis(temp);
            }
            return result;
        }

        public Vector MultiplyTransposedFIrc(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(Order, vector.Length);
                       
            // FIrc[s]^T * x = sum_over_s( Lc[s]^T * (Krc[s]^T * (inv(Krr[s]) * (Br[s]^T * x))) )
            var result = Vector.CreateZero(numCornerDofs);
            foreach (int s in Br.Keys)
            {
                Vector temp = Br[s].Multiply(vector, true);
                temp = factorizedKrr[s].SolveLinearSystem(temp);
                temp = Krc[s].Multiply(temp, true);
                temp = Lc[s].Multiply(temp, true);
                result.AddIntoThis(temp);
            }
            return result;
        }
    }
}
