using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPFlexibilityMatrix
    {
        private readonly Dictionary<int, SignedBooleanMatrixColMajor> Br;
        private readonly Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers;
        private readonly Dictionary<int, UnsignedBooleanMatrix> Lc;
        private readonly int numCornerDofs;

        public FetiDPFlexibilityMatrix(FetiDPDofSeparator dofSeparator, FetiDPLagrangeMultipliersEnumerator lagrangeEnumerator, 
            Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers)
        {
            this.Br = lagrangeEnumerator.BooleanMatrices;
            this.Lc = dofSeparator.CornerBooleanMatrices;
            this.matrixManagers = matrixManagers;
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
                IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];
                Vector temp = Br[s].Multiply(lhs, true);
                temp = matrices.MultiplyInverseKrrTimes(temp);
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
                IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];
                Vector temp = Lc[s].Multiply(vector);
                temp = matrices.MultiplyKrcTimes(temp);
                temp = matrices.MultiplyInverseKrrTimes(temp);
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
                IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];
                Vector temp = Br[s].Multiply(vector, true);
                temp = matrices.MultiplyInverseKrrTimes(temp);
                temp = matrices.MultiplyKcrTimes(temp);
                temp = Lc[s].Multiply(temp, true);
                result.AddIntoThis(temp);
            }
            return result;
        }
    }
}
