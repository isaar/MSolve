using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public class HomogeneousStiffnessDistribution : IStiffnessDistribution
    {
        private readonly Dictionary<int, int[]> boundaryDofMultiplicities;
        private readonly IDofSeparator dofSeparator;

        public HomogeneousStiffnessDistribution(IStructuralModel model, IDofSeparator dofSeparator)
        {
            this.dofSeparator = dofSeparator;
            this.boundaryDofMultiplicities = FindBoundaryDofMultiplicities(dofSeparator);
        }

        public double[] CalcBoundaryDofCoefficients(ISubdomain subdomain)
        {
            int[] multiplicites = boundaryDofMultiplicities[subdomain.ID];
            var inverseMultiplicites = new double[multiplicites.Length];
            for (int i = 0; i < multiplicites.Length; ++i) inverseMultiplicites[i] = 1.0 / multiplicites[i];
            return inverseMultiplicites;
        }

        public Dictionary<int, double> CalcBoundaryDofCoefficients(INode node, IDofType dofType)
        {
            var coeffs = new Dictionary<int, double>();
            double inverseMultiplicity = 1.0 / node.SubdomainsDictionary.Count;
            foreach (int subdomainID in node.SubdomainsDictionary.Keys) coeffs[subdomainID] = inverseMultiplicity;
            return coeffs;
        }

        public Dictionary<int, IMappingMatrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            ILagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, Matrix> boundarySignedBooleanMatrices)
        {
            return ScalingBooleanMatrixExplicit.CreateBpbOfSubdomains(this, lagrangeEnumerator, boundarySignedBooleanMatrices);
        }

        private static Dictionary<int, int[]> FindBoundaryDofMultiplicities(IDofSeparator dofSeparator)
        {
            var allMultiplicities = new Dictionary<int, int[]>();
            foreach (var idBoundaryDofs in dofSeparator.BoundaryDofs)
            {
                int subdomainID = idBoundaryDofs.Key;
                (INode node, IDofType dofType)[] boundaryDofs = idBoundaryDofs.Value;
                var multiplicities = new int[boundaryDofs.Length];
                for (int i = 0; i < boundaryDofs.Length; ++i)
                {
                    multiplicities[i] = boundaryDofs[i].node.SubdomainsDictionary.Count;
                }
                allMultiplicities[subdomainID] = multiplicities;
            }
            return allMultiplicities;
        }

        //TODO: Perhaps I could use int[] -> double[] -> DiagonalMatrix -> .Invert()
        //TODO: This can be parallelized OpenMP style.
        private static DiagonalMatrix InvertDofMultiplicities(int[] multiplicities)
        {
            var invMb = new double[multiplicities.Length];
            for (int i = 0; i < multiplicities.Length; ++i) invMb[i] = 1.0 / multiplicities[i];
            return DiagonalMatrix.CreateFromArray(invMb, false);
        }

        //TODO: This should be modified to CSR or CSC format and then benchmarked against the implicit alternative.
        /// <summary>
        /// Calculates the product Bpb = Bb * inv(Mb) explicitly, stores it and uses it for multiplications.
        /// </summary>
        private class ScalingBooleanMatrixExplicit : IMappingMatrix
        {
            private readonly Matrix explicitBpb;

            private ScalingBooleanMatrixExplicit(Matrix explicitBpb)
            {
                this.explicitBpb = explicitBpb;
            }

            public int NumColumns => explicitBpb.NumColumns;

            public int NumRows => explicitBpb.NumRows;

            public double this[int rowIdx, int colIdx] => explicitBpb[rowIdx, colIdx];

            internal static Dictionary<int, IMappingMatrix> CreateBpbOfSubdomains(
                HomogeneousStiffnessDistribution stiffnessDistribution, ILagrangeMultipliersEnumerator lagrangeEnumerator,
                Dictionary<int, Matrix> boundarySignedBooleanMatrices)
            {
                var matricesBpb = new Dictionary<int, IMappingMatrix>();
                foreach (int id in boundarySignedBooleanMatrices.Keys)
                {
                    Matrix Bb = boundarySignedBooleanMatrices[id];
                    DiagonalMatrix invMb = InvertDofMultiplicities(stiffnessDistribution.boundaryDofMultiplicities[id]);
                    Matrix Bpb = Bb.MultiplyRight(invMb.CopyToFullMatrix());
                    matricesBpb[id] = new ScalingBooleanMatrixExplicit(Bpb);
                }
                return matricesBpb;
            }

            public bool Equals(IIndexable2D other, double tolerance = 1E-13)
                => explicitBpb.Equals(other, tolerance);

            public Vector Multiply(Vector vector, bool transposeThis = false)
                => explicitBpb.Multiply(vector, transposeThis);

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
                => explicitBpb.MultiplyRight(other, transposeThis);
        }

        /// <summary>
        /// Stores the matrices Bb and inv(Mb). Matrix-vector and matrix-matrix multiplications with Bpb = Bb * inv(Mb) are
        /// performed implicitly, e.g. Bpb * x = Bb * (inv(Mb) * x).
        /// </summary>
        private class ScalingBooleanMatrixImplicit : IMappingMatrix
        {
            /// <summary>
            /// Signed boolean matrix with only the boundary dofs of the subdomain as columns. 
            /// </summary>
            private readonly Matrix Bb;

            /// <summary>
            /// Inverse of the diagonal matrix that stores the multiplicity of each boundary dof of the subdomain.
            /// </summary>
            private readonly DiagonalMatrix invMb;

            private ScalingBooleanMatrixImplicit(Matrix Bb, DiagonalMatrix invMb)
            {
                this.Bb = Bb;
                this.invMb = invMb;
            }

            public int NumColumns => invMb.NumColumns;

            public int NumRows => Bb.NumRows;

            internal static Dictionary<int, IMappingMatrix> CreateBpbOfSubdomains(
                HomogeneousStiffnessDistribution stiffnessDistribution, ILagrangeMultipliersEnumerator lagrangeEnumerator,
                Dictionary<int, Matrix> boundarySignedBooleanMatrices)
            {
                var matricesBpb = new Dictionary<int, IMappingMatrix>();
                foreach (int id in boundarySignedBooleanMatrices.Keys)
                {
                    Matrix Bb = boundarySignedBooleanMatrices[id];
                    DiagonalMatrix invMb = InvertDofMultiplicities(stiffnessDistribution.boundaryDofMultiplicities[id]);
                    matricesBpb[id] = new ScalingBooleanMatrixImplicit(Bb, invMb);
                }
                return matricesBpb;
            }

            public Vector Multiply(Vector vector, bool transposeThis = false)
            {
                if (transposeThis)
                {
                    // Bpb^T * x = (Bb * inv(Mb))^T * x = inv(Mb)^T * Bb^T * x = inv(Mb) * (Bb^T * x);
                    return invMb * Bb.Multiply(vector, true);
                }
                else
                {
                    // Bpb * x = Bb * (inv(Mb) * x)
                    return Bb.Multiply(invMb * vector);
                }
            }

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
            {
                if (transposeThis)
                {
                    // Bpb^T * X = (Bb * inv(Mb))^T * X = inv(Mb)^T * Bb^T * X = inv(Mb) * (Bb^T * X);
                    return invMb * Bb.MultiplyRight(other, true);
                }
                else
                {
                    // Bpb * X = Bb * (inv(Mb) * X)
                    return Bb.MultiplyRight(invMb * other);
                }
            }
                
        }

        #region older code
        ////TODO: Add this to linear algebra: ScaleRow, ScaleRows, ScaleColumn, ScaleColumns, plus ~IntoThis versions.
        ////      Alternatively only add left and right multiplication with diagonal matrices.
        //private static Matrix ScaleMatrixColumnsIntoThis(Matrix matrix, double[] scalars)
        //{
        //    var result = Matrix.CreateZero(matrix.NumRows, matrix.NumColumns);
        //    for (int j = 0; j < matrix.NumColumns; ++j)
        //    {
        //        for (int i = 0; i < matrix.NumRows; ++i) result[i, j] = matrix[i, j] * scalars[j];
        //    }
        //    return result;
        //}

        //public Matrix CalcBoundaryPreconditioningSignedBooleanMatrix(Matrix boundarySignedBooleanMatrix,
        //    int[] boundaryDofsMultiplicity)
        //{
        //    double[] inverseBoundaryDofsMultiplicity = InvertDofMultiplicities(boundaryDofsMultiplicity);
        //    return ScaleMatrixColumnsIntoThis(boundarySignedBooleanMatrix, inverseBoundaryDofsMultiplicity);
        //}

        ////TODO: Perhaps I could use int[] -> double[] -> DiagonalMatrix -> .Invert()
        ////TODO: This can be parallelized OpenMP style.
        //private DiagonalMatrix InvertDofMultiplicities(int[] multiplicities)
        //{
        //    var inverse = new double[multiplicities.Length];
        //    for (int i = 0; i < multiplicities.Length; ++i) inverse[i] = 1.0 / multiplicities[i];
        //    return DiagonalMatrix.CreateFromArray(inverse);
        //}

        ////TODO: Add this to linear algebra: ScaleRow, ScaleRows, ScaleColumn, ScaleColumns, plus ~IntoThis versions.
        ////      Alternatively only add left and right multiplication with diagonal matrices.
        //private Matrix ScaleMatrixColumnsIntoThis(Matrix matrix, double[] scalars)
        //{
        //    var result = Matrix.CreateZero(matrix.NumRows, matrix.NumColumns);
        //    for (int j = 0; j < matrix.NumColumns; ++j)
        //    {
        //        for (int i = 0; i < matrix.NumRows; ++i) result[i, j] = matrix[i, j] * scalars[j];
        //    }
        //    return result;
        //}
        #endregion
    }
}
