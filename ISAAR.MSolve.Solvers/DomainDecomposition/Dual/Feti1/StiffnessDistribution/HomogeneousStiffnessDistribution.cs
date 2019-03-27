using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    public class HomogeneousStiffnessDistribution : IFeti1StiffnessDistribution
    {
        //TODO: perhaps the should be removed from the interface methods
        private readonly Feti1DofSeparator dofSeparator;

        public HomogeneousStiffnessDistribution(IStructuralModel_v2 model, Feti1DofSeparator dofSeparator)
        {
            this.dofSeparator = dofSeparator;
            this.SubdomainGlobalConversion = new HomogeneousSubdomainGlobalConversion(model, dofSeparator);
        }

        public ISubdomainGlobalConversion SubdomainGlobalConversion { get; }

        public Dictionary<int, Matrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            LagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, Matrix> boundarySignedBooleanMatrices)
        {
            var matricesBpb = new Dictionary<int, Matrix>();
            foreach (int id in boundarySignedBooleanMatrices.Keys)
            {
                Matrix inverseMultiplicities = InvertDofMultiplicities(dofSeparator.BoundaryDofMultiplicities[id]);
                matricesBpb[id] = boundarySignedBooleanMatrices[id].MultiplyRight(inverseMultiplicities);
            }
            return matricesBpb;
        }

        //TODO: Perhaps I could use int[] -> double[] -> DiagonalMatrix -> .Invert()
        //TODO: This can be parallelized OpenMP style.
        private Matrix InvertDofMultiplicities(int[] multiplicities)
        {
            var inverse = Matrix.CreateZero(multiplicities.Length, multiplicities.Length);
            for (int i = 0; i < multiplicities.Length; ++i) inverse[i, i] = 1.0 / multiplicities[i];
            return inverse;
        }

        //TODO: Add this to linear algebra: ScaleRow, ScaleRows, ScaleColumn, ScaleColumns, plus ~IntoThis versions.
        //      Alternatively only add left and right multiplication with diagonal matrices.
        private Matrix ScaleMatrixColumnsIntoThis(Matrix matrix, double[] scalars)
        {
            var result = Matrix.CreateZero(matrix.NumRows, matrix.NumColumns);
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i) result[i, j] = matrix[i, j] * scalars[j];
            }
            return result;
        }

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
    }
}
