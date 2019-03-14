using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public class HeterogeneousStiffnessDistribution //: IStiffnessDistribution
    {
        public INodalLoadDistributor NodalLoadDistributor { get; } = new HeterogeneousNodalLoadDistributor();

        public Matrix CalcBoundaryPreconditioningSignedBooleanMatrix(int subdomainID, Matrix boundarySignedBooleanMatrix,
            int[] boundaryDofsMultiplicity, int[] boundaryDofs, Dictionary<int, IMatrixView> stiffnesses)
        {
            throw new NotImplementedException();
        }

        //TODO: this is also done when distributing the nodal loads. Do it here only and use the inv(Db) matrix there.
        //TODO: Kbb is also calculated for most preconditioners. Just take its diagonal and invert.
        private Matrix InvertBoundaryDofStiffnesses(int[] boundaryDofs, IMatrixView stiffness)
        {
            var inverse = Matrix.CreateZero(boundaryDofs.Length, boundaryDofs.Length);
            for (int i = 0; i < boundaryDofs.Length; ++i) inverse[i, i] = 1.0 / stiffness[i, i];
            return inverse;
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
