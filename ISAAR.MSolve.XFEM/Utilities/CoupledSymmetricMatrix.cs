using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Utilities
{
    internal class CoupledSymmetricMatrix : IMatrix
    {
        private readonly int dofs2Start;
        private readonly IMatrix matrix11;
        private readonly IMatrix matrix12;
        private readonly IMatrix matrix22;

        public CoupledSymmetricMatrix(IMatrix matrix11, IMatrix matrix12, IMatrix matrix22)
        {
#if DEBUG
            Preconditions.CheckSquare(matrix11);
            Preconditions.CheckSquare(matrix22);
            if (matrix12.NumRows != matrix11.NumRows) throw new NonMatchingDimensionsException(
                "Matrices 11 and 12 must have the same number of rows");
            if (matrix22.NumColumns != matrix22.NumColumns) throw new NonMatchingDimensionsException(
                "Matrices 22 and 12 must have the same number of columns");
#endif
            this.matrix11 = matrix11;
            this.matrix12 = matrix12;
            this.matrix22 = matrix22;
            this.dofs2Start = matrix11.NumRows;
            this.NumRows = matrix11.NumRows + matrix22.NumRows;
            this.NumColumns = matrix11.NumColumns + matrix22.NumColumns;
        }

        public double this[int rowIdx, int colIdx]
        {
            get
            {
                if (colIdx < dofs2Start)
                {
                    if (rowIdx < dofs2Start) return matrix11[rowIdx, colIdx];
                    else return matrix12[colIdx, rowIdx - dofs2Start];
                }
                else
                {
                    if (rowIdx < dofs2Start) return matrix12[rowIdx, colIdx - dofs2Start];
                    else return matrix22[rowIdx - dofs2Start, colIdx - dofs2Start];
                }
            }
        }

        public int NumColumns { get; }

        public int NumRows { get; }

        public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
            => LinearCombination(1.0, otherMatrix, otherCoefficient);

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
            => LinearCombinationIntoThis(1.0, otherMatrix, otherCoefficient);

        public void Clear()
        {
            matrix11.Clear();
            matrix12.Clear();
            matrix22.Clear();
        }

        public IMatrix Copy(bool copyIndexingData = false)
            => new CoupledSymmetricMatrix(matrix11.Copy(), matrix12.Copy(), matrix22.Copy());

        public Matrix CopyToFullMatrix() => DenseStrategies.CopyToFullMatrix(this);

        public IMatrix DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
            => DenseStrategies.DoEntrywise(this, matrix, binaryOperation);

        public void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation)
            => throw new NotImplementedException();

        public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new CoupledSymmetricMatrix(matrix11.DoToAllEntries(unaryOperation), matrix12.DoToAllEntries(unaryOperation),
                matrix22.DoToAllEntries(unaryOperation));
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            matrix11.DoToAllEntriesIntoThis(unaryOperation);
            matrix12.DoToAllEntriesIntoThis(unaryOperation);
            matrix22.DoToAllEntriesIntoThis(unaryOperation);
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13) => DenseStrategies.AreEqual(this, other, tolerance);

        public Vector GetColumn(int colIndex) => DenseStrategies.GetColumn(this, colIndex);

        public Vector GetRow(int rowIndex) => DenseStrategies.GetRow(this, rowIndex);

        public IMatrix GetSubmatrix(int[] rowIndices, int[] colIndices) 
            => DenseStrategies.GetSubmatrix(this, rowIndices, colIndices);

        public IMatrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
            => DenseStrategies.GetSubmatrix(this, rowStartInclusive, rowEndExclusive, colStartInclusive, colEndExclusive);

        public IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
            => DenseStrategies.LinearCombination(this, thisCoefficient, otherMatrix, otherCoefficient);

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
            => throw new NotImplementedException();

        public IVector Multiply(IVectorView vector, bool transposeThis = false)
            => throw new NotImplementedException(); //TODO: Implement multiply subvector or even better GetSubvectorView, which just offsets indices.

        public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false) //TODO: implement this properly. It is used for Kfc * uc
        {
            //TODO: implement efficient multiplication with subvectors/offsets.
            if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
            {
                Preconditions.CheckMultiplicationDimensions(this.NumColumns, lhsDense.Length);
                Preconditions.CheckSystemSolutionDimensions(this.NumRows, rhsDense.Length);

                int n1 = matrix11.NumColumns;
                int n2 = matrix22.NumColumns;
                rhsDense.Clear();
                rhsDense.AddSubvectorIntoThis(0, matrix11.Multiply(lhsDense.GetSubvector(0, n1)), 0, n1);
                rhsDense.AddSubvectorIntoThis(0, matrix12.Multiply(lhsDense.GetSubvector(n1, n1 + n2)), 0, n1);
                rhsDense.AddSubvectorIntoThis(n1, matrix12.Multiply(lhsDense.GetSubvector(0, n1), true), 0, n2);
                rhsDense.AddSubvectorIntoThis(n1, matrix22.Multiply(lhsDense.GetSubvector(n1, n1 + n2)), 0, n2);
            }
            else DenseStrategies.MultiplyIntoResult(this, lhsVector, rhsVector, transposeThis);
        }

        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
            => throw new NotImplementedException();

        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
            => throw new NotImplementedException();

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
            => throw new NotImplementedException();

        public IMatrix Scale(double scalar)
            => new CoupledSymmetricMatrix(matrix11.Scale(scalar), matrix12.Scale(scalar), matrix22.Scale(scalar));

        public void ScaleIntoThis(double scalar)
        {
            matrix11.ScaleIntoThis(scalar);
            matrix12.ScaleIntoThis(scalar);
            matrix22.ScaleIntoThis(scalar);
        }

        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            if (colIdx < dofs2Start)
            {
                if (rowIdx < dofs2Start) matrix11.SetEntryRespectingPattern(rowIdx, colIdx, value);
                else matrix12.SetEntryRespectingPattern(colIdx, rowIdx - dofs2Start, value);
            }
            else
            {
                if (rowIdx < dofs2Start) matrix12.SetEntryRespectingPattern(rowIdx, colIdx - dofs2Start, value);
                else matrix22.SetEntryRespectingPattern(rowIdx - dofs2Start, colIdx - dofs2Start, value);
            }
        }

        public IMatrix Transpose() => Copy();
    }
}
