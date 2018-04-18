using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Vectors
{
    public class SparseVector: IVectorView
    {
        private readonly double[] values;

        /// <summary>
        /// They must be sorted
        /// </summary>
        private readonly int[] rowIndices;

        private SparseVector(int length, double[] values, int[] rowIndices)
        {
            this.Length = length;
            this.values = values;
            this.rowIndices = rowIndices;
        }

        public int[] InternalRowIndices { get { return rowIndices; } }

        public double[] InternalValues { get { return values; } }

        public int Length { get; }

        public double this[int idx]
        {
            get
            {
                int sparseIdx = FindSparseIndexOf(idx);
                if (sparseIdx < 0) return 0.0;
                else return values[sparseIdx];
            }
        }

        public static SparseVector CreateFromArrays(int length, double[] values, int[] rowIndices, 
            bool checkInput, bool sortInput)
        {
            bool verifiedSorted = false;
            if (sortInput)
            {
                Array.Sort<int, double>(rowIndices, values);
                verifiedSorted = true;
            }

            if (checkInput)
            {
                int nnz = rowIndices.Length;
                if (values.Length != nnz)
                {
                    throw new ArgumentException("The length of the values and row indices arrays must be equal (and equal"
                        + $" to the number of non zero entries), but were {values.Length} and {rowIndices.Length} respectively");
                }
                if (!verifiedSorted) // if they have already been sorted we do not need to check them again.
                {
                    for (int i = 1; i < nnz; ++i)
                    {
                        if (rowIndices[i] <= rowIndices[i - 1]) throw new ArgumentException("The row indices must be sorted");
                    }
                }
                int maxIndex = rowIndices[nnz - 1]; // first sort them
                if (maxIndex >= length) 
                {
                    throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
                }
            }
            return new SparseVector(length, values, rowIndices);
        }

        public static SparseVector CreateFromDictionary(int length, Dictionary<int, double> rowsValues)
        {
            double[] values = new double[rowsValues.Count];
            int[] rowIndices = new int[rowsValues.Count];
            int nnz = 0;
            foreach (var rowValPair in rowsValues.OrderBy(pair => pair.Key))
            {
                rowIndices[nnz] = rowValPair.Key;
                values[nnz] = rowValPair.Value;
                ++nnz;
            }

            int maxIndex = rowIndices[nnz - 1]; // first sort them
            if (maxIndex >= length)
            {
                throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
            }
            return new SparseVector(length, values, rowIndices);
        }

        public static SparseVector CreateFromDictionary(int length, SortedDictionary<int, double> rowsValues)
        {
            double[] values = new double[rowsValues.Count];
            int[] rowIndices = new int[rowsValues.Count];
            int nnz = 0;
            foreach (var rowValPair in rowsValues)
            {
                rowIndices[nnz] = rowValPair.Key;
                values[nnz] = rowValPair.Value;
                ++nnz;
            }
            
            int maxIndex = rowIndices[nnz - 1]; // first sort them
            if (maxIndex >= length)
            {
                throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
            }
            return new SparseVector(length, values, rowIndices);
        }

        public double[] CopyToArray()
        {
            double[] result = new double[Length];
            for (int i = 0; i < values.Length; ++i) result[rowIndices[i]] = values[i];
            return result;
        }

        public VectorMKL CopyToFullVector()
        {
            return VectorMKL.CreateFromArray(CopyToArray(), false);
        }

        public int CountNonZeros()
        {
            return values.Length;
        }

        public IVectorView DoEntrywise(IVectorView other, Func<double, double, double> binaryOperation)
        {
            if (other is SparseVector otherSparse) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherSparse))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherSparse.values[i]);
                    }
                    return new SparseVector(Length, resultValues, rowIndices);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        public IVectorView DoToAllEntries(Func<double, double> unaryOperation)
        {
            // Only apply the operation on non zero entries
            double[] newValues = new double[values.Length];
            for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
                int[] rowIndicesCopy = new int[rowIndices.Length];
                Array.Copy(rowIndices, rowIndicesCopy, rowIndices.Length);
                return new SparseVector(Length, newValues, rowIndicesCopy);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new SparseVector(Length, newValues, rowIndices).CopyToFullVector();
            }
        }

        public double DotProduct(IVectorView other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            double sum = 0;
            for (int i = 0; i < values.Length; ++i) sum += values[i] * other[rowIndices[i]];
            return sum;
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = values.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(values[i], aggregator);
            aggregator = processZeros(Length - nnz, aggregator);
            return finalize(aggregator);
        }

        private int FindSparseIndexOf(int denseIdx)
        {
            Preconditions.CheckIndex1D(this, denseIdx);
            return Array.BinarySearch<int>(rowIndices, denseIdx); // only works if rowIndices are sorted!!!
        }

        private bool HasSameIndexer(SparseVector other)
        {
            return this.rowIndices == other.rowIndices;
        }
    }
}
