using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public class SparseVector: IVectorView
    {
        private readonly double[] values;

        /// <summary>
        /// They must be sorted
        /// </summary>
        private readonly int[] indices;

        private SparseVector(int length, double[] values, int[] indices)
        {
            this.Length = length;
            this.values = values;
            this.indices = indices;
        }

        internal int[] InternalIndices { get { return indices; } }

        internal double[] InternalValues { get { return values; } }

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

        public static SparseVector CreateFromArrays(int length, double[] values, int[] indices, 
            bool checkInput, bool sortInput)
        {
            bool verifiedSorted = false;
            if (sortInput)
            {
                Array.Sort<int, double>(indices, values);
                verifiedSorted = true;
            }

            if (checkInput)
            {
                int nnz = indices.Length;
                if (values.Length != nnz)
                {
                    throw new ArgumentException("The length of the values and indices arrays must be equal (and equal"
                        + $" to the number of non zero entries), but were {values.Length} and {indices.Length} respectively");
                }
                if (!verifiedSorted) // if they have already been sorted we do not need to check them again.
                {
                    for (int i = 1; i < nnz; ++i)
                    {
                        if (indices[i] <= indices[i - 1]) throw new ArgumentException("The indices array must be sorted");
                    }
                }
                int maxIndex = indices[nnz - 1]; // first sort them
                if (maxIndex >= length) 
                {
                    throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
                }
            }
            return new SparseVector(length, values, indices);
        }

        // TODO: add a version with tolerance
        public static SparseVector CreateFromDense(double[] denseArray)
        {
            var indicesValues = new SortedDictionary<int, double>();
            for (int i = 0; i < denseArray.Length; ++i)
            {
                if (denseArray[i] != 0) indicesValues.Add(i, denseArray[i]);
            }
            return CreateFromDictionary(denseArray.Length, indicesValues);
        }

        public static SparseVector CreateFromDense(double[] denseArray, double tolerance)
        {
            var indicesValues = new SortedDictionary<int, double>();
            for (int i = 0; i < denseArray.Length; ++i)
            {
                if (Math.Abs(denseArray[i]) > tolerance) indicesValues.Add(i, denseArray[i]);
            }
            return CreateFromDictionary(denseArray.Length, indicesValues);
        }

        public static SparseVector CreateFromDense(Vector denseVector)
        {
            return CreateFromDense(denseVector.InternalData);
        }

        public static SparseVector CreateFromDense(Vector denseVector, double tolerance)
        {
            return CreateFromDense(denseVector.InternalData, tolerance);
        }

        public static SparseVector CreateFromDictionary(int length, Dictionary<int, double> indicesValues)
        {
            double[] values = new double[indicesValues.Count];
            int[] indices = new int[indicesValues.Count];
            int nnz = 0;
            foreach (var idxValPair in indicesValues.OrderBy(pair => pair.Key))
            {
                indices[nnz] = idxValPair.Key;
                values[nnz] = idxValPair.Value;
                ++nnz;
            }

            int maxIndex = indices[nnz - 1]; // first sort them
            if (maxIndex >= length)
            {
                throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
            }
            return new SparseVector(length, values, indices);
        }

        public static SparseVector CreateFromDictionary(int length, SortedDictionary<int, double> indicesValues)
        {
            double[] values = new double[indicesValues.Count];
            int[] indices = new int[indicesValues.Count];
            int nnz = 0;
            foreach (var idxValPair in indicesValues)
            {
                indices[nnz] = idxValPair.Key;
                values[nnz] = idxValPair.Value;
                ++nnz;
            }
            
            int maxIndex = indices[nnz - 1]; // first sort them
            if (maxIndex >= length)
            {
                throw new ArgumentException($"This sparse vector contains {maxIndex + 1} entries, while its length is {length}");
            }
            return new SparseVector(length, values, indices);
        }

        public IVectorView Axpy(IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is SparseVector otherSparse) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherSparse))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] result = new double[this.values.Length];
                    Array.Copy(this.values, result, this.values.Length);
                    CBlas.Daxpy(values.Length, otherCoefficient, ref otherSparse.values[0], 1, ref result[0], 1);
                    return new SparseVector(Length, result, indices);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this vector)
            return DenseStrategies.LinearCombination(this, 1.0, otherVector, otherCoefficient);
        }

        public SparseVector Copy()
        {
            int n = values.Length;
            double[] valuesCopy = new double[n];
            Array.Copy(values, valuesCopy, n);
            int[] indicesCopy = new int[n];
            Array.Copy(indices, indicesCopy, n);
            return new SparseVector(n, valuesCopy, indicesCopy);
        }

        public double[] CopyToArray()
        {
            double[] result = new double[Length];
            for (int i = 0; i < values.Length; ++i) result[indices[i]] = values[i];
            return result;
        }

        public Vector CopyToFullVector()
        {
            return Vector.CreateFromArray(CopyToArray(), false);
        }

        public int CountNonZeros()
        {
            return values.Length;
        }

        public IVectorView DoEntrywise(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            if (vector is SparseVector otherSparse) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherSparse))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherSparse.values[i]);
                    }
                    return new SparseVector(Length, resultValues, indices);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this vector)
            return DenseStrategies.DoEntrywise(this, vector, binaryOperation);
        }

        public IVectorView DoToAllEntries(Func<double, double> unaryOperation)
        {
            // Only apply the operation on non zero entries
            double[] newValues = new double[values.Length];
            for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
                int[] indicesCopy = new int[indices.Length];
                Array.Copy(indices, indicesCopy, indices.Length);
                return new SparseVector(Length, newValues, indicesCopy);
            }
            else // The sparsity is destroyed. Revert to a full vector.
            {
                return new SparseVector(Length, newValues, indices).CopyToFullVector();
            }
        }

        public double DotProduct(IVectorView vector)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double sum = 0;
            for (int i = 0; i < values.Length; ++i) sum += values[i] * vector[indices[i]];
            return sum;
        }

        public bool Equals(IIndexable1D other, double tolerance = 1e-13)
        {
            if (this.Length != other.Length) return false;
            var comparer = new ValueComparer(tolerance);
            int previousIndex = 0;
            for (int i = 0; i < indices.Length; ++i)
            {
                int index = indices[i];
                for (int j = previousIndex; j < index; ++j) // zero entries between the stored ones
                {
                    if (!comparer.AreEqual(0.0, other[j])) return false; 
                }
                if (!comparer.AreEqual(values[i], other[index])) return false; // Non zero entry
                previousIndex = index + 1;
            }
            return true;
        }

        public IEnumerable<(int idx, double value)> EnumerateNonZeros()
        {
            for (int i = 0; i < values.Length; ++i)
            {
                yield return (indices[i], values[i]);
            }
        }

        public IVectorView LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is SparseVector otherSparse) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherSparse))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] result = new double[this.values.Length];
                    Array.Copy(this.values, result, this.values.Length);
                    CBlas.Daxpby(values.Length, otherCoefficient, ref otherSparse.values[0], 1, 
                        thisCoefficient, ref result[0], 1);
                    return new SparseVector(Length, result, indices);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this vector)
            return DenseStrategies.LinearCombination(this, thisCoefficient, otherVector, otherCoefficient);
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
            return Array.BinarySearch<int>(indices, denseIdx); // only works if indices are sorted!!!
        }

        private bool HasSameIndexer(SparseVector other)
        {
            return this.indices == other.indices;
        }
    }
}
