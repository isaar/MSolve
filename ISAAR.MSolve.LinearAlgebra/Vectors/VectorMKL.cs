using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;

//TODO: align data using mkl_malloc
//TODO: tensor product, vector2D, vector3D
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public class VectorMKL: IVectorView, ISliceable1D
    {
        private readonly double[] data;

        private VectorMKL(double[] data)
        {
            this.data = data;
            this.Length = data.Length;
        }

        public int Length { get; }

        /// <summary>
        /// TODO: make this package-private. It should only be used for passing raw arrays to linear algebra libraries
        /// </summary>
        internal double[] InternalData { get { return data; } }

        public double this[int i]
        {
            get { return data[i]; }
            set { data[i] = value; }
        }

        public static VectorMKL CreateFromArray(double[] data, bool copyArray = false)
        {
            if (copyArray)
            {
                double[] clone = new double[data.Length];
                Array.Copy(data, clone, data.Length);
                return new VectorMKL(clone);
            }
            else return new VectorMKL(data);
        }

        /// <summary>
        /// The original vector will be copied.
        /// </summary>
        /// <param name="original"></param>
        /// <returns></returns>
        public static VectorMKL CreateFromVector(IVectorView original)
        {
            if (original is VectorMKL casted) return CreateFromVector(casted);
            double[] clone = new double[original.Length];
            for (int i = 0; i < clone.Length; ++i) clone[i] = original[i];
            return new VectorMKL(clone);
        }

        public static VectorMKL CreateWithValue(int length, double value)
        {
            double[] data = new double[length];
            for (int i = 0; i < length; ++i) data[i] = value;
            return new VectorMKL(data);
        }

        public static VectorMKL CreateZero(int length)
        {
            return new VectorMKL(new double[length]);
        }

        #region operators 
        public static VectorMKL operator +(VectorMKL vector1, VectorMKL vector2)
        {
            return vector1.Axpy(1.0, vector2);
        }

        public static VectorMKL operator -(VectorMKL vector1, VectorMKL vector2)
        {
            return vector1.Axpy(-1.0, vector2); //The order is important
        }

        public static VectorMKL operator *(double scalar, VectorMKL vector)
        {
            return vector.Scale(scalar);
        }

        public static VectorMKL operator *(VectorMKL vector, double scalar)
        {
            return vector.Scale(scalar);
        }

        public static double operator *(VectorMKL vector1, VectorMKL vector2)
        {
            return vector1.DotProduct(vector2); //Perhaps call BLAS directly
        }
        #endregion

        /// <summary>
        /// result = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        public VectorMKL Axpy(double scalar, VectorMKL other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpy(Length, scalar, ref other.data[0], 1, ref result[0], 1);
            return new VectorMKL(result);
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        public void AxpyInPlace(double scalar, VectorMKL other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            CBlas.Daxpy(Length, scalar, ref other.data[0], 1, ref this.data[0], 1);
        }

        public VectorMKL Copy()
        {
            //TODO: Perhaps this should use BLAS
            return VectorMKL.CreateFromArray(data, true);
        }

        //Perhaps this should use BLAS
        public double[] CopyToArray()
        {
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return clone;
        }

        public void CopyToArray(int sourceIndex, double[] destinationArray, int destinationIndex, int length)
        {
            Array.Copy(data, sourceIndex, destinationArray, destinationIndex, length);
        }

        /// <summary>
        /// Copy length consecutive entries from sourceVector starting at sourceIndex, to this object starting at 
        /// destinationIndex.
        /// </summary>
        /// <param name="destinationIndex">Where to start writing in this object.</param>
        /// <param name="sourceVector">The vector containing the entries to be copied.</param>
        /// <param name="sourceIndex">Where to start reading in the source vector.</param>
        /// <param name="length">The number of entries to be copied.</param>
        public void CopyFromVector(int destinationIndex, IVectorView sourceVector, int sourceIndex, int length) 
        {
            //TODO: Perhaps a syntax closer to Array, 
            // e.g. Vector.Copy(sourceVector, sourceIndex, destinationVector, destinationIndex, length)
            for (int i = 0; i < length; ++i) data[i + destinationIndex] = sourceVector[i + sourceIndex];
        }

        public IVectorView DoEntrywise(IVectorView other, Func<double, double, double> binaryOperation)
        {
            if (other is VectorMKL casted) return DoEntrywise(other, binaryOperation);
            else return other.DoEntrywise(this, binaryOperation);
        }

        public VectorMKL DoEntrywise(VectorMKL other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, other);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(this.data[i], other.data[i]);
            return new VectorMKL(result);
        }

        public void DoEntrywiseIntoThis(VectorMKL other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, other);
            for (int i = 0; i < data.Length; ++i) data[i] = binaryOperation(data[i], other[i]);
        }

        IVectorView IVectorView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public VectorMKL DoToAllEntries(Func<double, double> unaryOperation)
        {
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = unaryOperation(data[i]);
            return new VectorMKL(result);
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = unaryOperation(data[i]);
        }

        public double DotProduct(IVectorView other)
        {
            if (other is VectorMKL)
            {
                Preconditions.CheckVectorDimensions(this, other);
                double[] rawDataOther = ((VectorMKL)other).InternalData;
                return CBlas.Ddot(Length, ref this.data[0], 1, ref rawDataOther[0], 1);
            }
            else return other.DotProduct(this); // Let the more complex/efficient object operate.
        }

        public bool Equals(VectorMKL other, ValueComparer comparer = null)
        {
            if (this.Length != other.Length) return false;
            if (comparer == null) comparer = new ValueComparer(1e-13);
            for (int i = 0; i < Length; ++i)
            {
                if (!comparer.AreEqual(this.data[i], other.data[i])) return false;
            }
            return true;
        }

        /// <summary>
        /// result = thisScalar * this + otherScalar * otherVector
        /// </summary>
        /// <returns></returns>
        public VectorMKL LinearCombination(double thisScalar, double otherScalar, VectorMKL otherVector)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpby(Length, otherScalar, ref otherVector.data[0], 1, thisScalar, ref result[0], 1);
            return new VectorMKL(result);
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <returns></returns>
        public void LinearCombinationInPlace(double thisScalar, double otherScalar, VectorMKL otherVector)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            CBlas.Daxpby(Length, otherScalar, ref otherVector.data[0], 1, thisScalar, ref this.data[0], 1);
        }

        public double Norm2()
        {
            return CBlas.Dnrm2(Length, ref data[0], 1);
        }

        /// <summary>
        /// This method is used to remove duplicate values of a Knot Value Vector and return the multiplicity up to
        /// the requested Knot. The multiplicity of a single Knot can be derived using the exported multiplicity vector. 
        /// The entries of this <see cref="VectorMKL"/> will be sorted.
        /// </summary>
        /// <returns></returns>
        public VectorMKL[] RemoveDuplicatesFindMultiplicity()
        {
            Array.Sort(data);
            HashSet<double> set = new HashSet<double>();
            int indexSingles = 0;
            double[] singles = new double[data.Length];

            int[] multiplicity = new int[data.Length];
            int counterMultiplicity = 0;

            for (int i = 0; i < data.Length; i++)
            {
                // If same integer is already present then add method will return
                // FALSE
                if (set.Add(data[i]) == true)
                {
                    singles[indexSingles] = data[i];

                    multiplicity[indexSingles] = counterMultiplicity;
                    indexSingles++;

                }
                else
                {
                    counterMultiplicity++;
                }
            }
            int numberOfZeros = 0;
            for (int i = data.Length - 1; i >= 0; i--)
            {
                if (singles[i] == 0)
                {
                    numberOfZeros++;
                }
                else
                {
                    break;
                }
            }
            VectorMKL[] singlesMultiplicityVectors = new VectorMKL[2];

            singlesMultiplicityVectors[0] = VectorMKL.CreateZero(data.Length - numberOfZeros);
            for (int i = 0; i < data.Length - numberOfZeros; i++)
            {
                singlesMultiplicityVectors[0][i] = singles[i];
            }

            singlesMultiplicityVectors[1] = VectorMKL.CreateZero(data.Length - numberOfZeros);
            for (int i = 0; i < data.Length - numberOfZeros; i++)
            {
                singlesMultiplicityVectors[1][i] = multiplicity[i];
            }

            return singlesMultiplicityVectors;
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue;
            for (int i = 0; i < data.Length; ++i) accumulator = processEntry(data[i], accumulator);
            // no zeros implied
            return finalize(accumulator);
        }

        //TODO: perhaps I should transfer this to a permutation matrix (implemented as a vector)
        public VectorMKL Reorder(IReadOnlyList<int> permutation, bool oldToNew)
        {
            if (permutation.Count != Length)
            {
                throw new NonMatchingDimensionsException($"This vector has length = {Length}, while the permutation vector has"
                    + $" {permutation.Count} entries");
            }
            double[] reordered = new double[Length];
            if (oldToNew)
            {
                for (int i = 0; i < Length; ++i) reordered[permutation[i]] = data[i];
            }
            else // TODO: can they be written as one in a smarter way?
            {
                for (int i = 0; i < Length; ++i) reordered[i] = data[permutation[i]];
            }
            return new VectorMKL(reordered);
        }

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public VectorMKL Scale(double scalar)
        {
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Dscal(Length, scalar, ref result[0], 1);
            return new VectorMKL(result);
        }

        /// <summary>
        /// this = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public void ScaleIntoThis(double scalar)
        {
            CBlas.Dscal(Length, scalar, ref data[0], 1);
        }

        public void SetAll(double value)
        {
            for (int i = 0; i < Length; ++i) data[i] = value;
        }

        /// <summary>
        /// Returns a subvector containing only the entries at the provided indices
        /// </summary>
        /// <param name="indices">Indices of the entries to be returned. They must be 0 &lt; = i &lt; <see cref="Length"/>.</param>
        /// <returns></returns>
        public VectorMKL Slice(int[] indices)
        {
            double[] subvector = new double[indices.Length];
            for (int i = 0; i < indices.Length; ++i) subvector[i] = data[indices[i]];
            return new VectorMKL(subvector);
        }

        /// <summary>
        /// Returns a subvector containing the entries at the indices between the provided start (inclusive) and end (exclusive).
        /// </summary>
        /// <param name="startInclusive">The first index from which to copy entries.</param>
        /// <param name="endExclusive">The index after the last one until which to copy entries.</param>
        /// <returns></returns>
        public VectorMKL Slice(int startInclusive, int endExclusive)
        {
            int newLength = endExclusive - startInclusive;
            double[] subvector = new double[newLength];
            for (int i = 0; i < newLength; ++i) subvector[i] = data[startInclusive + i];
            return new VectorMKL(subvector);
        }

        /// <summary>
        /// Doesn't copy anything. Remove this once the design is cleaned. 
        /// </summary>
        /// <returns></returns>
        public IVector ToLegacyVector()
        {
            return new Vector(data);
        }
    }
}
