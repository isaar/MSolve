using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;

//TODO: align data using mkl_malloc
//TODO: tensor product, vector2D, vector3D
//TODO: rename all slice operations to GetSubvector, SetSubvector, etc.
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public class Vector: IVectorView, ISliceable1D
    {
        private readonly double[] data;

        private Vector(double[] data)
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

        public static Vector CreateFromArray(double[] data, bool copyArray = false)
        {
            if (copyArray)
            {
                double[] clone = new double[data.Length];
                Array.Copy(data, clone, data.Length);
                return new Vector(clone);
            }
            else return new Vector(data);
        }

        /// <summary>
        /// The original vector will be copied.
        /// </summary>
        /// <param name="original"></param>
        /// <returns></returns>
        public static Vector CreateFromVector(IVectorView original)
        {
            if (original is Vector casted) return casted.Copy();
            double[] clone = new double[original.Length];
            for (int i = 0; i < clone.Length; ++i) clone[i] = original[i];
            return new Vector(clone);
        }

        public static Vector CreateWithValue(int length, double value)
        {
            double[] data = new double[length];
            for (int i = 0; i < length; ++i) data[i] = value;
            return new Vector(data);
        }

        public static Vector CreateZero(int length)
        {
            return new Vector(new double[length]);
        }

        #region operators 
        public static Vector operator +(Vector vector1, Vector vector2)
        {
            return vector1.Axpy(1.0, vector2);
        }

        public static Vector operator -(Vector vector1, Vector vector2)
        {
            return vector1.Axpy(-1.0, vector2); //The order is important
        }

        public static Vector operator *(double scalar, Vector vector)
        {
            return vector.Scale(scalar);
        }

        public static Vector operator *(Vector vector, double scalar)
        {
            return vector.Scale(scalar);
        }

        public static double operator *(Vector vector1, Vector vector2)
        {
            return vector1.DotProduct(vector2); //Perhaps call BLAS directly
        }
        #endregion

        public void AddIntoThis(Vector other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            CBlas.Daxpy(Length, 1.0, ref other.data[0], 1, ref this.data[0], 1);
        }

        public Vector Append(Vector last)
        {
            int n1 = this.data.Length;
            int n2 = last.data.Length;
            var result = new double[n1 + n2];
            Array.Copy(this.data, result, n1);
            Array.Copy(last.data, 0, result, n1, n2);
            return new Vector(result);
        }

        /// <summary>
        /// result = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        public Vector Axpy(double scalar, Vector other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpy(Length, scalar, ref other.data[0], 1, ref result[0], 1);
            return new Vector(result);
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        public void AxpyIntoThis(double scalar, Vector other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            CBlas.Daxpy(Length, scalar, ref other.data[0], 1, ref this.data[0], 1);
        }

        public Vector Copy()
        {
            //TODO: Perhaps this should use BLAS
            return Vector.CreateFromArray(data, true);
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
            if (other is Vector casted) return DoEntrywise(other, binaryOperation);
            else return other.DoEntrywise(this, binaryOperation);
        }

        public Vector DoEntrywise(Vector other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, other);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(this.data[i], other.data[i]);
            return new Vector(result);
        }

        public void DoEntrywiseIntoThis(Vector other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, other);
            for (int i = 0; i < data.Length; ++i) data[i] = binaryOperation(data[i], other.data[i]);
        }

        IVectorView IVectorView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public Vector DoToAllEntries(Func<double, double> unaryOperation)
        {
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = unaryOperation(data[i]);
            return new Vector(result);
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = unaryOperation(data[i]);
        }

        public double DotProduct(IVectorView other)
        {
            if (other is Vector)
            {
                Preconditions.CheckVectorDimensions(this, other);
                double[] rawDataOther = ((Vector)other).InternalData;
                return CBlas.Ddot(Length, ref this.data[0], 1, ref rawDataOther[0], 1);
            }
            else return other.DotProduct(this); // Let the more complex/efficient object operate.
        }

        public bool Equals(Vector other, ValueComparer comparer = null)
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
        public Vector LinearCombination(double thisScalar, double otherScalar, Vector otherVector)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpby(Length, otherScalar, ref otherVector.data[0], 1, thisScalar, ref result[0], 1);
            return new Vector(result);
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <returns></returns>
        public void LinearCombinationIntoThis(double thisScalar, double otherScalar, Vector otherVector)
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
        /// The entries of this <see cref="Vector"/> will be sorted.
        /// </summary>
        /// <returns></returns>
        public Vector[] RemoveDuplicatesFindMultiplicity()
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
            Vector[] singlesMultiplicityVectors = new Vector[2];

            singlesMultiplicityVectors[0] = Vector.CreateZero(data.Length - numberOfZeros);
            for (int i = 0; i < data.Length - numberOfZeros; i++)
            {
                singlesMultiplicityVectors[0][i] = singles[i];
            }

            singlesMultiplicityVectors[1] = Vector.CreateZero(data.Length - numberOfZeros);
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
        public Vector Reorder(IReadOnlyList<int> permutation, bool oldToNew)
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
            return new Vector(reordered);
        }

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public Vector Scale(double scalar)
        {
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Dscal(Length, scalar, ref result[0], 1);
            return new Vector(result);
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

        public void SetSubvector(Vector subvector, int destinationIndex)
        {
            if (destinationIndex + subvector.Length > this.Length) throw new NonMatchingDimensionsException(
                "The entries to set exceed this vector's length");
            Array.Copy(subvector.data, 0, this.data, destinationIndex, subvector.Length);
        }

        /// <summary>
        /// Returns a subvector containing only the entries at the provided indices
        /// </summary>
        /// <param name="indices">Indices of the entries to be returned. They must be 0 &lt; = i &lt; <see cref="Length"/>.</param>
        /// <returns></returns>
        public Vector Slice(int[] indices)
        {
            double[] subvector = new double[indices.Length];
            for (int i = 0; i < indices.Length; ++i) subvector[i] = data[indices[i]];
            return new Vector(subvector);
        }

        /// <summary>
        /// Returns a subvector containing the entries at the indices between the provided start (inclusive) and end (exclusive).
        /// </summary>
        /// <param name="startInclusive">The first index from which to copy entries.</param>
        /// <param name="endExclusive">The index after the last one until which to copy entries.</param>
        /// <returns></returns>
        public Vector Slice(int startInclusive, int endExclusive)
        {
            Preconditions.CheckIndex1D(this, startInclusive);
            Preconditions.CheckIndex1D(this, endExclusive - 1);
            if (endExclusive < startInclusive) throw new ArgumentException(
                $"Exclusive end = {endExclusive} must be >= inclusive start = {startInclusive}");

            int subvectorLength = endExclusive - startInclusive;
            double[] subvector = new double[subvectorLength];
            Array.Copy(this.data, startInclusive, subvector, 0, subvectorLength);
            return new Vector(subvector);
        }

        public void SubtractIntoThis(Vector other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            CBlas.Daxpy(Length, -1.0, ref other.data[0], 1, ref this.data[0], 1);
        }

        /// <summary>
        /// Doesn't copy anything. Remove this once the design is cleaned. 
        /// </summary>
        /// <returns></returns>
        public Numerical.LinearAlgebra.Interfaces.IVector ToLegacyVector()
        {
            return new Numerical.LinearAlgebra.Vector(data);
        }
    }
}
