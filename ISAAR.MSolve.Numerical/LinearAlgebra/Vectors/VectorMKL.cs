using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

//TODO: align data using mkl_malloc
//TODO: tensor product, vector2D, vector3D
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Vectors
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

        public static VectorMKL CreateFromArray(double[] data, bool copyArray = true)
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

        /// <summary>
        /// The original vector will be copied.
        /// </summary>
        /// <param name="original"></param>
        /// <returns></returns>
        public static VectorMKL CreateFromVector(VectorMKL original)
        {
            //TODO: Perhaps this should use BLAS
            //TODO: Perhaps it should be an instance method CopyToMatrix(). Or the instance method would return an interface.
            double[] data = original.data;
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
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

        public VectorMKL DoPointwise(IVectorView other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, other);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(data[i], other[i]);
            return new VectorMKL(result);
        }

        public void DoPointwiseInPlace(VectorMKL other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, other);
            for (int i = 0; i < data.Length; ++i) data[i] = binaryOperation(data[i], other[i]);
        }

        public VectorMKL DoToAllEntries(Func<double, double> unaryOperation)
        {
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = unaryOperation(data[i]);
            return new VectorMKL(result);
        }

        public void DoToAllEntriesInPlace(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = unaryOperation(data[i]);
        }

        public double DotProduct(IVectorView other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            if (other is VectorMKL)
            {
                double[] rawDataOther = ((VectorMKL)other).InternalData;
                return CBlas.Ddot(Length, ref this.data[0], 1, ref rawDataOther[0], 1);
            }
            else
            {
                double result = 0.0;
                for (int i = 0; i < data.Length; ++i) result += data[i] * other[i];
                return result;
            }
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

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue;
            for (int i = 0; i < data.Length; ++i) accumulator = processEntry(data[i], accumulator);
            // no zeros implied
            return finalize(accumulator);
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
        public IVectorOLD ToLegacyVector()
        {
            return new Vector(data);
        }
    }
}
