using System;
using System.Collections.Generic;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Reduction;

//TODO: align data using mkl_malloc
//TODO: tensor product, vector2D, vector3D
//TODO: have an AxpyToSubvector.
//TODO: GetSubvector, SetSubvector, AddSubvector, etc. Group them together and have a common naming/param convention.
//TODO: remove legacy vector conversions
namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// General purpose vector class with more functionality than other vectors. No sparsity is assumed.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Vector : IVector, ISliceable1D
    {
        private readonly double[] data;

        private Vector(double[] data)
        {
            this.data = data;
            this.Length = data.Length;
        }

        /// <summary>
        /// See <see cref="IIndexable1D.Length"/>.
        /// </summary>
        public int Length { get; }

        /// <summary>
        /// The internal array that stores the entries of the vector. It should only be used for passing the raw array to linear 
        /// algebra libraries.
        /// </summary>
        internal double[] InternalData { get { return data; } }

        /// <summary>
        /// See <see cref="IIndexable1D.this[int]"/>.
        /// </summary>
        public double this[int index]
        {
            get { return data[index]; }
            set { data[index] = value; }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector"/> with <paramref name="data"/> or a clone as its internal array.
        /// </summary>
        /// <param name="data">The entries of the vector to create.</param>
        /// <param name="copyArray">If true, <paramref name="data"/> will be copied and the new <see cref="Vector"/> instance 
        ///     will have a reference to the copy, which is safer. If false, the new vector will have a reference to 
        ///     <paramref name="data"/> itself, which is faster.</param>
        /// <returns></returns>
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
        /// Initializes a new instance of <see cref="Vector"/> by copying the entries of an existing vector: 
        /// <paramref name="original"/>.
        /// </summary>
        /// <param name="original">The original vector to copy.</param>
        public static Vector CreateFromVector(IVectorView original)
        {
            if (original is Vector casted) return casted.Copy();
            double[] clone = new double[original.Length];
            for (int i = 0; i < clone.Length; ++i) clone[i] = original[i];
            return new Vector(clone);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector"/> with all entries being equal to <paramref name="value"/>.
        /// </summary>
        /// <param name="length">The number of entries of the new <see cref="Vector"/> instance.</param>
        /// <param name="value">The value that all entries of the new vector will be initialized to.</param>
        public static Vector CreateWithValue(int length, double value)
        {
            double[] data = new double[length];
            for (int i = 0; i < length; ++i) data[i] = value;
            return new Vector(data);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector"/> with all entries being equal to 0.
        /// </summary>
        /// <param name="length">The number of entries of the new <see cref="Vector"/> instance.</param>
        public static Vector CreateZero(int length) => new Vector(new double[length]);

        #region operators 
        /// <summary>
        /// Performs the operation: result[i] = <paramref name="vector1"/>[i] + <paramref name="vector2"/>[i], 
        /// for 0 &lt;= i &lt; <paramref name="vector1"/>.<see cref="Length"/> = <paramref name="vector2"/>.<see cref="Length"/>.
        /// The resulting entries are written to a new <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vector1">The first <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector2"/>.</param>
        /// <param name="vector2">The second <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="vector1"/> and <paramref name="vector2"/>
        ///     have different <see cref="Length"/>.</exception>
        public static Vector operator +(Vector vector1, Vector vector2) => vector1.Axpy(vector2, 1.0);

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="vector1"/>[i] - <paramref name="vector2"/>[i], 
        /// for 0 &lt;= i &lt; <paramref name="vector1"/>.<see cref="Length"/> = <paramref name="vector2"/>.<see cref="Length"/>.
        /// The resulting entries are written to a new <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vector1">The first <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector2"/>.</param>
        /// <param name="vector2">The second <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="vector1"/> and <paramref name="vector2"/>
        ///     have different <see cref="Length"/>.</exception>
        public static Vector operator -(Vector vector1, Vector vector2) => vector1.Axpy(vector2, -1.0);

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * <paramref name="vector"/>[i],
        /// for 0 &lt;= i &lt; <paramref name="vector"/>.<see cref="Length"/>. The resulting entries are written to a new 
        /// <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        /// <param name="vector">The vector to multiply.</param>
        public static Vector operator *(double scalar, Vector vector) => vector.Scale(scalar);

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * <paramref name="vector"/>[i],
        /// for 0 &lt;= i &lt; <paramref name="vector"/>.<see cref="Length"/>. The resulting entries are written to a new 
        /// <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vector">The vector to multiply.</param>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        public static Vector operator *(Vector vector, double scalar) => vector.Scale(scalar);

        /// <summary>
        /// Performs the operation: scalar = sum(<paramref name="vector1"/>[i] * <paramref name="vector2"/>[i]), 
        /// for 0 &lt;= i &lt; <paramref name="vector1"/>.<see cref="Length"/> = <paramref name="vector2"/>.<see cref="Length"/>.
        /// </summary>
        /// <param name="vector1">The first <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector2"/>.</param>
        /// <param name="vector2">The second <see cref="Vector"/> operand. It must have the same <see cref="Length"/> as 
        ///     <paramref name="vector1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="vector1"/> and <paramref name="vector2"/>
        ///     have different <see cref="Length"/>.</exception>
        public static double operator *(Vector vector1, Vector vector2) => vector1.DotProduct(vector2); //TODO: Perhaps call BLAS directly
        #endregion

        /// <summary>
        /// Performs the operation: this[<paramref name="destinationIndex"/> + i] = this[<paramref name="destinationIndex"/> + i]
        /// + <paramref name="subvector"/>[i], for 0 &lt;= i &lt; <paramref name="subvector"/>.<see cref="Length"/>.
        /// </summary>
        /// <param name="subvector">The vector that will be added to a part of this <see cref="Vector"/> instance.</param>
        /// <param name="destinationIndex">The index into this <see cref="Vector"/> instance, at which to start the vector 
        ///     addition. Constraints: 0 &lt;= <paramref name="destinationIndex"/>, <paramref name="destinationIndex"/> + 
        ///     <paramref name="subvector"/>.<see cref="Length"/> &lt;= this.<see cref="Length"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <paramref name="destinationIndex"/> violates the 
        ///     described constraint.</exception>
        public void AddSubvector(Vector subvector, int destinationIndex)
        {
            if (destinationIndex + subvector.Length > this.Length) throw new NonMatchingDimensionsException(
                "The entries to set exceed this vector's length");
            CBlas.Daxpy(subvector.Length, 1.0, ref subvector.data[0], 1, ref this.data[destinationIndex], 1);
        }

        /// <summary>
        /// Creates a new <see cref="Vector"/> that contains all entries of this followed by all entries of 
        /// <paramref name="last"/>.
        /// </summary>
        /// <param name="last">The vector whose entries will be appended after all entries of this vector.</param>
        public Vector Append(Vector last)
        {
            //TODO: Move this to an ArrayManipulations utility class.
            int n1 = this.data.Length;
            int n2 = last.data.Length;
            var result = new double[n1 + n2];
            Array.Copy(this.data, result, n1);
            Array.Copy(last.data, 0, result, n1, n2);
            return new Vector(result);
        }

        /// <summary>
        /// See <see cref="IVectorView.Axpy(IVectorView, double)"/>.
        /// </summary>
        public IVectorView Axpy(IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector casted) return Axpy(casted, otherCoefficient);
            else return otherVector.LinearCombination(otherCoefficient, this, 1.0); // To avoid accessing zero entries
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>:
        /// result[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i]. 
        /// The resulting vector is written to a new <see cref="Vector"/> and then returned.
        /// </summary>
        /// <param name="otherVector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="Length"/> than this.</exception>
        public Vector Axpy(Vector otherVector, double otherCoefficient)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpy(Length, otherCoefficient, ref otherVector.data[0], 1, ref result[0], 1);
            return new Vector(result);
        }

        /// <summary>
        /// See <see cref="IVector.AxpyIntoThis(IVectorView, double)"/>.
        /// </summary>
        public void AxpyIntoThis(IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector casted) AxpyIntoThis(casted, otherCoefficient);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                for (int i = 0; i < Length; ++i)
                {
                    this.data[i] += otherCoefficient * otherVector[i];
                }
            }
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>:
        /// this[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i]. 
        /// The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="otherVector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="Length"/> than this.</exception>
        public void AxpyIntoThis(Vector otherVector, double otherCoefficient)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            CBlas.Daxpy(Length, otherCoefficient, ref otherVector.data[0], 1, ref this.data[0], 1);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector"/> by copying the entries of this instance.
        /// </summary>
        public Vector Copy() => Vector.CreateFromArray(data, true); //TODO: Perhaps this should use BLAS

        /// <summary>
        /// See <see cref="IVectorView.CopyToArray"/>.
        /// </summary>
        public double[] CopyToArray() //Perhaps this should use BLAS
        {
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return clone;
        }

        /// <summary>
        /// Copies <paramref name="length"/> consecutive entries from this <see cref="Vector"/> to a System array starting from 
        /// the provided indices.
        /// </summary>
        /// <param name="sourceIndex">The index into this <see cref="Vector"/> where to start copying from.</param>
        /// <param name="destinationArray">The System array where entries of this vector will be copied to.</param>
        /// <param name="destinationIndex">The index into this <paramref name="destinationArray"/> where to start copying 
        ///     to.</param>
        /// <param name="length">The number of entries to copy.</param>
        public void CopyToArray(int sourceIndex, double[] destinationArray, int destinationIndex, int length)
        {
            Array.Copy(data, sourceIndex, destinationArray, destinationIndex, length);
        }

        /// <summary>
        /// Copies <paramref name="length"/> consecutive entries from <paramref name="sourceVector"/> to this 
        /// <see cref="Vector"/> starting from the provided indices.
        /// </summary>
        /// <param name="destinationIndex">The index into this <see cref="Vector"/> where to start copying to.</param>
        /// <param name="sourceVector">The vector containing the entries to be copied.</param>
        /// <param name="sourceIndex">The index into this <paramref name="sourceVector"/> where to start copying from.</param>
        /// <param name="length">The number of entries to copy.</param>
        public void CopyFromVector(int destinationIndex, IVectorView sourceVector, int sourceIndex, int length)
        {
            //TODO: Perhaps a syntax closer to Array, 
            // e.g. Vector.Copy(sourceVector, sourceIndex, destinationVector, destinationIndex, length)
            for (int i = 0; i < length; ++i) data[i + destinationIndex] = sourceVector[i + sourceIndex];
        }

        /// <summary>
        /// See <see cref="IVectorView.DoEntrywise(IVectorView, Func{double, double, double})"/>.
        /// </summary>
        public IVectorView DoEntrywise(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            if (vector is Vector casted) return DoEntrywise(vector, binaryOperation);
            else return vector.DoEntrywise(this, binaryOperation);
        }

        /// <summary>
        /// Performs a binary operation on each pair of entries, for 0 &lt;= i &lt; this.<see cref="Length"/>:
        /// result[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector is written to a new <see cref="Vector"/> and then returned.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="vector"/> has different 
        ///     <see cref="Length"/> than this.</exception>
        public Vector DoEntrywise(Vector vector, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(this.data[i], vector.data[i]);
            return new Vector(result);
        }

        /// <summary>
        /// See <see cref="IVector.DoEntrywiseIntoThis(IVectorView, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            if (vector is Vector casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else
            {
                Preconditions.CheckVectorDimensions(this, vector);
                for (int i = 0; i < data.Length; ++i) data[i] = binaryOperation(data[i], vector[i]);
            }
        }

        /// <summary>
        /// Performs a binary operation on each pair of entries, for 0 &lt;= i &lt; this.<see cref="Length"/>:
        /// result[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="vector"/> has different 
        ///     <see cref="Length"/> than this.</exception>
        public void DoEntrywiseIntoThis(Vector vector, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            for (int i = 0; i < data.Length; ++i) data[i] = binaryOperation(data[i], vector.data[i]);
        }

        /// <summary>
        /// See <see cref="IVectorView.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        IVectorView IVectorView.DoToAllEntries(Func<double, double> unaryOperation) => DoToAllEntries(unaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i]), for 0 &lt;=
        /// i &lt; this.<see cref="Length"/>. The resulting vector is written to a new <see cref="Vector"/> and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        public Vector DoToAllEntries(Func<double, double> unaryOperation)
        {
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = unaryOperation(data[i]);
            return new Vector(result);
        }

        /// <summary>
        /// See <see cref="IVector.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        void IVector.DoToAllEntriesIntoThis(Func<double, double> unaryOperation) => DoToAllEntriesIntoThis(unaryOperation);

        /// <summary>
        /// Performs a unary operation on each entry: this[i] = <paramref name="unaryOperation"/>(this[i]), for 0 &lt;= i
        /// &lt; this.<see cref="Length"/>. The resulting vector overwrites the entries of this <see cref= "Vector" /> instance.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = unaryOperation(data[i]);
        }

        /// <summary>
        /// See <see cref="IVectorView.DotProduct(IVectorView)"/>.
        /// </summary>
        public double DotProduct(IVectorView vector)
        {
            if (vector is Vector casted) return DotProduct(casted);
            else return vector.DotProduct(this); // Let the more complex/efficient object operate.
        }

        /// <summary>
        /// Calculates the dot (or inner/scalar) product of this vector with <paramref name="vector"/>:
        /// result = sum over all i of this[i] * <paramref name="vector"/>[i]).
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="Length"/> as this.</param>
        public double DotProduct(Vector vector)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            return CBlas.Ddot(Length, ref this.data[0], 1, ref vector.data[0], 1);
        }

        /// <summary>
        /// See <see cref="IIndexable1D.Equals(IIndexable1D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable1D other, double tolerance = 1e-13)
        {
            if (other is Vector casted)
            {
                if (this.Length != other.Length) return false;
                var comparer = new ValueComparer(tolerance);
                for (int i = 0; i < Length; ++i)
                {
                    if (!comparer.AreEqual(this.data[i], casted.data[i])) return false;
                }
                return true;
            }
            else return other.Equals(this, tolerance); // To avoid accessing zero entries
        }

        /// <summary>
        /// See <see cref="ISliceable1D.GetSubvector(int[])"/>.
        /// </summary>
        public Vector GetSubvector(int[] indices)
        {
            double[] subvector = new double[indices.Length];
            for (int i = 0; i < indices.Length; ++i) subvector[i] = data[indices[i]];
            return new Vector(subvector);
        }

        /// <summary>
        /// See <see cref="ISliceable1D.GetSubvector(int, int)"/>.
        /// </summary>
        public Vector GetSubvector(int startInclusive, int endExclusive)
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

        /// <summary>
        /// Returns true if this[i] &lt;= <paramref name="tolerance"/> for 0 &lt;= i &lt; this.<see cref="Length"/>. 
        /// Otherwise false is returned.
        /// </summary>
        /// <param name="tolerance">The tolerance under which a vector entry is considered to be 0. It can be set to 0, to check 
        ///     if the entries are exactly 0.</param>
        public bool IsZero(double tolerance) => DenseStrategies.IsZero(data, tolerance);

        /// <summary>
        /// See <see cref="IVectorView.LinearCombination(double, IVectorView, double)"/>.
        /// </summary>
        public IVectorView LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else return otherVector.LinearCombination(otherCoefficient, this, thisCoefficient); // To avoid accessing zero entries
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>:
        /// result[i] = <paramref name="thisCoefficient"/> * this[i] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherVector"/>[i]. The resulting vector is written to a new <see cref="Vector"/> and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Vector"/>.</param>
        /// <param name="otherVector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="Length"/> than this.</exception>
        public Vector LinearCombination(double thisCoefficient, Vector otherVector, double otherCoefficient)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpby(Length, otherCoefficient, ref otherVector.data[0], 1, thisCoefficient, ref result[0], 1);
            return new Vector(result);
        }

        /// <summary>
        /// See <see cref="IVector.LinearCombinationIntoThis(double, IVectorView, double)"/>.
        /// </summary>
        public void LinearCombinationIntoThis(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                for (int i = 0; i < Length; ++i)
                {
                    this.data[i] = thisCoefficient * this.data[i] + otherCoefficient * otherVector[i];
                }
            }
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>:
        /// this[i] = <paramref name="thisCoefficient"/> * this[i] + <paramref name="otherCoefficient"/> * 
        /// <paramref name="otherVector"/>[i]. The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Vector"/>.</param>
        /// <param name="otherVector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherVector"/> has different 
        ///     <see cref="Length"/> than this.</exception> 
        public void LinearCombinationIntoThis(double thisCoefficient, Vector otherVector, double otherCoefficient)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            CBlas.Daxpby(Length, otherCoefficient, ref otherVector.data[0], 1, thisCoefficient, ref this.data[0], 1);
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>: 
        /// result[i] = this[i] * <paramref name="vector"/>[i]. 
        /// The resulting vector is written to a new <see cref="Vector"/> and then returned.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        public Vector MultiplyEntrywise(Vector vector) //TODO: use MKL's vector math
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = this.data[i] * vector.data[i];
            return new Vector(result);
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>: 
        /// this[i] = this[i] * <paramref name="vector"/>[i]. 
        /// The resulting vector overwrites the entries of this <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vector">A vector with the same <see cref="Length"/> as this <see cref="Vector"/> instance.</param>
        public void MultiplyEntrywiseIntoThis(Vector vector) //TODO: use MKL's vector math
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) this.data[i] *= vector.data[i];
        }

        /// <summary>
        /// Calculates the Euclidian norm or 2-norm of this vector. For more see 
        /// https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm.
        /// </summary>
        public double Norm2() => CBlas.Dnrm2(Length, ref data[0], 1);

        /// <summary>
        /// This method is used to remove duplicate values of a Knot Value Vector and return the multiplicity up to
        /// the requested Knot. The multiplicity of a single Knot can be derived using the exported multiplicity vector. 
        /// The entries of this <see cref="Vector"/> will be sorted.
        /// </summary>
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

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue;
            for (int i = 0; i < data.Length; ++i) accumulator = processEntry(data[i], accumulator);
            // no zeros implied
            return finalize(accumulator);
        }

        /// <summary>
        /// Creates a new <see cref="Vector"/> that contains the entries of this <see cref="Vector"/> with a different order,
        /// which is specified by the provided <paramref name="permutation"/> and <paramref name="oldToNew"/>.
        /// </summary>
        /// <param name="permutation">An array that contains the indices of this <see cref="Vector"/> in a different 
        ///     order.</param>
        /// <param name="oldToNew">If true, reordered[<paramref name="permutation"/>[i]] = original[i]. If false, 
        ///     reordered[i] = original[<paramref name="permutation"/>[i]].</param>
        public Vector Reorder(IReadOnlyList<int> permutation, bool oldToNew)
        {
            //TODO: perhaps I should transfer this to a permutation matrix (implemented as a vector)
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
        /// See <see cref="IVectorView.Scale(double)"/>.
        /// </summary>
        IVectorView IVectorView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; this.<see cref="Length"/>: 
        /// result[i] = <paramref name="scalar"/> * this[i]. 
        /// The resulting vector is written to a new <see cref="Vector"/> and then returned.
        /// </summary>
        /// <param name="scalar">The scalar value that multiplies all entries of the vector.</param>
        public Vector Scale(double scalar)
        {
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Dscal(Length, scalar, ref result[0], 1);
            return new Vector(result);
        }

        /// <summary>
        /// See <see cref="IVector.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => CBlas.Dscal(Length, scalar, ref data[0], 1);

        /// <summary>
        /// Sets all entries of this vector to be equal to <paramref name="value"/>.
        /// </summary>
        /// <param name="value">The value that all entries of the this vector will be equal to.</param>
        public void SetAll(double value)
        {
            for (int i = 0; i < Length; ++i) data[i] = value;
        }

        /// <summary>
        /// See <see cref="IVector.SetEntryRespectingPattern(int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int index, double value) => data[index] = value;

        /// <summary>
        /// Performs the operation: this[<paramref name="destinationIndex"/> + i] = <paramref name="subvector"/>[i], for
        /// 0 &lt;= i &lt; <paramref name="subvector"/>.<see cref="Length"/>.
        /// </summary>
        /// <param name="subvector">The vector that will be copied to a part of this <see cref="Vector"/> instance.</param>
        /// <param name="destinationIndex">The index into this <see cref="Vector"/> instance, at which to start the vector 
        ///     copying. Constraints: 0 &lt;= <paramref name="destinationIndex"/>, <paramref name="destinationIndex"/> + 
        ///     <paramref name="subvector"/>.<see cref="Length"/> &lt;= this.<see cref="Length"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <paramref name="destinationIndex"/> violates the 
        ///     described constraint.</exception>
        public void SetSubvector(Vector subvector, int destinationIndex) //TODO: this has common functionality with CopyFromVector
        {
            if (destinationIndex + subvector.Length > this.Length) throw new NonMatchingDimensionsException(
                "The entries to set exceed this vector's length");
            Array.Copy(subvector.data, 0, this.data, destinationIndex, subvector.Length);
        }

        /// <summary>
        /// Creates a new instance of the legacy vector class <see cref="Numerical.LinearAlgebra.Vector"/> with the same internal
        /// array as this <see cref="Vector"/> instance. Doesn't copy anything.
        /// </summary>
        public Numerical.LinearAlgebra.Interfaces.IVector ToLegacyVector() => new Numerical.LinearAlgebra.Vector(data);
    }
}
