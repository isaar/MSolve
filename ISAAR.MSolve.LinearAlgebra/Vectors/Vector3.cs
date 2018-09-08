using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Reduction;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    /// <summary>
    /// A vector with 3 entries. Optimized version of <see cref="Vector"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Vector3: IVectorView
    {
        private readonly double[] data;

        private Vector3(double[] data)
        {
            this.data = data;
        }

        /// <summary>
        /// See <see cref="IIndexable1D.Length"/>.
        /// </summary>
        public int Length { get { return 3; } }

        /// <summary>
        /// The internal array that stores the entries of the vector. It should only be used for passing the raw array to linear 
        /// algebra libraries.
        /// </summary>
        internal double[] InternalData { get { return data; } }

        /// <summary>
        /// See <see cref="IIndexable1D.this[int]"/>.
        /// </summary>
        public double this[int i]
        {
            get { return data[i]; }
            set { data[i] = value; }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector3"/> that contains the provided entries: 
        /// {<paramref name="entry0"/>, <paramref name="entry1"/>, <paramref name="entry2"/>}.
        /// </summary>
        /// <param name="entry0">The first entry of the new vector.</param>
        /// <param name="entry1">The second entry of the new vector.</param>
        /// <param name="entry2">The third entry of the new vector.</param>
        public static Vector3 Create(double entry0, double entry1, double entry2)
        {
            return new Vector3(new double[] { entry0, entry1, entry2 });
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector3"/> with <paramref name="data"/> or a clone as its internal array.
        /// </summary>
        /// <param name="data">The entries of the vector to create. Constraints: <paramref name="data"/>.Length == 3.</param>
        /// <param name="copyArray">If true, <paramref name="data"/> will be copied and the new <see cref="Vector3"/> instance 
        ///     will have a reference to the copy, which is safer. If false, the new vector will have a reference to 
        ///     <paramref name="data"/> itself, which is faster.</param>
        /// <returns></returns>
        public static Vector3 CreateFromArray(double[] data, bool copyArray = false)
        {
            if (data.Length != 3) throw new NonMatchingDimensionsException(
                $"The provided array had length = {data.Length} instead of 3");
            if (copyArray) return new Vector3(new double[] { data[0], data[1] });
            else return new Vector3(data);
        }

        /// <summary>
        /// Creates an new instance of <see cref="Vector3"/> with all entries being equal to 0.
        /// </summary>
        public static Vector3 CreateZero() => new Vector3(new double[3]);

        #region operators
        /// <summary>
        /// Performs the operation: result[i] = <paramref name="v1"/>[i] + <paramref name="v2"/>[i], 
        /// for 0 &lt;= i &lt; 3. The resulting entries are written to a new <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="v1">The first <see cref="Vector2"/> operand.</param>
        /// <param name="v2">The second <see cref="Vector2"/> operand.</param>
        public static Vector3 operator +(Vector3 v1, Vector3 v2)
        {
            return new Vector3(new double[] { v1.data[0] + v2.data[0], v1.data[1] + v2.data[1], v1.data[2] + v2.data[2] });
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="v1"/>[i] - <paramref name="v2"/>[i], 
        /// for 0 &lt;= i &lt; 3. The resulting entries are written to a new <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="v1">The first <see cref="Vector2"/> operand.</param>
        /// <param name="v2">The second <see cref="Vector2"/> operand.</param>
        public static Vector3 operator -(Vector3 v1, Vector3 v2)
        {
            return new Vector3(new double[] { v1.data[0] - v2.data[0], v1.data[1] - v2.data[1], v1.data[2] - v2.data[2] });
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * <paramref name="vector"/>[i],
        /// for 0 &lt;= i &lt; 3. The resulting entries are written to a new <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        /// <param name="vector">The vector to multiply.</param>
        public static Vector3 operator *(double scalar, Vector3 vector)
        {
            return new Vector3(new double[] { scalar * vector.data[0], scalar * vector.data[1], scalar * vector.data[2] });
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * <paramref name="vector"/>[i],
        /// for 0 &lt;= i &lt; 3. The resulting entries are written to a new <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="vector">The vector to multiply.</param>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        public static Vector3 operator *(Vector3 vector, double scalar)
        {
            return new Vector3(new double[] { scalar * vector.data[0], scalar * vector.data[1], scalar * vector.data[2] });
        }

        /// <summary>
        /// Performs the operation: scalar = sum(<paramref name="v1"/>[i] * <paramref name="v2"/>[i]), 
        /// for 0 &lt;= i &lt; 3. The resulting entries are written to a new <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="v1">The first <see cref="Vector2"/> operand.</param>
        /// <param name="v2">The second <see cref="Vector2"/> operand.</param>
        public static double operator *(Vector3 v1, Vector3 v2)
        {
            return v1.data[0] * v2.data[0] + v1.data[1] * v2.data[1] + v1.data[2] * v2.data[2];
        }
        #endregion

        /// <summary>
        /// Performs the operation: this[i] = this[i] + <paramref name="vector"/>[i], 
        /// for 0 &lt;= i &lt; 3. The resulting vector overwrites the entries of this <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="vector">A vector with three entries.</param>
        public void AddIntoThis(Vector3 vector)
        {
            this.data[0] += vector.data[0];
            this.data[1] += vector.data[1];
            this.data[2] += vector.data[2];
        }

        /// <summary>
        /// See <see cref="IVectorView.Axpy(IVectorView, double)"/>.
        /// </summary>
        public IVectorView Axpy(IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector3 casted) return Axpy(casted, otherCoefficient);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                return new Vector3(new double[]
                {
                    data[0] + otherCoefficient * otherVector[0],
                    data[1] + otherCoefficient * otherVector[1],
                    data[2] + otherCoefficient * otherVector[2]
                });
            }
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i], 
        /// for 0 &lt;= i &lt; 3. The resulting vector is written to a new <see cref="Vector3"/> and then returned.
        /// </summary>
        /// <param name="otherVector">A vector with three entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public Vector3 Axpy(Vector3 otherVector, double otherCoefficient)
        {
            return new Vector3(new double[] 
            {
                this.data[0] + otherCoefficient * otherVector.data[0],
                this.data[1] + otherCoefficient * otherVector.data[1],
                this.data[2] + otherCoefficient * otherVector.data[2]
            });
        }

        /// <summary>
        /// Performs the operation: this[i] = <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i] + this[i], 
        /// for 0 &lt;= i &lt; 3. The resulting vector overwrites the entries of this <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="otherVector">A vector with three entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public void AxpyIntoThis(Vector3 otherVector, double otherCoefficient)
        {
            this.data[0] += otherCoefficient * otherVector.data[0];
            this.data[1] += otherCoefficient * otherVector.data[1];
            this.data[2] += otherCoefficient * otherVector.data[2];
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Vector3"/> by deep copying the entries of this instance.
        /// </summary>
        public Vector3 Copy() => new Vector3(new double[] { data[0], data[1], data[2] });

        /// <summary>
        /// See <see cref="IVectorView.CopyToArray"/>.
        /// </summary>
        public double[] CopyToArray() => new double[] { data[0], data[1], data[2] };

        /// <summary>
        /// Performs the operation: result = { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}.
        /// The result is a vector. Also note that: other.Cross(this) = - this.Cross(other).
        /// </summary>
        /// <param name="vector">A vector with three entries.</param>
        public Vector3 CrossProduct(Vector3 vector)
        {
            return new Vector3(new double[]
            {
                this.data[1] * vector.data[2] - this.data[2] * vector.data[1],
                this.data[2] * vector.data[0] - this.data[0] * vector.data[2],
                this.data[0] * vector.data[1] - this.data[1] * vector.data[0]
            });
        }

        /// <summary>
        /// See <see cref="IVectorView.DoEntrywise(IVectorView, Func{double, double, double})"/>.
        /// </summary>
        public IVectorView DoEntrywise(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            if (vector is Vector3 casted) return DoEntrywise(vector, binaryOperation);
            else
            {
                Preconditions.CheckVectorDimensions(this, vector);
                return new Vector3(new double[] { binaryOperation(this.data[0], vector[0]),
                    binaryOperation(this.data[1], vector[1]), binaryOperation(this.data[2], vector[2]) });
            }
        }

        /// <summary>
        /// Performs a binary operation on each pair of entries, for 0 &lt;= i &lt; 3: 
        /// result[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector is written to a new <see cref="Vector3"/> and then returned.
        /// </summary>
        /// <param name="vector">A vector with three entries.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        public Vector3 DoEntrywise(Vector3 vector, Func<double, double, double> binaryOperation)
        {
            return new Vector3(new double[] { binaryOperation(this.data[0], vector.data[0]),
                binaryOperation(this.data[1], vector.data[1]), binaryOperation(this.data[2], vector.data[2]) });
        }

        /// <summary>
        /// Performs a binary operation on each pair of entries, for 0 &lt;= i &lt; 3: 
        /// result[i] = <paramref name="binaryOperation"/>(this[i], <paramref name="vector"/>[i]). 
        /// The resulting vector overwrites the entries of this <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="vector">A vector with three entries.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        public void DoEntrywiseIntoThis(Vector3 vector, Func<double, double, double> binaryOperation)
        {
            this.data[0] = binaryOperation(this.data[0], vector.data[0]);
            this.data[1] = binaryOperation(this.data[1], vector.data[1]);
            this.data[2] = binaryOperation(this.data[2], vector.data[2]);
        }

        /// <summary>
        /// See <see cref="IVectorView.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        IVectorView IVectorView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i]), 
        /// for 0 &lt;= i &lt; 3. The resulting vector is written to a new <see cref="Vector3"/> and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        public Vector3 DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new Vector3(new double[] { unaryOperation(data[0]), unaryOperation(data[1]), unaryOperation(data[2]) });
        }

        /// <summary>
        /// Performs a unary operation on each entry: result[i] = <paramref name="unaryOperation"/>(this[i]), 
        /// for 0 &lt;= i &lt; 3. The resulting vector overwrites the entries of this <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            this.data[0] = unaryOperation(this.data[0]);
            this.data[1] = unaryOperation(this.data[1]);
            this.data[2] = unaryOperation(this.data[2]);
        }

        /// <summary>
        /// See <see cref="IVectorView.DotProduct(IVectorView)"/>.
        /// </summary>
        public double DotProduct(IVectorView vector)
        {
            if (vector is Vector3 casted) return DotProduct(casted);
            else
            {
                Preconditions.CheckVectorDimensions(this, vector);
                return data[0] * vector[0] + data[1] * vector[1] + data[2] * vector[2];
            }
        }

        /// <summary>
        /// Calculates the dot (or inner/scalar) product of this vector with <paramref name="vector"/>:
        /// result = this[0] * <paramref name="vector"/>[0] + this[1] * <paramref name="vector"/>[1] 
        ///     + this[2] * <paramref name="vector"/>[2].
        /// </summary>
        /// <param name="vector">A vector with three entries</param>
        public double DotProduct(Vector3 vector)
        {
            return this.data[0] * vector.data[0] + this.data[1] * vector.data[1] + this.data[2] * vector.data[2];
        }

        /// <summary>
        /// See <see cref="IIndexable1D.Equals(IIndexable1D, double)"/>.
        /// </summary>
        bool IIndexable1D.Equals(IIndexable1D other, double tolerance)
        {
            if (other.Length != 3) return false;
            else
            {
                var comparer = new ValueComparer(tolerance);
                return comparer.AreEqual(this.data[0], other[0]) && comparer.AreEqual(this.data[1], other[1])
                    && comparer.AreEqual(this.data[2], other[2]);
            }
        }

        /// <summary>
        /// Returns true if this[i] - <paramref name="other"/>[i] is within the acceptable <paramref name="tolerance"/> for all
        /// 0 &lt;= i &lt; 3. 
        /// </summary>
        /// <param name="other">The other vector that this <see cref="Vector3"/> instance will be compared to.</param>
        /// <param name="tolerance">The entries at index i of the two vectors will be considered equal, if
        ///     (<paramref name="other"/>[i] - this[i]) / this[i] &lt;= <paramref name="tolerance"/>. Setting 
        ///     <paramref name="tolerance"/> = 0, will check if these entries are exactly the same.</param>
        public bool Equals(Vector3 other, double tolerance = 1e-13)
        {
            var comparer = new ValueComparer(tolerance);
            return comparer.AreEqual(this.data[0], other.data[0]) && comparer.AreEqual(this.data[1], other.data[1])
                && comparer.AreEqual(this.data[2], other.data[2]);
        }

        /// <summary>
        /// See <see cref="IVectorView.LinearCombination(double, IVectorView, double)"/>.
        /// </summary>
        public IVectorView LinearCombination(double thisCoefficient, IVectorView otherVector, double otherCoefficient)
        {
            if (otherVector is Vector3 casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else
            {
                Preconditions.CheckVectorDimensions(this, otherVector);
                return new Vector3(new double[]
                {
                    thisCoefficient * data[0] + otherCoefficient * otherVector[0],
                    thisCoefficient * data[1] + otherCoefficient * otherVector[1],
                    thisCoefficient * data[2] + otherCoefficient * otherVector[2]
                });
            }
        }

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="thisCoefficient"/> * this[i] + 
        ///     <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i], 
        /// for 0 &lt;= i &lt; 3. The resulting vector is written to a new <see cref="Vector3"/> and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this vector.</param>
        /// <param name="otherVector">A vector with three entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public Vector3 LinearCombination(double thisCoefficient, Vector3 otherVector, double otherCoefficient)
        {
            return new Vector3(new double[] { thisCoefficient * this.data[0] + otherCoefficient * otherVector.data[0],
                thisCoefficient * this.data[1] + otherCoefficient * otherVector.data[1],
                thisCoefficient * this.data[2] + otherCoefficient * otherVector.data[2] });
        }

        /// <summary>
        /// Performs the operation: this[i] = <paramref name="thisCoefficient"/> * this[i] + 
        ///     <paramref name="otherCoefficient"/> * <paramref name="otherVector"/>[i], 
        /// for 0 &lt;= i &lt; 3. The resulting vector overwrites the entries of this <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this vector.</param>
        /// <param name="otherVector">A vector with three entries.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherVector"/>.</param>
        public void LinearCombinationIntoThis(double thisCoefficient, Vector3 otherVector, double otherCoefficient)
        {
            this.data[0] = thisCoefficient * data[0] * otherCoefficient * otherVector.data[0];
            this.data[1] = thisCoefficient * data[1] * otherCoefficient * otherVector.data[1];
            this.data[2] = thisCoefficient * data[2] * otherCoefficient * otherVector.data[2];
        }

        /// <summary>
        /// Calculates the Euclidian norm or 2-norm of this vector. For more see 
        /// https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm.
        /// </summary>
        public double Norm2() => Math.Sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            // no zeros implied
            double accumulator = identityValue;
            accumulator = processEntry(data[0], accumulator);
            accumulator = processEntry(data[1], accumulator);
            accumulator = processEntry(data[2], accumulator);
            return finalize(accumulator);
        }

        /// <summary>
        /// See <see cref="IVectorView.Scale(double)"/>.
        /// </summary>
        IVectorView IVectorView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the operation: result[i] = <paramref name="scalar"/> * this[i],
        /// for 0 &lt;= i &lt; 23. The resulting vector is written to a new <see cref="Vector3"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        public Vector3 Scale(double scalar) => new Vector3(new double[] { scalar * data[0], scalar * data[1], scalar * data[2] });

        /// <summary>
        /// Performs the operation: this[i] = <paramref name="scalar"/> * this[i],
        /// for 0 &lt;= i &lt; 3. The resulting vector overwrites the entries of this <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this vector.</param>
        public void ScaleIntoThis(double scalar)
        {
            data[0] *= scalar;
            data[1] *= scalar;
            data[2] *= scalar;
        }

        ///<summary>
        /// Sets all entries of this vector to be equal to <paramref name="value"/>.
        /// </summary>
        /// <param name="value">The value that all entries of the this vector will be equal to.</param>
        public void SetAll(double value)
        {
            data[0] = value;
            data[1] = value;
            data[2] = value;
        }

        /// <summary>
        /// Performs the operation: this[i] = this[i] - <paramref name="vector"/>[i], 
        /// for 0 &lt;= i &lt; 3. The resulting vector overwrites the entries of this <see cref="Vector3"/> instance.
        /// </summary>
        /// <param name="vector">A vector with three entries.</param>
        public void SubtractIntoThis(Vector3 vector)
        {
            this.data[0] -= vector.data[0];
            this.data[1] -= vector.data[1];
            this.data[2] -= vector.data[3];
        }
    }
}
