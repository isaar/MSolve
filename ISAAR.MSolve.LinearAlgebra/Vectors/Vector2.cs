using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.LinearAlgebra.Vectors
{
    public class Vector2: IVectorView
    {
        private readonly double[] data;

        private Vector2(double[] data)
        {
            this.data = data;
        }

        public int Length { get { return 2; } }

        /// <summary>
        /// TODO: make this package-private. It should only be used for passing raw arrays to linear algebra libraries
        /// </summary>
        internal double[] InternalData { get { return data; } }

        public double this[int i]
        {
            get { return data[i]; }
            set { data[i] = value; }
        }

        public static Vector2 Create(double entry0, double entry1)
        {
            return new Vector2(new double[] { entry0, entry1 });
        }

        public static Vector2 CreateFromArray(double[] data, bool copyArray = false)
        {
            if (data.Length != 2) throw new NonMatchingDimensionsException(
                $"The provided array had length = {data.Length} instead of 2");
            if (copyArray) return new Vector2(new double[] { data[0], data[1] });
            else return new Vector2(data);
        }

        public static Vector2 CreateZero()
        {
            return new Vector2(new double[2]);
        }

        #region operators 
        public static Vector2 operator +(Vector2 v1, Vector2 v2)
        {
            return new Vector2(new double[] {v1.data[0] + v2.data[0], v1.data[1] + v2.data[1] });
        }

        public static Vector2 operator -(Vector2 v1, Vector2 v2)
        {
            return new Vector2(new double[] { v1.data[0] - v2.data[0], v1.data[1] - v2.data[1] });
        }

        public static Vector2 operator *(double scalar, Vector2 vector)
        {
            return new Vector2(new double[] { scalar * vector.data[0], scalar * vector.data[1] });
        }

        public static Vector2 operator *(Vector2 vector, double scalar)
        {
            return new Vector2(new double[] { scalar * vector.data[0], scalar * vector.data[1] });
        }

        public static double operator *(Vector2 v1, Vector2 v2)
        {
            return v1.data[0] * v2.data[0] + v1.data[1] * v2.data[1];
        }
        #endregion

        public void AddIntoThis(Vector2 other)
        {
            this.data[0] += other.data[0];
            this.data[1] += other.data[1];
        }

        /// <summary>
        /// result = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        public Vector2 Axpy(double scalar, Vector2 other)
        {
            return new Vector2(new double[] { this.data[0] + scalar * other.data[0], this.data[1] + scalar * other.data[1] });
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        public void AxpyIntoThis(double scalar, Vector2 other)
        {
            this.data[0] += scalar * other.data[0];
            this.data[1] += scalar * other.data[1];
        }

        public Vector2 Copy()
        {
            return new Vector2(new double[] { data[0], data[1] });
        }

        public double[] CopyToArray()
        {
            return new double[] { data[0], data[1] };
        }

        /// <summary>
        /// Defined as: cross = this[0] * other[1] - this[1] * other[0]. Also: other.Cross(this) = - this.Cross(other)
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public double CrossProduct(Vector2 other)
        {
            return this.data[0] * other.data[1] - this.data[1] * other.data[0];
        }

        public IVectorView DoEntrywise(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            if (vector is Vector2 casted) return DoEntrywise(vector, binaryOperation);
            else
            {
                Preconditions.CheckVectorDimensions(this, vector);
                return new Vector2(new double[] { binaryOperation(this.data[0], vector[0]),
                    binaryOperation(this.data[1], vector[1]) });
            }
        }

        public Vector2 DoEntrywise(Vector2 other, Func<double, double, double> binaryOperation)
        {
            return new Vector2(new double[] { binaryOperation(this.data[0], other.data[0]),
                binaryOperation(this.data[1], other.data[1]) });
        }

        public void DoEntrywiseIntoThis(Vector2 other, Func<double, double, double> binaryOperation)
        {
            this.data[0] = binaryOperation(this.data[0], other.data[0]);
            this.data[1] = binaryOperation(this.data[1], other.data[1]);
        }

        IVectorView IVectorView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public Vector2 DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new Vector2(new double[] { unaryOperation(data[0]), unaryOperation(data[1]) });
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            this.data[0] = unaryOperation(this.data[0]);
            this.data[1] = unaryOperation(this.data[1]);
        }

        public double DotProduct(IVectorView other)
        {
            if (other is Vector2 casted) return this.data[0] * casted.data[0] + this.data[1] * casted.data[1];
            else
            {
                {
                    Preconditions.CheckVectorDimensions(this, other);
                    return data[0] * other[0] + data[1] * other[1];
                }
            }
        }

        bool IIndexable1D.Equals(IIndexable1D other, double tolerance)
        {
            if (other.Length != 2) return false;
            else
            {
                var comparer = new ValueComparer(tolerance);
                return comparer.AreEqual(this.data[0], other[0]) && comparer.AreEqual(this.data[1], other[1]);
            }
        }

        public bool Equals(Vector2 other, double tolerance = 1e-13)
        {
            var comparer = new ValueComparer(tolerance);
            return comparer.AreEqual(this.data[0], other.data[0]) && comparer.AreEqual(this.data[1], other.data[1]);
        }

        /// <summary>
        /// result = thisScalar * this + otherScalar * otherVector2
        /// </summary>
        /// <returns></returns>
        public Vector2 LinearCombination(double thisScalar, double otherScalar, Vector2 otherVector2)
        {
            return new Vector2(new double[] { thisScalar * this.data[0] + otherScalar * otherVector2.data[0],
                thisScalar * this.data[1] + otherScalar * otherVector2.data[1] });
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <returns></returns>
        public void LinearCombinationIntoThis(double thisScalar, double otherScalar, Vector2 otherVector2)
        {
            this.data[0] = thisScalar * data[0] * otherScalar * otherVector2.data[0];
            this.data[1] = thisScalar * data[1] * otherScalar * otherVector2.data[1];
        }

        public double Norm2()
        {
            return Math.Sqrt(data[0] * data[0] + data[1] * data[1]);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue; // no zeros implied
            accumulator = processEntry(data[0], accumulator);
            accumulator = processEntry(data[1], accumulator);
            return finalize(accumulator);
        }

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public Vector2 Scale(double scalar)
        {
            return new Vector2(new double[] { scalar * data[0], scalar * data[1] });
        }

        /// <summary>
        /// this = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public void ScaleIntoThis(double scalar)
        {
            data[0] *= scalar;
            data[1] *= scalar;
        }

        public void SetAll(double value)
        {
            data[0] = value;
            data[1] = value;
        }

        public void SubtractIntoThis(Vector2 other)
        {
            this.data[0] -= other.data[0];
            this.data[1] -= other.data[1];
        }
    }
}
