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
    public class Vector3: IVectorView
    {
        private readonly double[] data;

        private Vector3(double[] data)
        {
            this.data = data;
        }

        public int Length { get { return 3; } }

        /// <summary>
        /// TODO: make this package-private. It should only be used for passing raw arrays to linear algebra libraries
        /// </summary>
        internal double[] InternalData { get { return data; } }

        public double this[int i]
        {
            get { return data[i]; }
            set { data[i] = value; }
        }

        public static Vector3 Create(double entry0, double entry1, double entry2)
        {
            return new Vector3(new double[] { entry0, entry1, entry2 });
        }

        public static Vector3 CreateFromArray(double[] data, bool copyArray = false)
        {
            if (data.Length != 3) throw new NonMatchingDimensionsException(
                $"The provided array had length = {data.Length} instead of 3");
            if (copyArray) return new Vector3(new double[] { data[0], data[1] });
            else return new Vector3(data);
        }

        public static Vector3 CreateZero()
        {
            return new Vector3(new double[3]);
        }

        #region operators 
        public static Vector3 operator +(Vector3 v1, Vector3 v2)
        {
            return new Vector3(new double[] { v1.data[0] + v2.data[0], v1.data[1] + v2.data[1], v1.data[2] + v2.data[2] });
        }

        public static Vector3 operator -(Vector3 v1, Vector3 v2)
        {
            return new Vector3(new double[] { v1.data[0] - v2.data[0], v1.data[1] - v2.data[1], v1.data[2] - v2.data[2] });
        }

        public static Vector3 operator *(double scalar, Vector3 vector)
        {
            return new Vector3(new double[] { scalar * vector.data[0], scalar * vector.data[1], scalar * vector.data[2] });
        }

        public static Vector3 operator *(Vector3 vector, double scalar)
        {
            return new Vector3(new double[] { scalar * vector.data[0], scalar * vector.data[1], scalar * vector.data[2] });
        }

        public static double operator *(Vector3 v1, Vector3 v2)
        {
            return v1.data[0] * v2.data[0] + v1.data[1] * v2.data[1] + v1.data[2] * v2.data[2];
        }
        #endregion

        public void AddIntoThis(Vector3 other)
        {
            this.data[0] += other.data[0];
            this.data[1] += other.data[1];
            this.data[2] += other.data[2];
        }

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
        /// result = this + scalar * other
        /// </summary>
        /// <param name="otherVector"></param>
        /// <param name="otherCoefficient"></param>
        /// <returns></returns>
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
        /// this = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="otherCoefficient"></param>
        public void AxpyIntoThis(Vector3 otherVector, double otherCoefficient)
        {
            this.data[0] += otherCoefficient * otherVector.data[0];
            this.data[1] += otherCoefficient * otherVector.data[1];
            this.data[2] += otherCoefficient * otherVector.data[2];
        }

        public Vector3 Copy()
        {
            return new Vector3(new double[] { data[0], data[1], data[2] });
        }

        public double[] CopyToArray()
        {
            return new double[] { data[0], data[1], data[2] };
        }

        /// <summary>
        /// Defined as: cross = { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}
        /// Also: other.Cross(this) = - this.Cross(other)
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public Vector3 CrossProduct(Vector3 vector)
        {
            return new Vector3(new double[]
            {
                this.data[1] * vector.data[2] - this.data[2] * vector.data[1],
                this.data[2] * vector.data[0] - this.data[0] * vector.data[2],
                this.data[0] * vector.data[1] - this.data[1] * vector.data[0]
            });
        }

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

        public Vector3 DoEntrywise(Vector3 vector, Func<double, double, double> binaryOperation)
        {
            return new Vector3(new double[] { binaryOperation(this.data[0], vector.data[0]),
                binaryOperation(this.data[1], vector.data[1]), binaryOperation(this.data[2], vector.data[2]) });
        }

        public void DoEntrywiseIntoThis(Vector3 vector, Func<double, double, double> binaryOperation)
        {
            this.data[0] = binaryOperation(this.data[0], vector.data[0]);
            this.data[1] = binaryOperation(this.data[1], vector.data[1]);
            this.data[2] = binaryOperation(this.data[2], vector.data[2]);
        }

        IVectorView IVectorView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public Vector3 DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new Vector3(new double[] { unaryOperation(data[0]), unaryOperation(data[1]), unaryOperation(data[2]) });
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            this.data[0] = unaryOperation(this.data[0]);
            this.data[1] = unaryOperation(this.data[1]);
            this.data[2] = unaryOperation(this.data[2]);
        }

        public double DotProduct(IVectorView vector)
        {
            if (vector is Vector3 casted) return DotProduct(casted);
            else
            {
                Preconditions.CheckVectorDimensions(this, vector);
                return data[0] * vector[0] + data[1] * vector[1] + data[2] * vector[2];
            }
        }

        public double DotProduct(Vector3 vector)
        {
            return this.data[0] * vector.data[0] + this.data[1] * vector.data[1] + this.data[2] * vector.data[2];
        }

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

        public bool Equals(Vector3 other, double tolerance = 1e-13)
        {
            var comparer = new ValueComparer(tolerance);
            return comparer.AreEqual(this.data[0], other.data[0]) && comparer.AreEqual(this.data[1], other.data[1])
                && comparer.AreEqual(this.data[2], other.data[2]);
        }

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
        /// result = thisScalar * this + otherScalar * otherVector3
        /// </summary>
        /// <returns></returns>
        public Vector3 LinearCombination(double thisCoefficient, Vector3 otherVector, double otherCoefficient)
        {
            return new Vector3(new double[] { thisCoefficient * this.data[0] + otherCoefficient * otherVector.data[0],
                thisCoefficient * this.data[1] + otherCoefficient * otherVector.data[1],
                thisCoefficient * this.data[2] + otherCoefficient * otherVector.data[2] });
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <returns></returns>
        public void LinearCombinationIntoThis(double thisCoefficient, Vector3 otherVector, double otherCoefficient)
        {
            this.data[0] = thisCoefficient * data[0] * otherCoefficient * otherVector.data[0];
            this.data[1] = thisCoefficient * data[1] * otherCoefficient * otherVector.data[1];
            this.data[2] = thisCoefficient * data[2] * otherCoefficient * otherVector.data[2];
        }

        public double Norm2()
        {
            return Math.Sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
        }

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
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public Vector3 Scale(double scalar)
        {
            return new Vector3(new double[] { scalar * data[0], scalar * data[1], scalar * data[2] });
        }

        /// <summary>
        /// this = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public void ScaleIntoThis(double scalar)
        {
            data[0] *= scalar;
            data[1] *= scalar;
            data[2] *= scalar;
        }

        public void SetAll(double value)
        {
            data[0] = value;
            data[1] = value;
            data[2] = value;
        }
    }
}
