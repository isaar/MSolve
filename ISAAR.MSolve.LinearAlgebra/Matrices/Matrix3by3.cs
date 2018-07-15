using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public class Matrix3by3 : IMatrix
    {
        private readonly double[,] data;

        private Matrix3by3(double[,] data)
        {
            this.data = data;
        }

        public bool IsSquare { get { return true; } }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get { return 3; } }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return 3; } }

        /// <summary>
        /// It should only be used for passing raw arrays to linear algebra libraries.
        /// </summary>
        internal double[,] InternalData { get { return data; } }

        /// <summary>
        /// The entry with row index = i and column index = j. 
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= i &lt; <see cref="NumRows"/></param>
        /// <param name="colIdx">The column index: 0 &lt;= j &lt; <see cref="NumColumns"/></param>
        /// <returns>The entry with indices i, j</returns>
        public double this[int rowIdx, int colIdx] //TODO: Should I add bound checking?
        {
            get { return data[rowIdx, colIdx]; }
            set { data[rowIdx, colIdx] = value; }
        }

        public static Matrix3by3 Create(double x00, double x01, double x02, double x10, double x11, double x12, 
            double x20, double x21, double x22)
        {
            return new Matrix3by3(new double[,] { { x00, x01, x02 }, { x10, x11, x12 }, { x20, x21, x22 } });
        }

        /// <summary>
        /// Create a new <see cref="Matrix3by3"/> from a provided array. The array will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional array containing the elements of the matrix</param>
        /// <returns></returns>
        public static Matrix3by3 CreateFromArray(double[,] array2D, bool copyArray)
        {
            // TODO: more efficient checking
            if ((array2D.GetLength(0) != 3) || (array2D.GetLength(1) != 3)) throw new NonMatchingDimensionsException(
                $"The provided array was {array2D.GetLength(0)}-by{array2D.GetLength(1)} instead of 2-by-2");
            if (copyArray)
            {
                double[,] copy = new double[,] 
                { 
                    { array2D[0, 0], array2D[0, 1], array2D[0, 2] }, 
                    { array2D[1, 0], array2D[1, 1], array2D[1, 2] },
                    { array2D[2, 0], array2D[2, 1], array2D[2, 2] }
                };
                return new Matrix3by3(copy);
            }
            else return new Matrix3by3(array2D);
        }

        public static Matrix3by3 CreateIdentity()
        {
            return new Matrix3by3(new double[,] { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } });
        }

        public static Matrix3by3 CreateWithValue(double value)
        {
            return new Matrix3by3(new double[,] { { value, value, value }, { value, value, value }, { value, value, value } });
        }

        /// <summary>
        /// Create a new <see cref="Matrix3by3"/> with all entries equal to 0.
        /// </summary> 
        /// <param name="numRows">The number of rows of the matrix.</param>
        /// <param name="numColumns">The number of rows of the matrix.</param>
        /// <returns></returns>
        public static Matrix3by3 CreateZero()
        {
            return new Matrix3by3(new double[3, 3]);
        }

        #region operators (use extension operators when they become available)
        public static Matrix3by3 operator +(Matrix3by3 m1, Matrix3by3 m2)
        {
            return new Matrix3by3(new double[,]
            {
                { m1.data[0, 0] + m2.data[0, 0], m1.data[0, 1] + m2.data[0, 1], m1.data[0, 2] + m2.data[0, 2] },
                { m1.data[1, 0] + m2.data[1, 0], m1.data[1, 1] + m2.data[1, 1], m1.data[1, 2] + m2.data[1, 2] },
                { m1.data[2, 0] + m2.data[2, 0], m1.data[2, 1] + m2.data[2, 1], m1.data[2, 2] + m2.data[2, 2] }
            });
        }

        public static Matrix3by3 operator -(Matrix3by3 m1, Matrix3by3 m2)
        {
            return new Matrix3by3(new double[,]
            {
                { m1.data[0, 0] - m2.data[0, 0], m1.data[0, 1] - m2.data[0, 1], m1.data[0, 2] - m2.data[0, 2] },
                { m1.data[1, 0] - m2.data[1, 0], m1.data[1, 1] - m2.data[1, 1], m1.data[1, 2] - m2.data[1, 2] },
                { m1.data[2, 0] - m2.data[2, 0], m1.data[2, 1] - m2.data[2, 1], m1.data[2, 2] - m2.data[2, 2] }
            });
        }

        public static Matrix3by3 operator *(double scalar, Matrix3by3 matrix) => matrix.Scale(scalar);

        public static Matrix3by3 operator *(Matrix3by3 matrix, double scalar) => matrix.Scale(scalar);

        public static Matrix3by3 operator *(Matrix3by3 matrixLeft, Matrix3by3 matrixRight)
            => matrixLeft.MultiplyRight(matrixRight, false, false);

        public static Vector3 operator *(Matrix3by3 matrixLeft, Vector3 vectorRight)
            => matrixLeft.MultiplyRight(vectorRight, false);

        public static Vector3 operator *(Vector3 vectorLeft, Matrix3by3 matrixRight)
            => matrixRight.MultiplyRight(vectorLeft, true);
        #endregion

        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix3by3 casted) return Axpy(casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                return new Matrix3by3(new double[,]
                {
                    {
                        data[0, 0] + otherCoefficient * otherMatrix[0, 0],
                        data[0, 1] + otherCoefficient * otherMatrix[0, 1],
                        data[0, 2] + otherCoefficient * otherMatrix[0, 2]
                    },
                    {
                        data[1, 0] + otherCoefficient * otherMatrix[1, 0],
                        data[1, 1] + otherCoefficient * otherMatrix[1, 1],
                        data[1, 2] + otherCoefficient * otherMatrix[1, 2]
                    },
                    {
                        data[2, 0] + otherCoefficient * otherMatrix[2, 0],
                        data[2, 1] + otherCoefficient * otherMatrix[2, 1],
                        data[2, 2] + otherCoefficient * otherMatrix[2, 2]
                    },
                });
            }
        }

        public Matrix3by3 Axpy(Matrix3by3 otherMatrix, double otherCoefficient)
        {
            return new Matrix3by3(new double[,]
            {
                {
                    this.data[0, 0] + otherCoefficient * otherMatrix.data[0, 0],
                    this.data[0, 1] + otherCoefficient * otherMatrix.data[0, 1],
                    this.data[0, 2] + otherCoefficient * otherMatrix.data[0, 2]
                },
                {
                    this.data[1, 0] + otherCoefficient * otherMatrix.data[1, 0],
                    this.data[1, 1] + otherCoefficient * otherMatrix.data[1, 1],
                    this.data[1, 2] + otherCoefficient * otherMatrix.data[1, 2]
                },
                {
                    this.data[2, 0] + otherCoefficient * otherMatrix.data[2, 0],
                    this.data[2, 1] + otherCoefficient * otherMatrix.data[2, 1],
                    this.data[2, 2] + otherCoefficient * otherMatrix.data[2, 2]
                }
            });
        }

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix3by3 casted) AxpyIntoThis(casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                data[0, 0] += otherCoefficient * otherMatrix[0, 0];
                data[0, 1] += otherCoefficient * otherMatrix[0, 1];
                data[0, 2] += otherCoefficient * otherMatrix[0, 2];
                data[1, 0] += otherCoefficient * otherMatrix[1, 0];
                data[1, 1] += otherCoefficient * otherMatrix[1, 1];
                data[1, 2] += otherCoefficient * otherMatrix[1, 2];
                data[2, 0] += otherCoefficient * otherMatrix[2, 0];
                data[2, 1] += otherCoefficient * otherMatrix[2, 1];
                data[2, 2] += otherCoefficient * otherMatrix[2, 2];
            }
        }

        public void AxpyIntoThis(Matrix3by3 other, double otherCoefficient)
        {
            this.data[0, 0] += otherCoefficient * other.data[0, 0];
            this.data[0, 1] += otherCoefficient * other.data[0, 1];
            this.data[0, 2] += otherCoefficient * other.data[0, 2];
            this.data[1, 0] += otherCoefficient * other.data[1, 0];
            this.data[1, 1] += otherCoefficient * other.data[1, 1];
            this.data[1, 2] += otherCoefficient * other.data[1, 2];
            this.data[2, 0] += otherCoefficient * other.data[2, 0];
            this.data[2, 1] += otherCoefficient * other.data[2, 1];
            this.data[2, 2] += otherCoefficient * other.data[2, 2];
        }

        public double CalcDeterminant()
        {
            // Laplace formula: 
            // det = a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20)  + a02 * (a10 * a21 - a11 * a20)
            return data[0, 0] * (data[1, 1] * data[2, 2] - data[1, 2] * data[2, 1])
                - data[0, 1] * (data[1, 0] * data[2, 2] - data[1, 2] * data[2, 0])
                + data[0, 2] * (data[1, 0] * data[2, 1] - data[1, 1] * data[2, 0]);
        }

        public Matrix3by3 Copy()
        {
            return new Matrix3by3(new double[,] 
            { 
                { data[0, 0], data[0, 1], data[0, 2] }, 
                { data[1, 0], data[1, 1], data[1, 2] }, 
                { data[2, 0], data[2, 1], data[2, 2] }
            });
        }

        /// <summary>
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = 3
        /// and length(1) = 3. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="NumRows"/>, <see cref="NumRows"/>] array 
        /// with the entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            return new double[,]
            {
                { data[0, 0], data[0, 1], data[0, 2] },
                { data[1, 0], data[1, 1], data[1, 2] },
                { data[2, 0], data[2, 1], data[2, 2] }
            };
        }

        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is Matrix3by3 casted) return DoEntrywise(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, other);
                return new Matrix3by3(new double[,]
                {
                    {
                        binaryOperation(data[0, 0], other[0, 0]),
                        binaryOperation(data[0, 1], other[0, 1]),
                        binaryOperation(data[0, 2], other[0, 2])
                    },
                    {
                        binaryOperation(data[1, 0], other[1, 0]),
                        binaryOperation(data[1, 1], other[1, 1]),
                        binaryOperation(data[1, 2], other[1, 2])
                    },
                    {
                        binaryOperation(data[2, 0], other[2, 0]),
                        binaryOperation(data[2, 1], other[2, 1]),
                        binaryOperation(data[2, 2], other[2, 2])
                    }
                });
            }
        }

        public Matrix3by3 DoEntrywise(Matrix3by3 other, Func<double, double, double> binaryOperation)
        {
            return new Matrix3by3(new double[,]
            {
                {
                    binaryOperation(this.data[0, 0], other.data[0, 0]),
                    binaryOperation(this.data[0, 1], other.data[0, 1]),
                    binaryOperation(this.data[0, 2], other.data[0, 2])
                },
                {
                    binaryOperation(this.data[1, 0], other.data[1, 0]),
                    binaryOperation(this.data[1, 1], other.data[1, 1]),
                    binaryOperation(this.data[1, 2], other.data[1, 2])
                },
                {
                    binaryOperation(this.data[2, 0], other.data[2, 0]),
                    binaryOperation(this.data[2, 1], other.data[2, 1]),
                    binaryOperation(this.data[2, 2], other.data[2, 2])
                }
            });
        }

        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is Matrix3by3 casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, other);
                data[0, 0] = binaryOperation(data[0, 0], other[0, 0]);
                data[0, 1] = binaryOperation(data[0, 1], other[0, 1]);
                data[0, 2] = binaryOperation(data[0, 2], other[0, 2]);
                data[1, 0] = binaryOperation(data[1, 0], other[1, 0]);
                data[1, 1] = binaryOperation(data[1, 1], other[1, 1]);
                data[1, 2] = binaryOperation(data[1, 2], other[1, 2]);
                data[2, 0] = binaryOperation(data[2, 0], other[2, 0]);
                data[2, 1] = binaryOperation(data[2, 1], other[2, 1]);
                data[2, 2] = binaryOperation(data[2, 2], other[2, 2]);
            }
        }

        public void DoEntrywiseIntoThis(Matrix3by3 other, Func<double, double, double> binaryOperation)
        {
            this.data[0, 0] = binaryOperation(this.data[0, 0], other.data[0, 0]);
            this.data[0, 1] = binaryOperation(this.data[0, 1], other.data[0, 1]);
            this.data[0, 2] = binaryOperation(this.data[0, 2], other.data[0, 2]);
            this.data[1, 0] = binaryOperation(this.data[1, 0], other.data[1, 0]);
            this.data[1, 1] = binaryOperation(this.data[1, 1], other.data[1, 1]);
            this.data[1, 2] = binaryOperation(this.data[1, 2], other.data[1, 2]);
            this.data[2, 0] = binaryOperation(this.data[2, 0], other.data[2, 0]);
            this.data[2, 1] = binaryOperation(this.data[2, 1], other.data[2, 1]);
            this.data[2, 2] = binaryOperation(this.data[2, 2], other.data[2, 2]);
        }

        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public Matrix3by3 DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new Matrix3by3(new double[,]
            {
                { unaryOperation(data[0, 0]), unaryOperation(data[0, 1]), unaryOperation(data[0, 2]) },
                { unaryOperation(data[1, 0]), unaryOperation(data[1, 1]), unaryOperation(data[1, 2]) },
                { unaryOperation(data[2, 0]), unaryOperation(data[2, 1]), unaryOperation(data[2, 2]) }
            });
        }

        void IMatrix.DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            DoToAllEntriesIntoThis(unaryOperation);
        }

        // Ok for a DenseMatrix, but for sparse formats some operation (e.g scale) maintain the sparsity pattern,
        // while others don't
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            data[0, 0] = unaryOperation(data[0, 0]);
            data[0, 1] = unaryOperation(data[0, 1]);
            data[0, 2] = unaryOperation(data[0, 2]);
            data[1, 0] = unaryOperation(data[1, 0]);
            data[1, 1] = unaryOperation(data[1, 1]);
            data[1, 2] = unaryOperation(data[1, 2]);
            data[2, 0] = unaryOperation(data[2, 0]);
            data[2, 1] = unaryOperation(data[2, 1]);
            data[2, 2] = unaryOperation(data[2, 2]);
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            var comparer = new ValueComparer(1e-13);
            if (other is Matrix3by3 casted)
            {
                return comparer.AreEqual(this.data[0, 0], casted.data[0, 0])
                    && comparer.AreEqual(this.data[0, 1], casted.data[0, 1])
                    && comparer.AreEqual(this.data[0, 2], casted.data[0, 2])
                    && comparer.AreEqual(this.data[1, 0], casted.data[1, 0])
                    && comparer.AreEqual(this.data[1, 1], casted.data[1, 1])
                    && comparer.AreEqual(this.data[1, 2], casted.data[1, 2])
                    && comparer.AreEqual(this.data[2, 0], casted.data[2, 0])
                    && comparer.AreEqual(this.data[2, 1], casted.data[2, 1])
                    && comparer.AreEqual(this.data[2, 2], casted.data[2, 2]);
            }
            else
            {
                return comparer.AreEqual(this.data[0, 0], other[0, 0])
                    && comparer.AreEqual(this.data[0, 1], other[0, 1])
                    && comparer.AreEqual(this.data[0, 2], other[0, 2])
                    && comparer.AreEqual(this.data[1, 0], other[1, 0])
                    && comparer.AreEqual(this.data[1, 1], other[1, 1])
                    && comparer.AreEqual(this.data[1, 2], other[1, 2])
                    && comparer.AreEqual(this.data[2, 0], other[2, 0])
                    && comparer.AreEqual(this.data[2, 1], other[2, 1])
                    && comparer.AreEqual(this.data[2, 2], other[2, 2]);
            }
        }

        public Matrix3by3 Invert()
        {
            (Matrix3by3 inverse, double det) = InvertAndDetermninant();
            return inverse;
        }

        public (Matrix3by3 inverse, double determinant) InvertAndDetermninant()
        {
            // Minors of the first row (with signs)
            double c00 = data[1, 1] * data[2, 2] - data[1, 2] * data[2, 1];     // c00 = + a11*a22 - a12*a21
            double c01 = -data[1, 0] * data[2, 2] + data[1, 2] * data[2, 0];    // c01 = - a10*a22 + a12*a20
            double c02 = data[1, 0] * data[2, 1] - data[1, 1] * data[2, 0];     // c02 = + a10*a21 - a11*a20
            double det = data[0, 0] * c00 + data[0, 1] * c01 + data[0, 2] * c02; // Laplace: det = a00*c00 + a01*c01 + a02*c02 
            if (Math.Abs(det) < AnalyticFormulas.determinantTolerance)
            {
                throw new SingularMatrixException($"Determinant == {det}");
            }

            // Cramer's rule: inverse = 1/det * transpose(C), C = matrix of minors
            double[,] inverse = new double[3, 3];
            inverse[0, 0] = c00 / det; // inv[0,0]: c10
            inverse[1, 0] = c01 / det; // inv[1,0]: c01
            inverse[2, 0] = c02 / det; // inv[2,0]: c02
            inverse[0, 1] = (-data[0, 1] * data[2, 2] + data[0, 2] * data[2, 1]) / det;    // inv[0,1]: c10 = - a01*a22 + a02*a21
            inverse[1, 1] = (data[0, 0] * data[2, 2] - data[0, 2] * data[2, 0]) / det;     // inv[1,1]: c11 = + a00*a22 - a02*a20
            inverse[2, 1] = (-data[0, 0] * data[2, 1] + data[0, 1] * data[2, 0]) / det;    // inv[2,1]: c12 = - a00*a21 + a01*a20
            inverse[0, 2] = (data[0, 1] * data[1, 2] - data[0, 2] * data[1, 1]) / det;     // inv[0,2]: c20 = + a01*a12 - a02*a11
            inverse[1, 2] = (-data[0, 0] * data[1, 2] + data[0, 2] * data[1, 2]) / det;    // inv[1,2]: c21 = - a00*a12 + a02*a12
            inverse[2, 2] = (data[0, 0] * data[1, 1] - data[0, 1] * data[1, 0]) / det;     // inv[2,2]: c22 = + a00*a11 - a01*a10

            return (new Matrix3by3(inverse), det);
        }

        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix3by3 casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                return new Matrix3by3(new double[,]
                {
                    {
                        thisCoefficient * data[0, 0] + otherCoefficient * otherMatrix[0, 0],
                        thisCoefficient * data[0, 1] + otherCoefficient * otherMatrix[0, 1],
                        thisCoefficient * data[0, 2] + otherCoefficient * otherMatrix[0, 2]

                    },
                    {
                        thisCoefficient * data[1, 0] + otherCoefficient * otherMatrix[1, 0],
                        thisCoefficient * data[1, 1] + otherCoefficient * otherMatrix[1, 1],
                        thisCoefficient * data[1, 2] + otherCoefficient * otherMatrix[1, 2]
                    },
                    {
                        thisCoefficient * data[2, 0] + otherCoefficient * otherMatrix[2, 0],
                        thisCoefficient * data[2, 1] + otherCoefficient * otherMatrix[2, 1],
                        thisCoefficient * data[2, 2] + otherCoefficient * otherMatrix[2, 2]
                    }
                });
            }
        }

        public Matrix3by3 LinearCombination(double thisCoefficient, Matrix3by3 otherMatrix, double otherCoefficient)
        {
            return new Matrix3by3(new double[,]
            {
                {
                    thisCoefficient * this.data[0, 0] + otherCoefficient * otherMatrix.data[0, 0],
                    thisCoefficient * this.data[0, 1] + otherCoefficient * otherMatrix.data[0, 1],
                    thisCoefficient * this.data[0, 2] + otherCoefficient * otherMatrix.data[0, 2]
                },
                {
                    thisCoefficient * this.data[1, 0] + otherCoefficient * otherMatrix.data[1, 0],
                    thisCoefficient * this.data[1, 1] + otherCoefficient * otherMatrix.data[1, 1],
                    thisCoefficient * this.data[1, 2] + otherCoefficient * otherMatrix.data[1, 2]
                },
                {
                    thisCoefficient * this.data[2, 0] + otherCoefficient * otherMatrix.data[2, 0],
                    thisCoefficient * this.data[2, 1] + otherCoefficient * otherMatrix.data[2, 1],
                    thisCoefficient * this.data[2, 2] + otherCoefficient * otherMatrix.data[2, 2]
                }
            });
        }

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix3by3 casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                data[0, 0] = thisCoefficient * data[0, 0] + otherCoefficient * otherMatrix[0, 0];
                data[0, 1] = thisCoefficient * data[0, 1] + otherCoefficient * otherMatrix[0, 1];
                data[0, 2] = thisCoefficient * data[0, 2] + otherCoefficient * otherMatrix[0, 2];
                data[1, 0] = thisCoefficient * data[1, 0] + otherCoefficient * otherMatrix[1, 0];
                data[1, 1] = thisCoefficient * data[1, 1] + otherCoefficient * otherMatrix[1, 1];
                data[1, 2] = thisCoefficient * data[1, 2] + otherCoefficient * otherMatrix[1, 2];
                data[2, 0] = thisCoefficient * data[2, 0] + otherCoefficient * otherMatrix[2, 0];
                data[2, 1] = thisCoefficient * data[2, 1] + otherCoefficient * otherMatrix[2, 1];
                data[2, 2] = thisCoefficient * data[2, 2] + otherCoefficient * otherMatrix[2, 2];
            }
        }

        public void LinearCombinationIntoThis(double thisCoefficient, Matrix3by3 otherMatrix, double otherCoefficient)
        {
            this.data[0, 0] = thisCoefficient * this.data[0, 0] + otherCoefficient * otherMatrix.data[0, 0];
            this.data[0, 1] = thisCoefficient * this.data[0, 1] + otherCoefficient * otherMatrix.data[0, 1];
            this.data[0, 2] = thisCoefficient * this.data[0, 2] + otherCoefficient * otherMatrix.data[0, 2];
            this.data[1, 0] = thisCoefficient * this.data[1, 0] + otherCoefficient * otherMatrix.data[1, 0];
            this.data[1, 1] = thisCoefficient * this.data[1, 1] + otherCoefficient * otherMatrix.data[1, 1];
            this.data[1, 2] = thisCoefficient * this.data[1, 2] + otherCoefficient * otherMatrix.data[1, 2];
            this.data[2, 0] = thisCoefficient * this.data[2, 0] + otherCoefficient * otherMatrix.data[2, 0];
            this.data[2, 1] = thisCoefficient * this.data[2, 1] + otherCoefficient * otherMatrix.data[2, 1];
            this.data[2, 2] = thisCoefficient * this.data[2, 2] + otherCoefficient * otherMatrix.data[2, 2];
        }

        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return other.MultiplyRight(this, transposeOther, transposeThis); //TODO: optimize this
        }

        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            // For now piggy back on Matrix. TODO: Optimize this
            double[] colMajor = new double[] {
                data[0, 0], data[1, 0], data[2, 0], data[0, 1], data[1, 1], data[2, 1], data[0, 2], data[1, 2], data[2, 2] };
            return Matrix.CreateFromArray(colMajor, 3, 3, false).MultiplyRight(other, transposeThis, transposeOther);
        }

        /// <summary>
        /// Matrix-matrix multiplication, with the other matrix on the right: this [2-by-2] * other [2-by-2] 
        /// or transpose(this [2-by-2]) * other [2-by-2].
        /// </summary>
        /// <param name="other">A matrix with as many rows as the column of this matrix.</param>
        /// <param name="transposeThis">Set to true to transpose this (the left matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <param name="transposeOther">Set to true to transpose other (the right matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <returns>A matrix with dimensions (m x n)</returns>
        public Matrix3by3 MultiplyRight(Matrix3by3 other, bool transposeThis = false, bool transposeOther = false)
        {
            double[,] result = new double[3, 3];
            if (transposeThis)
            {
                if (transposeOther)
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[1, 0] * other.data[0, 1] + this.data[2, 0] * other.data[0, 2];
                    result[0, 1] = this.data[0, 0] * other.data[1, 0] + this.data[1, 0] * other.data[1, 1] + this.data[2, 0] * other.data[1, 2];
                    result[0, 2] = this.data[0, 0] * other.data[2, 0] + this.data[1, 0] * other.data[2, 1] + this.data[2, 0] * other.data[2, 2];
                    result[1, 0] = this.data[0, 1] * other.data[0, 0] + this.data[1, 1] * other.data[0, 1] + this.data[2, 1] * other.data[0, 2];
                    result[1, 1] = this.data[0, 1] * other.data[1, 0] + this.data[1, 1] * other.data[1, 1] + this.data[2, 1] * other.data[1, 2];
                    result[1, 2] = this.data[0, 1] * other.data[2, 0] + this.data[1, 1] * other.data[2, 1] + this.data[2, 1] * other.data[2, 2];
                    result[2, 0] = this.data[0, 2] * other.data[0, 0] + this.data[1, 2] * other.data[0, 1] + this.data[2, 2] * other.data[0, 2];
                    result[2, 1] = this.data[0, 2] * other.data[1, 0] + this.data[1, 2] * other.data[1, 1] + this.data[2, 2] * other.data[1, 2];
                    result[2, 2] = this.data[0, 2] * other.data[2, 0] + this.data[1, 2] * other.data[2, 1] + this.data[2, 2] * other.data[2, 2];
                }
                else
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[1, 0] * other.data[1, 0] + this.data[2, 0] * other.data[2, 0];
                    result[0, 1] = this.data[0, 0] * other.data[0, 1] + this.data[1, 0] * other.data[1, 1] + this.data[2, 0] * other.data[2, 1];
                    result[0, 2] = this.data[0, 0] * other.data[0, 2] + this.data[1, 0] * other.data[1, 2] + this.data[2, 0] * other.data[2, 2];
                    result[1, 0] = this.data[0, 1] * other.data[0, 0] + this.data[1, 1] * other.data[1, 0] + this.data[2, 1] * other.data[2, 0];
                    result[1, 1] = this.data[0, 1] * other.data[0, 1] + this.data[1, 1] * other.data[1, 1] + this.data[2, 1] * other.data[2, 1];
                    result[1, 2] = this.data[0, 1] * other.data[0, 2] + this.data[1, 1] * other.data[1, 2] + this.data[2, 1] * other.data[2, 2];
                    result[2, 0] = this.data[0, 2] * other.data[0, 0] + this.data[2, 2] * other.data[1, 0] + this.data[2, 2] * other.data[2, 0];
                    result[2, 1] = this.data[0, 2] * other.data[0, 1] + this.data[2, 2] * other.data[1, 1] + this.data[2, 2] * other.data[2, 1];
                    result[2, 2] = this.data[0, 2] * other.data[0, 2] + this.data[2, 2] * other.data[1, 2] + this.data[2, 2] * other.data[2, 2];
                }
            }
            else
            {
                if (transposeOther)
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[0, 1] * other.data[0, 1] + this.data[0, 2] * other.data[0, 2];
                    result[0, 1] = this.data[0, 0] * other.data[1, 0] + this.data[0, 1] * other.data[1, 1] + this.data[0, 2] * other.data[1, 2];
                    result[0, 2] = this.data[0, 0] * other.data[2, 0] + this.data[0, 1] * other.data[2, 1] + this.data[0, 2] * other.data[2, 2];
                    result[1, 0] = this.data[1, 0] * other.data[0, 0] + this.data[1, 1] * other.data[0, 1] + this.data[1, 2] * other.data[0, 2];
                    result[1, 1] = this.data[1, 0] * other.data[1, 0] + this.data[1, 1] * other.data[1, 1] + this.data[1, 2] * other.data[1, 2];
                    result[1, 2] = this.data[1, 0] * other.data[2, 0] + this.data[1, 1] * other.data[2, 1] + this.data[1, 2] * other.data[2, 2];
                    result[2, 0] = this.data[2, 0] * other.data[0, 0] + this.data[2, 1] * other.data[0, 1] + this.data[2, 2] * other.data[0, 2];
                    result[2, 1] = this.data[2, 0] * other.data[1, 0] + this.data[2, 1] * other.data[1, 1] + this.data[2, 2] * other.data[1, 2];
                    result[2, 2] = this.data[2, 0] * other.data[2, 0] + this.data[2, 1] * other.data[2, 1] + this.data[2, 2] * other.data[2, 2];
                }
                else
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[0, 1] * other.data[1, 0] + this.data[0, 2] * other.data[2, 0];
                    result[0, 1] = this.data[0, 0] * other.data[0, 1] + this.data[0, 1] * other.data[1, 1] + this.data[0, 2] * other.data[2, 1];
                    result[0, 2] = this.data[0, 0] * other.data[0, 2] + this.data[0, 1] * other.data[1, 2] + this.data[0, 2] * other.data[2, 2];
                    result[1, 0] = this.data[1, 0] * other.data[0, 0] + this.data[1, 1] * other.data[1, 0] + this.data[1, 2] * other.data[2, 0];
                    result[1, 1] = this.data[1, 0] * other.data[0, 1] + this.data[1, 1] * other.data[1, 1] + this.data[1, 2] * other.data[2, 1];
                    result[1, 2] = this.data[1, 0] * other.data[0, 2] + this.data[1, 1] * other.data[1, 2] + this.data[1, 2] * other.data[2, 2];
                    result[2, 0] = this.data[2, 0] * other.data[0, 0] + this.data[2, 1] * other.data[1, 0] + this.data[2, 2] * other.data[2, 0];
                    result[2, 1] = this.data[2, 0] * other.data[0, 1] + this.data[2, 1] * other.data[1, 1] + this.data[2, 2] * other.data[2, 1];
                    result[2, 2] = this.data[2, 0] * other.data[0, 2] + this.data[2, 1] * other.data[1, 2] + this.data[2, 2] * other.data[2, 2];
                }
            }
            return new Matrix3by3(result);
        }

        public Vector MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            Preconditions.CheckMultiplicationDimensions(2, vector.Length);
            if (transposeThis)
            {
                return Vector.CreateFromArray(new double[]
                {
                    data[0, 0] * vector[0] + data[1, 0] * vector[1] + data[2, 0] * vector[2],
                    data[0, 1] * vector[0] + data[1, 1] * vector[1] + data[2, 1] * vector[2],
                    data[0, 2] * vector[0] + data[1, 2] * vector[1] + data[2, 2] * vector[2]
                });
            }
            else
            {
                return Vector.CreateFromArray(new double[]
                {
                    data[0, 0] * vector[0] + data[0, 1] * vector[1] + data[0, 2] * vector[2],
                    data[1, 0] * vector[0] + data[1, 1] * vector[1] + data[1, 2] * vector[2],
                    data[2, 0] * vector[0] + data[2, 1] * vector[1] + data[2, 2] * vector[2]
                });
            }
        }

        /// <summary>
        /// Matrix-vector multiplication, with the vector on the right: matrix * vector or transpose(matrix) * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to 2.</param>
        /// <param name="transposeThis">Set to true to transpose this (the left matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <returns></returns>
        public Vector3 MultiplyRight(Vector3 vector, bool transposeThis = false)
        {
            if (transposeThis)
            {
                return Vector3.Create(
                    data[0, 0] * vector[0] + data[1, 0] * vector[1] + data[2, 0] * vector[2],
                    data[0, 1] * vector[0] + data[1, 1] * vector[1] + data[2, 1] * vector[2],
                    data[0, 2] * vector[0] + data[1, 2] * vector[1] + data[2, 2] * vector[2]);
            }
            else
            {
                return Vector3.Create(
                    data[0, 0] * vector[0] + data[0, 1] * vector[1] + data[0, 2] * vector[2],
                    data[1, 0] * vector[0] + data[1, 1] * vector[1] + data[1, 2] * vector[2],
                    data[2, 0] * vector[0] + data[2, 1] * vector[1] + data[2, 2] * vector[2]);
            }
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            double accumulator = identityValue; // no zeros implied
            accumulator = processEntry(data[0, 0], processEntry(data[0, 1], processEntry(data[0, 2],
                processEntry(data[1, 0], processEntry(data[1, 1], processEntry(data[1, 2],
                processEntry(data[2, 0], processEntry(data[2, 1], processEntry(data[2, 2], accumulator)))))))));
            return finalize(accumulator);
        }

        IMatrixView IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public Matrix3by3 Scale(double scalar)
        {
            return new Matrix3by3(new double[,]
            {
                { scalar * data[0, 0], scalar * data[0, 1], scalar * data[0, 2] },
                { scalar * data[1, 0], scalar * data[1, 1], scalar * data[1, 2] },
                { scalar * data[2, 0], scalar * data[2, 1], scalar * data[2, 2] }
            });
        }

        /// <summary>
        /// this = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public void ScaleIntoThis(double scalar)
        {
            data[0, 0] *= scalar; data[0, 1] *= scalar; data[0, 2] *= scalar;
            data[1, 0] *= scalar; data[1, 1] *= scalar; data[1, 2] *= scalar;
            data[2, 0] *= scalar; data[2, 1] *= scalar; data[2, 2] *= scalar;
        }

        public void SetAll(double value)
        {
            data[0, 0] = value; data[0, 1] = value; data[0, 2] = value;
            data[1, 0] = value; data[1, 1] = value; data[1, 2] = value;
            data[2, 0] = value; data[2, 1] = value; data[2, 2] = value;
        }

        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            data[rowIdx, colIdx] = value;
        }

        /// <summary>
        /// Doesn't copy anything. Remove this once the design is cleaned. 
        /// </summary>
        /// <returns></returns>
        public Numerical.LinearAlgebra.Interfaces.IMatrix2D ToLegacyMatrix()
        {
            return new Numerical.LinearAlgebra.Matrix2D(CopyToArray2D());
        }

        IMatrixView IMatrixView.Transpose()
        {
            return Transpose();
        }

        public Matrix3by3 Transpose()
        {
            return new Matrix3by3(new double[,] 
            {
                { data[0, 0], data[1, 0], data[2, 0] }, 
                { data[0, 1], data[1, 1], data[2, 1] },
                { data[0, 2], data[1, 2], data[2, 2] }
            });
        }
    }
}
