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
    public class Matrix2by2: IMatrix
    {
        private readonly double[,] data;

        private Matrix2by2(double[,] data)
        {
            this.data = data;
        }

        public bool IsSquare { get { return true; } }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get { return 2; } }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return 2; } }

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

        public static Matrix2by2 Create(double entry00, double entry01, double entry10, double entry11)
        {
            return new Matrix2by2(new double[,] { { entry00, entry01 }, { entry10, entry11 } });
        }

        /// <summary>
        /// Create a new <see cref="Matrix2by2"/> from a provided array. The array will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional array containing the elements of the matrix</param>
        /// <returns></returns>
        public static Matrix2by2 CreateFromArray(double[,] array2D, bool copyArray = false)
        {
            // TODO: more efficient checking
            if ((array2D.GetLength(0) != 2) || (array2D.GetLength(1) != 2)) throw new NonMatchingDimensionsException(
                $"The provided array was {array2D.GetLength(0)}-by{array2D.GetLength(1)} instead of 2-by-2");
            if (copyArray)
            {
                double[,] copy = new double[,] { { array2D[0, 0], array2D[0, 1] }, { array2D[1, 0], array2D[1, 1] } };
                return new Matrix2by2(copy);
            }
            else return new Matrix2by2(array2D);
        }

        public static Matrix2by2 CreateIdentity()
        {
            return new Matrix2by2(new double[,] { { 1.0, 0.0 }, { 0.0, 1.0 } });
        }

        public static Matrix2by2 CreateWithValue(double value)
        {
            return new Matrix2by2(new double[,] { { value, value }, { value, value } });
        }

        /// <summary>
        /// Create a new <see cref="Matrix2by2"/> with all entries equal to 0.
        /// </summary> 
        /// <param name="numRows">The number of rows of the matrix.</param>
        /// <param name="numColumns">The number of rows of the matrix.</param>
        /// <returns></returns>
        public static Matrix2by2 CreateZero()
        {
            return new Matrix2by2(new double[2, 2]);
        }

        #region operators (use extension operators when they become available)
        public static Matrix2by2 operator +(Matrix2by2 m1, Matrix2by2 m2)
        {
            return new Matrix2by2(new double[,]
            {
                { m1.data[0, 0] + m2.data[0, 0], m1.data[0, 1] + m2.data[0, 1] }, 
                { m1.data[1, 0] + m2.data[1, 0], m1.data[1, 1] + m2.data[1, 1] }
            });
        }

        public static Matrix2by2 operator -(Matrix2by2 m1, Matrix2by2 m2)
        {
            return new Matrix2by2(new double[,]
            {
                { m1.data[0, 0] - m2.data[0, 0], m1.data[0, 1] - m2.data[0, 1] },
                { m1.data[1, 0] - m2.data[1, 0], m1.data[1, 1] - m2.data[1, 1] }
            });
        }

        public static Matrix2by2 operator *(double scalar, Matrix2by2 matrix) => matrix.Scale(scalar);

        public static Matrix2by2 operator *(Matrix2by2 matrix, double scalar) => matrix.Scale(scalar);

        public static Matrix2by2 operator *(Matrix2by2 matrixLeft, Matrix2by2 matrixRight)
            => matrixLeft.MultiplyRight(matrixRight, false, false);

        public static Vector2 operator *(Matrix2by2 matrixLeft, Vector2 vectorRight)
            => matrixLeft.MultiplyRight(vectorRight, false);

        public static Vector2 operator *(Vector2 vectorLeft, Matrix2by2 matrixRight)
            => matrixRight.MultiplyRight(vectorLeft, true);
        #endregion

        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix2by2 casted) return Axpy(casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                return new Matrix2by2(new double[,]
                {
                    { data[0, 0] + otherCoefficient * otherMatrix[0, 0], data[0, 1] + otherCoefficient * otherMatrix[0, 1] },
                    { data[1, 0] + otherCoefficient * otherMatrix[1, 0], data[1, 1] + otherCoefficient * otherMatrix[1, 1] },
                });
            }
        }

        public Matrix2by2 Axpy(Matrix2by2 otherMatrix, double otherCoefficient)
        {
            return new Matrix2by2(new double[,]
            {
                {
                    this.data[0, 0] + otherCoefficient * otherMatrix.data[0, 0],
                    this.data[0, 1] + otherCoefficient * otherMatrix.data[0, 1]
                },
                {
                    this.data[1, 0] + otherCoefficient * otherMatrix.data[1, 0],
                    this.data[1, 1] + otherCoefficient * otherMatrix.data[1, 1]
                },
            });
        }

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix2by2 casted) AxpyIntoThis(casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                data[0, 0] += otherCoefficient * otherMatrix[0, 0];
                data[0, 1] += otherCoefficient * otherMatrix[0, 1];
                data[1, 0] += otherCoefficient * otherMatrix[1, 0];
                data[1, 1] += otherCoefficient * otherMatrix[1, 1];
            }
        }

        public void AxpyIntoThis(Matrix2by2 other, double otherCoefficient)
        {
            this.data[0, 0] += otherCoefficient * other.data[0, 0];
            this.data[0, 1] += otherCoefficient * other.data[0, 1];
            this.data[1, 0] += otherCoefficient * other.data[1, 0];
            this.data[1, 1] += otherCoefficient * other.data[1, 1];
        }

        public double CalcDeterminant()
        {
            // Leibniz formula:
            return data[0, 0] * data[1, 1] - data[0, 1] * data[1, 0];
        }

        public Matrix2by2 Copy()
        {
            return new Matrix2by2(new double[,] { { data[0, 0], data[0, 1] }, { data[1, 0], data[1, 1] } });
        }

        /// <summary>
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = 2 
        /// and length(1) = 2. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="NumRows"/>, <see cref="NumRows"/>] array 
        /// with the entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            return new double[,] { { data[0, 0], data[0, 1] }, { data[1, 0], data[1, 1] } };
        }

        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is Matrix2by2 casted) return DoEntrywise(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, other);
                return new Matrix2by2(new double[,]
                {
                    { binaryOperation(data[0, 0], other[0, 0]), binaryOperation(data[0, 1], other[0, 1]) },
                    { binaryOperation(data[1, 0], other[1, 0]), binaryOperation(data[1, 1], other[1, 1]) }
                });
            }
        }

        public Matrix2by2 DoEntrywise(Matrix2by2 other, Func<double, double, double> binaryOperation)
        {
            return new Matrix2by2(new double[,]
            {
                { binaryOperation(this.data[0, 0], other.data[0, 0]), binaryOperation(this.data[0, 1], other.data[0, 1]) },
                { binaryOperation(this.data[1, 0], other.data[1, 0]), binaryOperation(this.data[1, 1], other.data[1, 1]) }
            });
        }

        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is Matrix2by2 casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, other);
                data[0, 0] = binaryOperation(data[0, 0], other[0, 0]);
                data[0, 1] = binaryOperation(data[0, 1], other[0, 1]);
                data[1, 0] = binaryOperation(data[1, 0], other[1, 0]);
                data[1, 1] = binaryOperation(data[1, 1], other[1, 1]);
            }
        }

        public void DoEntrywiseIntoThis(Matrix2by2 other, Func<double, double, double> binaryOperation)
        {
            this.data[0, 0] = binaryOperation(this.data[0, 0], other.data[0, 0]);
            this.data[0, 1] = binaryOperation(this.data[0, 1], other.data[0, 1]);
            this.data[1, 0] = binaryOperation(this.data[1, 0], other.data[1, 0]);
            this.data[1, 1] = binaryOperation(this.data[1, 1], other.data[1, 1]);
        }

        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public Matrix2by2 DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new Matrix2by2(new double[,]
            {
                { unaryOperation(data[0, 0]), unaryOperation(data[0, 1]) },
                { unaryOperation(data[1, 0]), unaryOperation(data[1, 1]) }
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
            data[0, 0] = unaryOperation(data[0, 0]); data[0, 1] = unaryOperation(data[0, 1]);
            data[1, 0] = unaryOperation(data[1, 0]); data[1, 1] = unaryOperation(data[1, 1]);
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            var comparer = new ValueComparer(1e-13);
            if (other is Matrix2by2 casted)
            {
                return comparer.AreEqual(this.data[0, 0], casted.data[0, 0])
                    && comparer.AreEqual(this.data[0, 1], casted.data[0, 1])
                    && comparer.AreEqual(this.data[1, 0], casted.data[1, 0])
                    && comparer.AreEqual(this.data[1, 1], casted.data[1, 1]);
            }
            else
            {
                return comparer.AreEqual(this.data[0, 0], other[0, 0])
                    && comparer.AreEqual(this.data[0, 1], other[0, 1])
                    && comparer.AreEqual(this.data[1, 0], other[1, 0])
                    && comparer.AreEqual(this.data[1, 1], other[1, 1]);
            }            
        }

        public Matrix2by2 Invert()
        {
            (Matrix2by2 inverse, double det) = InvertAndDetermninant();
            return inverse;
        }

        public (Matrix2by2 inverse, double determinant) InvertAndDetermninant()
        {
            // Leibniz formula:
            double det = data[0, 0] * data[1, 1] - data[0, 1] * data[1, 0];
            if (Math.Abs(det) < AnalyticFormulas.determinantTolerance)
            {
                throw new SingularMatrixException($"Determinant == {det}");
            }

            // Cramer's rule: inverse = 1/det * [a11 -a01; -a10 a00]
            double[,] inverse = new double[,] {{ data[1, 1] / det, - data[0, 1] /det } , { -data[1, 0]/det, data[0, 0]/det }};
            return (new Matrix2by2(inverse), det);
        }

        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix2by2 casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                return new Matrix2by2(new double[,]
                {
                    {
                        thisCoefficient * data[0, 0] + otherCoefficient * otherMatrix[0, 0],
                        thisCoefficient * data[0, 1] + otherCoefficient * otherMatrix[0, 1]
                    },
                    {
                        thisCoefficient * data[1, 0] + otherCoefficient * otherMatrix[1, 0],
                        thisCoefficient * data[1, 1] + otherCoefficient * otherMatrix[1, 1]
                    },
                });
            }
        }

        public Matrix2by2 LinearCombination(double thisCoefficient, Matrix2by2 otherMatrix, double otherCoefficient)
        {
            return new Matrix2by2(new double[,]
            {
                {
                    thisCoefficient * this.data[0, 0] + otherCoefficient * otherMatrix.data[0, 0],
                    thisCoefficient * this.data[0, 1] + otherCoefficient * otherMatrix.data[0, 1]
                },
                {
                    thisCoefficient * this.data[1, 0] + otherCoefficient * otherMatrix.data[1, 0],
                    thisCoefficient * this.data[1, 1] + otherCoefficient * otherMatrix.data[1, 1]
                },
            });
        }

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix2by2 casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                data[0, 0] = thisCoefficient * data[0, 0] + otherCoefficient * otherMatrix[0, 0];
                data[0, 1] = thisCoefficient * data[0, 1] + otherCoefficient * otherMatrix[0, 1];
                data[1, 0] = thisCoefficient * data[1, 0] + otherCoefficient * otherMatrix[1, 0];
                data[1, 1] = thisCoefficient * data[1, 1] + otherCoefficient * otherMatrix[1, 1];
            }
        }

        public void LinearCombinationIntoThis(double thisCoefficient, Matrix2by2 otherMatrix, double otherCoefficient)
        {
            this.data[0, 0] = thisCoefficient * this.data[0, 0] + otherCoefficient * otherMatrix.data[0, 0];
            this.data[0, 1] = thisCoefficient * this.data[0, 1] + otherCoefficient * otherMatrix.data[0, 1];
            this.data[1, 0] = thisCoefficient * this.data[1, 0] + otherCoefficient * otherMatrix.data[1, 0];
            this.data[1, 1] = thisCoefficient * this.data[1, 1] + otherCoefficient * otherMatrix.data[1, 1];
        }

        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return other.MultiplyRight(this, transposeOther, transposeThis); //TODO: optimize this
        }

        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            // For now piggy back on Matrix. TODO: Optimize this
            double[] colMajor = new double[] { data[0, 0], data[1, 0], data[0, 1], data[1, 1] };
            return Matrix.CreateFromArray(colMajor, 2, 2, false).MultiplyRight(other, transposeThis, transposeOther);
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
        public Matrix2by2 MultiplyRight(Matrix2by2 other, bool transposeThis = false, bool transposeOther = false)
        {
            double[,] result = new double[2, 2];
            if (transposeThis)
            {
                if (transposeOther)
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[1, 0] * other.data[0, 1];
                    result[0, 1] = this.data[0, 0] * other.data[1, 0] + this.data[1, 0] * other.data[1, 1];
                    result[1, 0] = this.data[0, 1] * other.data[0, 0] + this.data[1, 1] * other.data[0, 1];
                    result[1, 1] = this.data[0, 1] * other.data[1, 0] + this.data[1, 1] * other.data[1, 1];
                }
                else
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[1, 0] * other.data[1, 0];
                    result[0, 1] = this.data[0, 0] * other.data[0, 1] + this.data[1, 0] * other.data[1, 1];
                    result[1, 0] = this.data[0, 1] * other.data[0, 0] + this.data[1, 1] * other.data[1, 0];
                    result[1, 1] = this.data[0, 1] * other.data[0, 1] + this.data[1, 1] * other.data[1, 1];
                }
            }
            else
            {
                if (transposeOther)
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[0, 1] * other.data[0, 1];
                    result[0, 1] = this.data[0, 0] * other.data[1, 0] + this.data[0, 1] * other.data[1, 1];
                    result[1, 0] = this.data[1, 0] * other.data[0, 0] + this.data[1, 1] * other.data[0, 1];
                    result[1, 1] = this.data[1, 0] * other.data[1, 0] + this.data[1, 1] * other.data[1, 1];
                }
                else
                {
                    result[0, 0] = this.data[0, 0] * other.data[0, 0] + this.data[0, 1] * other.data[1, 0];
                    result[0, 1] = this.data[0, 0] * other.data[0, 1] + this.data[0, 1] * other.data[1, 1];
                    result[1, 0] = this.data[1, 0] * other.data[0, 0] + this.data[1, 1] * other.data[1, 0];
                    result[1, 1] = this.data[1, 0] * other.data[0, 1] + this.data[1, 1] * other.data[1, 1];
                }
            }
            return new Matrix2by2(result);
        }

        public Vector MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            Preconditions.CheckMultiplicationDimensions(2, vector.Length);
            if (transposeThis)
            {
                return Vector.CreateFromArray(new double[] 
                {
                    data[0, 0] * vector[0] + data[1, 0] * vector[1],
                    data[0, 1] * vector[0] + data[1, 1] * vector[1]
                });
            }
            else
            {
                return Vector.CreateFromArray(new double[]
                {
                    data[0, 0] * vector[0] + data[0, 1] * vector[1],
                    data[1, 0] * vector[0] + data[1, 1] * vector[1]
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
        public Vector2 MultiplyRight(Vector2 vector, bool transposeThis = false)
        {
            if (transposeThis)
            {
                return Vector2.Create(data[0, 0] * vector[0] + data[1, 0] * vector[1],
                    data[0, 1] * vector[0] + data[1, 1] * vector[1]);
            }
            else
            {
                return Vector2.Create(data[0, 0] * vector[0] + data[0, 1] * vector[1],
                    data[1, 0] * vector[0] + data[1, 1] * vector[1]);
            }
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            double accumulator = identityValue; // no zeros implied
            accumulator = processEntry(data[0, 0], processEntry(data[0, 1], 
                processEntry(data[1, 0], processEntry(data[1, 1], accumulator))));
            return finalize(accumulator);
        }

        IMatrixView IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public Matrix2by2 Scale(double scalar)
        {
            return new Matrix2by2(new double[,] 
            { 
                { scalar * data[0, 0], scalar * data[0, 1] }, 
                { scalar * data[1, 0], scalar * data[1, 1] }
            });
        }

        /// <summary>
        /// this = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public void ScaleIntoThis(double scalar)
        {
            data[0, 0] *= scalar; data[0, 1] *= scalar;
            data[1, 0] *= scalar; data[1, 1] *= scalar;
        }

        public void SetAll(double value)
        {
            data[0, 0] = value; data[0, 1] = value;
            data[1, 0] = value; data[1, 1] = value;
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

        public Matrix2by2 Transpose()
        {
            return new Matrix2by2(new double[,] { { data[0, 0], data[1, 0] }, { data[0, 1], data[1, 1] } });
        }
    }
}
