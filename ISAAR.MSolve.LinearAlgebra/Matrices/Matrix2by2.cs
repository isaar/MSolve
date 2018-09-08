using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// A matrix with 2 rows and 2 columns. Optimized version of <see cref="Matrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Matrix2by2: IMatrix
    {
        private readonly double[,] data;

        private Matrix2by2(double[,] data)
        {
            this.data = data;
        }

        /// <summary>
        /// Returns true if <see cref="NumRows"/> == <see cref="NumColumns"/>.
        /// </summary>
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
        /// The internal array that stores the entries of the matrix. It should only be used for passing the raw array to linear 
        /// algebra libraries.
        /// </summary>
        internal double[,] InternalData { get { return data; } }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        public double this[int rowIdx, int colIdx]
        {
            get { return data[rowIdx, colIdx]; }
            set { data[rowIdx, colIdx] = value; }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix2by2"/> that contains the provided entries: 
        /// {{<paramref name="entry00"/>, <paramref name="entry01"/>}, 
        ///  {<paramref name="entry10"/>, <paramref name="entry11"/>}}.
        /// </summary>
        public static Matrix2by2 Create(double entry00, double entry01, double entry10, double entry11)
            => new Matrix2by2(new double[,] { { entry00, entry01 }, { entry10, entry11 } });

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix2by2"/> with <paramref name="array2D"/> or a clone as its internal 
        /// array.
        /// </summary>
        /// <param name="array2D">A 2-dimensional array containing the entries of the matrix. Constraints 
        ///     <paramref name="array2D"/>.GetLength(0) ==  <paramref name="array2D"/>.GetLength(0) == 2.</param>
        /// <param name="copyArray">If true, <paramref name="array2D"/> will be copied and the new <see cref="Matrix2by2"/>  
        ///     instance will have a reference to the copy, which is safer. If false, the new matrix will have a reference to 
        ///     <paramref name="array2D"/> itself, which is faster.</param>
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

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix2by2"/> that is equal to the identity matrix, namely a square matrix 
        /// with non-diagonal entries being equal to 0 and diagonal entries being equal to 1.
        /// </summary>
        public static Matrix2by2 CreateIdentity() => new Matrix2by2(new double[,] { { 1.0, 0.0 }, { 0.0, 1.0 } });

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix2by2"/> with all entries being equal to <paramref name="value"/>.
        /// </summary>
        public static Matrix2by2 CreateWithValue(double value) 
            => new Matrix2by2(new double[,] { { value, value }, { value, value } });

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix2by2"/> with all entries being equal to 0.
        /// </summary> 
        public static Matrix2by2 CreateZero() => new Matrix2by2(new double[2, 2]);

        #region operators (use extension operators when they become available)
        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="m1"/>[i, j] + <paramref name="m2"/>[i, j], 
        /// for 0 &lt;= i, j &lt; 2. The resulting entries are written to a new <see cref="Matrix2by2"/> instance.
        /// </summary>
        /// <param name="m1">The first <see cref="Matrix2by2"/> operand.</param>
        /// <param name="m2">The second <see cref="Matrix2by2"/> operand.</param>
        public static Matrix2by2 operator +(Matrix2by2 m1, Matrix2by2 m2)
        {
            return new Matrix2by2(new double[,]
            {
                { m1.data[0, 0] + m2.data[0, 0], m1.data[0, 1] + m2.data[0, 1] }, 
                { m1.data[1, 0] + m2.data[1, 0], m1.data[1, 1] + m2.data[1, 1] }
            });
        }

        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="m1"/>[i, j] - <paramref name="m2"/>[i, j], 
        /// for 0 &lt;= i, j &lt; 2. The resulting entries are written to a new <see cref="Matrix2by2"/> instance.
        /// </summary>
        /// <param name="m1">The first <see cref="Matrix2by2"/> operand.</param>
        /// <param name="m2">The second <see cref="Matrix2by2"/> operand.</param>
        public static Matrix2by2 operator -(Matrix2by2 m1, Matrix2by2 m2)
        {
            return new Matrix2by2(new double[,]
            {
                { m1.data[0, 0] - m2.data[0, 0], m1.data[0, 1] - m2.data[0, 1] },
                { m1.data[1, 0] - m2.data[1, 0], m1.data[1, 1] - m2.data[1, 1] }
            });
        }

        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="scalar"/> * <paramref name="matrix"/>[i, j], 
        /// for 0 &lt;= i, j &lt; 2. The resulting entries are written to a new <see cref="Matrix2by2"/> instance.
        /// </summary>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        /// <param name="matrix">The matrix to multiply.</param>
        public static Matrix2by2 operator *(double scalar, Matrix2by2 matrix) => matrix.Scale(scalar);

        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="scalar"/> * <paramref name="matrix"/>[i, j], 
        /// for 0 &lt;= i, j &lt; 2. The resulting entries are written to a new <see cref="Matrix2by2"/> instance.
        /// </summary>
        /// <param name="matrix">The matrix to multiply.</param>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        public static Matrix2by2 operator *(Matrix2by2 matrix, double scalar) => matrix.Scale(scalar);

        /// <summary>
        /// Performs the matrix-matrix multiplication: result = <paramref name="matrixLeft"/> * <paramref name="matrixRight"/>.
        /// </summary>
        /// <param name="matrixLeft">The <see cref="Matrix2by2"/> operand on the left.</param>
        /// <param name="matrixRight">The <see cref="Matrix2by2"/> operand on the right.</param>
        public static Matrix2by2 operator *(Matrix2by2 matrixLeft, Matrix2by2 matrixRight)
            => matrixLeft.MultiplyRight(matrixRight, false, false);

        /// <summary>
        /// Performs the matrix-vector multiplication: result = <paramref name="matrixLeft"/> * <paramref name="vectorRight"/>.
        /// </summary>
        /// <param name="matrixLeft">The <see cref="Matrix2by2"/> operand on the left.</param>
        /// <param name="vectorRight">The <see cref="Vector2"/> operand on the right. It can be considered as a column 
        ///     vector.</param>
        public static Vector2 operator *(Matrix2by2 matrixLeft, Vector2 vectorRight)
            => matrixLeft.MultiplyRight(vectorRight, false);

        /// <summary>
        /// Performs the matrix-vector multiplication: result = <paramref name="vectorLeft"/> * <paramref name="matrixRight"/>.
        /// </summary>
        /// <param name="vectorLeft">The <see cref="Vector2"/> operand on the left. It can be considered as a row vector.</param>
        /// <param name="matrixRight">The <see cref="Matrix2by2"/> operand on the right.</param>
        public static Vector2 operator *(Vector2 vectorLeft, Matrix2by2 matrixRight)
            => matrixRight.MultiplyRight(vectorLeft, true);
        #endregion

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
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

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix is written to a new <see cref="Matrix2by2"/> and then returned.
        /// </summary>
        /// <param name="otherMatrix">A matrix with 2 rows and 2 columns.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
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

        /// <summary>
        /// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
        /// </summary>
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

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="Matrix2by2"/> instance.
        /// </summary>
        /// <param name="otherMatrix">A matrix with 2 rows and 2 columns.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        public void AxpyIntoThis(Matrix2by2 other, double otherCoefficient)
        {
            this.data[0, 0] += otherCoefficient * other.data[0, 0];
            this.data[0, 1] += otherCoefficient * other.data[0, 1];
            this.data[1, 0] += otherCoefficient * other.data[1, 0];
            this.data[1, 1] += otherCoefficient * other.data[1, 1];
        }

        /// <summary>
        /// Calculates the determinant of this matrix. If the inverse matrix is also needed, use
        /// <see cref="InvertAndDetermninant"/> instead.
        /// </summary>
        public double CalcDeterminant()
        {
            // Leibniz formula:
            return data[0, 0] * data[1, 1] - data[0, 1] * data[1, 0];
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix2by2"/> by copying the entries of this instance.
        /// </summary>
        public Matrix2by2 Copy() => new Matrix2by2(new double[,] { { data[0, 0], data[0, 1] }, { data[1, 0], data[1, 1] } });

        /// <summary>
        /// Copies the entries of the matrix into a 2-dimensional array. The returned array has length(0) = length(1) = 2. 
        /// </summary>
        public double[,] CopyToArray2D() => new double[,] { { data[0, 0], data[0, 1] }, { data[1, 0], data[1, 1] } };

        /// <summary>
        /// See <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public IMatrixView DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
        {
            if (matrix is Matrix2by2 casted) return DoEntrywise(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, matrix);
                return new Matrix2by2(new double[,]
                {
                    { binaryOperation(data[0, 0], matrix[0, 0]), binaryOperation(data[0, 1], matrix[0, 1]) },
                    { binaryOperation(data[1, 0], matrix[1, 0]), binaryOperation(data[1, 1], matrix[1, 1]) }
                });
            }
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// result[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i, j])
        /// The resulting matrix is written to a new <see cref="Matrix2by2"/> and then returned.
        /// </summary>
        /// <param name="matrix">A matrix with 2 rows and 2 columns.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        public Matrix2by2 DoEntrywise(Matrix2by2 matrix, Func<double, double, double> binaryOperation)
        {
            return new Matrix2by2(new double[,]
            {
                { binaryOperation(this.data[0, 0], matrix.data[0, 0]), binaryOperation(this.data[0, 1], matrix.data[0, 1]) },
                { binaryOperation(this.data[1, 0], matrix.data[1, 0]), binaryOperation(this.data[1, 1], matrix.data[1, 1]) }
            });
        }

        /// <summary>
        /// See <see cref="IMatrix.DoEntrywiseIntoThis(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation)
        {
            if (matrix is Matrix2by2 casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, matrix);
                data[0, 0] = binaryOperation(data[0, 0], matrix[0, 0]);
                data[0, 1] = binaryOperation(data[0, 1], matrix[0, 1]);
                data[1, 0] = binaryOperation(data[1, 0], matrix[1, 0]);
                data[1, 1] = binaryOperation(data[1, 1], matrix[1, 1]);
            }
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i, j])
        /// The resulting matrix overwrites the entries of this <see cref="Matrix2by2"/> instance.
        /// </summary>
        /// <param name="matrix">A matrix with 2 rows and 2 columns.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        public void DoEntrywiseIntoThis(Matrix2by2 matrix, Func<double, double, double> binaryOperation)
        {
            this.data[0, 0] = binaryOperation(this.data[0, 0], matrix.data[0, 0]);
            this.data[0, 1] = binaryOperation(this.data[0, 1], matrix.data[0, 1]);
            this.data[1, 0] = binaryOperation(this.data[1, 0], matrix.data[1, 0]);
            this.data[1, 1] = binaryOperation(this.data[1, 1], matrix.data[1, 1]);
        }

        /// <summary>
        /// See <see cref="IMatrixView.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation) => DoToAllEntries(unaryOperation);

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// result[i, j] = <paramref name="unaryOperation"/>(this[i,j])
        /// The resulting matrix is written to a new <see cref="Matrix2by2"/> and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        public Matrix2by2 DoToAllEntries(Func<double, double> unaryOperation)
        {
            return new Matrix2by2(new double[,]
            {
                { unaryOperation(data[0, 0]), unaryOperation(data[0, 1]) },
                { unaryOperation(data[1, 0]), unaryOperation(data[1, 1]) }
            });
        }

        /// <summary>
        /// See <see cref="IMatrix.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            data[0, 0] = unaryOperation(data[0, 0]); data[0, 1] = unaryOperation(data[0, 1]);
            data[1, 0] = unaryOperation(data[1, 0]); data[1, 1] = unaryOperation(data[1, 1]);
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
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

        /// <summary>
        /// Calculates the inverse matrix and returns it in a new <see cref="Matrix2by2"/> instance. This only works if this 
        /// <see cref="Matrix2by2"/> is invertible. If the determinant matrix is also needed, use 
        /// <see cref="InvertAndDetermninant"/> instead.
        /// </summary>
        /// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
        public Matrix2by2 Invert()
        {
            (Matrix2by2 inverse, double det) = InvertAndDetermninant();
            return inverse;
        }

        /// <summary>
        /// Calculates the determinant and the inverse matrix and returns the latter in a new <see cref="Matrix"/> instance. 
        /// This only works if this <see cref="Matrix2by2"/> is invertible.
        /// </summary>
        /// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
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

        /// <summary>
        /// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
        /// </summary>
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

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// result[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix is written to a new <see cref="Matrix2by2"/> and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix2by2"/>.</param>
        /// <param name="otherMatrix">A matrix with 2 rows and 2 columns.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
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

        /// <summary>
        /// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
        /// </summary>
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

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// result[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="Matrix2by2"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix2by2"/>.</param>
        /// <param name="otherMatrix">A matrix with 2 rows and 2 columns.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        public void LinearCombinationIntoThis(double thisCoefficient, Matrix2by2 otherMatrix, double otherCoefficient)
        {
            this.data[0, 0] = thisCoefficient * this.data[0, 0] + otherCoefficient * otherMatrix.data[0, 0];
            this.data[0, 1] = thisCoefficient * this.data[0, 1] + otherCoefficient * otherMatrix.data[0, 1];
            this.data[1, 0] = thisCoefficient * this.data[1, 0] + otherCoefficient * otherMatrix.data[1, 0];
            this.data[1, 1] = thisCoefficient * this.data[1, 1] + otherCoefficient * otherMatrix.data[1, 1];
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyLeft(IMatrixView matrix, bool transposeThis = false, bool transposeOther = false)
        {
            return matrix.MultiplyRight(this, transposeOther, transposeThis); //TODO: optimize this
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyRight(IMatrixView matrix, bool transposeThis = false, bool transposeOther = false)
        {
            // For now piggy back on Matrix. TODO: Optimize this
            double[] colMajor = new double[] { data[0, 0], data[1, 0], data[0, 1], data[1, 1] };
            return Matrix.CreateFromArray(colMajor, 2, 2, false).MultiplyRight(matrix, transposeThis, transposeOther);
        }

        /// <summary>
        /// Performs the matrix-matrix multiplication: oper(this) * oper(<paramref name="matrix"/>).
        /// </summary>
        /// <param name="matrix">A matrix with 2 rows and 2 columns.</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <param name="transposeOther">If true, oper(<paramref name="matrix"/>) = transpose(<paramref name="matrix"/>). 
        ///     Otherwise oper(<paramref name="matrix"/>) = <paramref name="matrix"/>.</param>
        public Matrix2by2 MultiplyRight(Matrix2by2 matrix, bool transposeThis = false, bool transposeOther = false)
        {
            double[,] result = new double[2, 2];
            if (transposeThis)
            {
                if (transposeOther)
                {
                    result[0, 0] = this.data[0, 0] * matrix.data[0, 0] + this.data[1, 0] * matrix.data[0, 1];
                    result[0, 1] = this.data[0, 0] * matrix.data[1, 0] + this.data[1, 0] * matrix.data[1, 1];
                    result[1, 0] = this.data[0, 1] * matrix.data[0, 0] + this.data[1, 1] * matrix.data[0, 1];
                    result[1, 1] = this.data[0, 1] * matrix.data[1, 0] + this.data[1, 1] * matrix.data[1, 1];
                }
                else
                {
                    result[0, 0] = this.data[0, 0] * matrix.data[0, 0] + this.data[1, 0] * matrix.data[1, 0];
                    result[0, 1] = this.data[0, 0] * matrix.data[0, 1] + this.data[1, 0] * matrix.data[1, 1];
                    result[1, 0] = this.data[0, 1] * matrix.data[0, 0] + this.data[1, 1] * matrix.data[1, 0];
                    result[1, 1] = this.data[0, 1] * matrix.data[0, 1] + this.data[1, 1] * matrix.data[1, 1];
                }
            }
            else
            {
                if (transposeOther)
                {
                    result[0, 0] = this.data[0, 0] * matrix.data[0, 0] + this.data[0, 1] * matrix.data[0, 1];
                    result[0, 1] = this.data[0, 0] * matrix.data[1, 0] + this.data[0, 1] * matrix.data[1, 1];
                    result[1, 0] = this.data[1, 0] * matrix.data[0, 0] + this.data[1, 1] * matrix.data[0, 1];
                    result[1, 1] = this.data[1, 0] * matrix.data[1, 0] + this.data[1, 1] * matrix.data[1, 1];
                }
                else
                {
                    result[0, 0] = this.data[0, 0] * matrix.data[0, 0] + this.data[0, 1] * matrix.data[1, 0];
                    result[0, 1] = this.data[0, 0] * matrix.data[0, 1] + this.data[0, 1] * matrix.data[1, 1];
                    result[1, 0] = this.data[1, 0] * matrix.data[0, 0] + this.data[1, 1] * matrix.data[1, 0];
                    result[1, 1] = this.data[1, 0] * matrix.data[0, 1] + this.data[1, 1] * matrix.data[1, 1];
                }
            }
            return new Matrix2by2(result);
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IVectorView, bool)"/>.
        /// </summary>
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
        /// Performs the matrix-vector multiplication: oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// </summary>
        /// <param name="vector">A vector with 2 entries.</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
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

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            double accumulator = identityValue; // no zeros implied
            accumulator = processEntry(data[0, 0], processEntry(data[0, 1], 
                processEntry(data[1, 0], processEntry(data[1, 1], accumulator))));
            return finalize(accumulator);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Scale(double)"/>.
        /// </summary>
        IMatrixView IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the following operation for 0 &lt;= i, j &lt; 2:
        /// result[i, j] = <paramref name="scalar"/> * <paramref name="this"/>[i, j].
        /// The resulting matrix is written to a new <see cref="Matrix2by2"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        public Matrix2by2 Scale(double scalar)
        {
            return new Matrix2by2(new double[,] 
            { 
                { scalar * data[0, 0], scalar * data[0, 1] }, 
                { scalar * data[1, 0], scalar * data[1, 1] }
            });
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar)
        {
            data[0, 0] *= scalar; data[0, 1] *= scalar;
            data[1, 0] *= scalar; data[1, 1] *= scalar;
        }

        /// <summary>
        /// Sets all entries of this matrix to be equal to <paramref name="value"/>.
        /// </summary>
        /// <param name="value">The value that all entries of the this matrix will be equal to.</param>
        public void SetAll(double value)
        {
            data[0, 0] = value; data[0, 1] = value;
            data[1, 0] = value; data[1, 1] = value;
        }

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value) => data[rowIdx, colIdx] = value;

        /// <summary>
        /// Creates a new instance of the legacy matrix class <see cref="Numerical.LinearAlgebra.Matrix"/>, by copying the 
        /// entries of this <see cref="Matrix2by2"/> instance. 
        /// </summary>
        public Numerical.LinearAlgebra.Interfaces.IMatrix2D ToLegacyMatrix() 
            => new Numerical.LinearAlgebra.Matrix2D(CopyToArray2D());

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        IMatrixView IMatrixView.Transpose() => Transpose();

        /// <summary>
        /// Initializes a new <see cref="Matrix2by2"/> instance, that is transpose to this: result[i, j] = this[j, i].  
        /// The entries will be explicitly copied.
        /// </summary>
        public Matrix2by2 Transpose() 
            => new Matrix2by2(new double[,] { { data[0, 0], data[1, 0] }, { data[0, 1], data[1, 1] } });
    }
}
