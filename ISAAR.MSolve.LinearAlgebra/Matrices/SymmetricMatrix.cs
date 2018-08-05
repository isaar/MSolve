using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: align data using mkl_malloc
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Symmetric matrix. Only the upper triangle is stored in Packed format (only stores the n*(n+1)/2 non zeros) and column 
    /// major order. Uses Intel MKL. Do not use this, since it is an experimantal class, which will probably be removed.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SymmetricMatrix: IMatrix, ISymmetricMatrix
    {
        /// <summary>
        /// Packed storage, column major order, upper triangle: 
        /// A[i,j] = data[i + j*(j+1)/2] for 0 &lt;= i &lt;= j &lt; n.
        /// </summary>
        private readonly double[] data;

        private SymmetricMatrix(double[] data, int order, DefiniteProperty definiteness)
        {
            this.data = data;
            this.Definiteness = definiteness;
            this.Order = order;
            this.NumRows = order;
            this.NumColumns = order;
        }

        /// <summary>
        /// Used to query if the matrix is positive definite etc. Usually this is not known beforehand, which corresponds to
        /// <see cref="DefiniteProperty.Unknown"/>. Cholesky factorization reveals this and sets this property to
        /// <see cref="DefiniteProperty.PositiveDefinite"/> or <see cref="DefiniteProperty.Indefinite"/>. Mutating the matrix 
        /// will reset it to <see cref="DefiniteProperty.Unknown"/>. If the <see cref="SymmetricMatrix"/> is created from a 
        /// known matrix/array, the caller can assume responsibility for setting this property. WARNING: only set this propetry 
        /// if you are absolutely sure. 
        /// </summary>
        public DefiniteProperty Definiteness { get; set; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// The number of columns of the matrix.
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows or columns of the matrix.
        /// </summary>
        public int Order { get; }

        public int NumNonZeros => throw new NotImplementedException();

        /// <summary>
        /// The entry with row index = i and column index = j. Setting an entry A[i, j] = value, will also set A[j, i] = value. Therefore the matrix will stay symmetric 
        /// This property is not that efficient, due to the necessary bound checking.
        /// </summary>
        /// <param name="i">The row index: 0 &lt;= i &lt; <see cref="Order"/></param>
        /// <param name="j">The column index: 0 &lt;= j &lt; <see cref="Order"/></param>
        /// <returns>The entry with indices i, j</returns>
        public double this[int i, int j]
        {
            get //TODO: Perhaps keep the check in Debug mode only.
            {
                if ((i < 0) || (i >= Order) || (j < 0) || (j >= Order))
                {
                    throw new IndexOutOfRangeException($"Invalid indices: ({i}, {j})");
                }
                if (i <= j) return data[Find1DIndex(i, j)];
                else return data[Find1DIndex(j, i)];
            }
        }

        /// <summary>
        /// Create a new <see cref="SymmetricMatrix"/> from the lower (subdiagonal) or upper (superdiagonal) portion of the 
        /// provided array. The array entries will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional containing the elements of the whole matrix. Its lengths in both dimensions 
        ///     must be the same.</param>
        /// <param name="definiteness">If the caller knows that the matrix is positive definite, etc, he can set this property 
        ///     during creation of the <see cref="SymmetricMatrix"/> object.</param>
        /// <returns></returns>
        public static SymmetricMatrix CreateFromArray(double[,] array2D, 
            DefiniteProperty definiteness = DefiniteProperty.Unknown)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1);
            if (numRows != numCols)
            {
                string msg = string.Format("Provided array must have the same dimensions, but was ({0}x{1})", numRows, numCols);
                throw new NonMatchingDimensionsException(msg);
            }
            return new SymmetricMatrix(Conversions.Array2DToPackedUpperColMajor(array2D), numRows, definiteness);
        }

        /// <summary>
        /// Create a new <see cref="SymmetricMatrix"/> from a provided array. The array can be copied (for extra safety)
        /// or not (for extra performance).
        /// </summary>
        /// <param name="array1D">A 1-dimensional array containing the elements of the upper triangle of the matrix in column 
        ///     major order.</param>
        /// <param name="order"> The order of the matrix. It must be positive and match the length of <see cref="array1D"/>. If a 
        ///     value is provided, these will not be checked. If no value is provided, the order will be calculated from 
        ///     <see cref="array1D"/> instead.</param>
        /// <param name="definiteness">If the caller knows that the matrix is positive definite etc., he can set this property 
        ///     during creation of the <see cref="SymmetricMatrix"/> object.</param>
        /// <param name="copyArray">True to make a deep copy of <see cref="array1D"/>. 
        ///     False (default) to use <see cref="array1D"/> as its internal storage.</param>
        /// <returns></returns>
        public static SymmetricMatrix CreateFromArray(double[] array1D, int order = 0,
            DefiniteProperty definiteness = DefiniteProperty.Unknown, bool copyArray = false)
        {
            int n = (order == 0) ? Conversions.PackedLengthToOrder(array1D.Length) : order;
            if (copyArray)
            {
                var clone = new double[array1D.Length];
                Array.Copy(array1D, clone, array1D.Length);
                return new SymmetricMatrix(clone, n, definiteness);
            }
            else return new SymmetricMatrix(array1D, n, definiteness);
        }

        /// <summary>
        /// The caller is responsible for the original matrix being symmetric
        /// </summary>
        /// <param name="originalMatrix"></param>
        /// <returns></returns>
        public static SymmetricMatrix CreateFromMatrix(Matrix originalMatrix)
        {
            double[] data = Conversions.FullColMajorToPackedUpperColMajor(originalMatrix.InternalData, originalMatrix.NumColumns);
            return new SymmetricMatrix(data, originalMatrix.NumColumns, DefiniteProperty.Unknown);
        }

        /// <summary>
        /// Create a new <see cref="SymmetricMatrix"/> with the specified order and all entries equal to 0.
        /// </summary> 
        /// <param name="order">The number of rows or columns of the matrix.</param>
        /// <returns></returns>
        public static SymmetricMatrix CreateZero(int order)
        {
            double[] data = new double[order * order];
            //This matrix will be used as a canvas, thus we cannot infer that it is indefinite yet.
            return new SymmetricMatrix(data, order, DefiniteProperty.Unknown); 
        }

        #region operators (use extension operators when they become available)
        public static SymmetricMatrix operator +(SymmetricMatrix matrix1, SymmetricMatrix matrix2) 
            => matrix1.DoEntrywise(matrix2, (x, y) => x + y);

        public static SymmetricMatrix operator -(SymmetricMatrix matrix1, SymmetricMatrix matrix2) 
            => matrix1.DoEntrywise(matrix2, (x, y) => x - y);

        public static SymmetricMatrix operator *(double scalar, SymmetricMatrix matrix) 
            => matrix.DoToAllEntries(x => scalar * x);

        public static SymmetricMatrix operator *(SymmetricMatrix matrix, double scalar) 
            => matrix.DoToAllEntries(x => scalar * x);

        public static IMatrixView operator *(SymmetricMatrix matrixLeft, IMatrixView matrixRight)
            => matrixLeft.MultiplyRight(matrixRight, false, false);

        public static IMatrixView operator *(IMatrixView matrixLeft, SymmetricMatrix matrixRight)
            => matrixRight.MultiplyLeft(matrixLeft, false, false);

        public static Vector operator *(SymmetricMatrix matrixLeft, Vector vectorRight)
            => matrixLeft.MultiplyRight(vectorRight, false);

        public static Vector operator *(Vector vectorLeft, SymmetricMatrix matrixRight)
            => matrixRight.MultiplyRight(vectorLeft, true);

        #endregion

        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SymmetricMatrix casted) return Axpy(casted, otherCoefficient);
            else return DoEntrywise(otherMatrix, (x1, x2) => x1 + otherCoefficient * x2); //TODO: optimize this
        }

        public SymmetricMatrix Axpy(SymmetricMatrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref result[0], 1);
            return new SymmetricMatrix(result, NumColumns, DefiniteProperty.Unknown);
        }

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SymmetricMatrix casted) AxpyIntoThis(casted, otherCoefficient);
            else if (otherMatrix is ISymmetricMatrix otherSYM)
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                for (int j = 0; j < NumColumns; ++j)
                {
                    for (int i = 0; i <= j; ++i)
                    {
                        this.data[Find1DIndex(i, j)] += otherCoefficient * otherMatrix[i, j];
                    }
                }
            }
            else
            {
                throw new SymmetricPatternModifiedException("This operation is legal only if the other matrix is also symmetric.");
            }
        }

        public void AxpyIntoThis(SymmetricMatrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref this.data[0], 1);
            this.Definiteness = DefiniteProperty.Unknown;
        }

        /// <summary>
        /// Calculate the determinant of this matrix.
        /// If <see cref="Definiteness"/> != <see cref="DefiniteProperty.PositiveDefinite"/>, the calculation will be very 
        /// cumbersome and require an extra O(n^2) space: It will expand the packed storage format into a full data array, 
        /// factorize it using LU and then calculate the determinant. For positive definite matrices, it is more efficient,
        /// since it uses the Cholesky factorization directly.
        /// </summary>
        /// <returns></returns>
        public double CalcDeterminant()
        {
            if (Definiteness == DefiniteProperty.PositiveDefinite)
            {
                return FactorCholesky().CalcDeterminant();
            }
            else
            {
                return CopyToGeneralMatrix().FactorLU().CalcDeterminant(); //TODO: Find how to do it with Bunch-Kaufman
            }
        }

        public SymmetricMatrix Copy()
        {
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return new SymmetricMatrix(clone, Order, Definiteness);
        }

        /// <summary>
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="Order"/> 
        /// and length(1) = <see cref="Order"/>. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="Order"/>, <see cref="Order"/>] array 
        /// with the entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            return Conversions.PackedUpperColMajorToArray2DSymm(data, Order);
        }

        public Matrix CopyToGeneralMatrix()
        {
            double[] fullData = Conversions.PackedUpperColMajorToFullSymmColMajor(data, Order);
            return Matrix.CreateFromArray(fullData, Order, Order, false);
        }

        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is SymmetricMatrix casted) return DoEntrywise(casted, binaryOperation);
            else return DenseStrategies.DoEntrywise(this, other, binaryOperation); //TODO: optimize this
        }

        public SymmetricMatrix DoEntrywise(SymmetricMatrix other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(this.data[i], other.data[i]);
            return new SymmetricMatrix(result, Order, DefiniteProperty.Unknown);
        }

        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is SymmetricMatrix casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else if (other is ISymmetricMatrix otherSYM)
            {
                Preconditions.CheckSameMatrixDimensions(this, other);
                for (int j = 0; j < NumColumns; ++j)
                {
                    for (int i = 0; i <= j; ++i)
                    {
                        int index1D = Find1DIndex(i, j);
                        this.data[index1D] = binaryOperation(this.data[index1D], other[i, j]);
                    }
                }
            }
            else
            {
                throw new SymmetricPatternModifiedException("This operation is legal only if the other matrix is also symmetric.");
            }
        }

        public void DoEntrywiseIntoThis(SymmetricMatrix other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            for (int i = 0; i < data.Length; ++i) this.data[i] = binaryOperation(this.data[i], other.data[i]);
            Definiteness = DefiniteProperty.Unknown;
        }

        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public SymmetricMatrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            var result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i)
            {
                result[i] = unaryOperation(data[i]);
            }
            return new SymmetricMatrix(result, NumRows, DefiniteProperty.Unknown);
        }

        void IMatrix.DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            DoToAllEntriesIntoThis(unaryOperation);
        }

        void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < NumRows * NumColumns; ++i)
            {
                data[i] = unaryOperation(data[i]);
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        /// <summary>
        /// Calculates some factorization of the symmetric matrix.
        /// </summary>
        /// <param name="tryCholesky">True to apply the Cholesky factorization, before resorting to more expensive procedures.
        ///     Since Cholesky has a complecity of O(1/3*n^3), this will further increase the computational cost for matrices 
        ///     that are not positive definite. Therefore set it to false if you are sure that the matrix is not positive
        ///     definite.</param>
        /// <returns></returns>
        public ITriangulation Factorize(bool tryCholesky = true)
        {
            // This should throw an exception if the posdef assumption is wrong.
            if (Definiteness == DefiniteProperty.PositiveDefinite) return FactorCholesky();  

            if (tryCholesky)
            {
                try
                {
                    return FactorCholesky();
                }
                catch (IndefiniteMatrixException) { }
            }

            return FactorBunchKaufman();
        }

        public BunchKaufmanFactorization FactorBunchKaufman()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates the Cholesky factorization of the matrix. Will throw <see cref="IndefiniteMatrixException"/> if the 
        /// matrix is not positive definite.
        /// </summary>
        /// <returns></returns>
        public CholeskyPacked FactorCholesky()
        {
            var factor = CholeskyPacked.Factorize(Order, data);
            Definiteness = DefiniteProperty.PositiveDefinite; // An exception would have been thrown otherwise.
            return factor;
        }

        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SymmetricMatrix casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else return DoEntrywise(otherMatrix, (x1, x2) => thisCoefficient * x1 + otherCoefficient * x2); //TODO: optimize this
        }

        public SymmetricMatrix LinearCombination(double thisCoefficient, SymmetricMatrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref result[0], 1);
            return new SymmetricMatrix(result, NumColumns, DefiniteProperty.Unknown);
        }

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SymmetricMatrix casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else if (otherMatrix is ISymmetricMatrix otherSYM)
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                for (int j = 0; j < NumColumns; ++j)
                {
                    for (int i = 0; i <= j; ++i)
                    {
                        int index1D = Find1DIndex(i, j);
                        this.data[index1D] = thisCoefficient * this.data[index1D] + otherCoefficient * otherMatrix[i, j];
                    }
                }
            }
            else
            {
                throw new SymmetricPatternModifiedException("This operation is legal only if the other matrix is also symmetric.");
            }
        }

        public void LinearCombinationIntoThis(double thisCoefficient, SymmetricMatrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref this.data[0], 1);
            this.Definiteness = DefiniteProperty.Unknown;
        }

        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(other, this, transposeOther, transposeThis);
        }

        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(this, other, transposeThis, transposeOther);
        }

        public Vector MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            if (vector is Vector) return MultiplyRight((Vector)vector, transposeThis);
            else throw new NotImplementedException();
        }

        /// <summary>
        /// Matrix vector multiplication, with the vector on the right: matrix * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
        /// <returns></returns>
        public Vector MultiplyRight(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(this.NumColumns, vector.Length);
            double[] result = new double[NumRows];
            CBlas.Dspmv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, Order,
                1.0, ref data[0], ref vector.InternalData[0], 1, 0.0, ref result[0], 1);
            return Vector.CreateFromArray(result, false);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            for (int j = 0; j < data.Length; ++j)
            {
                aggregator = processEntry(data[Find1DIndex(j, j)], aggregator);
                for (int i = 0; i < j; ++j)
                {
                    // A[j,i] = A[i,j], but doubling it will not work for all reductions
                    int idx1D = Find1DIndex(i, j);
                    aggregator = processEntry(data[idx1D], aggregator);
                    aggregator = processEntry(data[idx1D], aggregator); 
                }
            }
            // no zeros implied
            return finalize(aggregator);
        }

        IMatrixView IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public SymmetricMatrix Scale(double scalar)
        {
            int numStoredEntries = this.data.Length;
            double[] result = new double[numStoredEntries];
            Array.Copy(this.data, result, numStoredEntries);
            CBlas.Dscal(numStoredEntries, scalar, ref result[0], 1);
            return new SymmetricMatrix(result, this.Order, this.Definiteness);
        }

        public void ScaleIntoThis(double scalar) => CBlas.Dscal(data.Length, scalar, ref data[0], 1);

        // Not very efficient
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            Definiteness = DefiniteProperty.Unknown;
            if ((rowIdx < 0) || (rowIdx >= Order) || (colIdx < 0) || (colIdx >= Order))
            {
                throw new IndexOutOfRangeException($"Invalid indices: ({rowIdx}, {colIdx})");
            }
            if (rowIdx <= colIdx) data[Find1DIndex(rowIdx, colIdx)] = value;
            else data[Find1DIndex(colIdx, rowIdx)] = value;
        }

        public IMatrixView Transpose()
        {
            return Transpose(true);
        }

        public SymmetricMatrix Transpose(bool copyInternalArray)
        {
            return SymmetricMatrix.CreateFromArray(data, Order, Definiteness, copyInternalArray);
        }


        public void WriteToFile(string path, bool append = false)
        {
            using (var writer = new System.IO.StreamWriter(path, append))
            {
#if DEBUG
                writer.AutoFlush = true; // To look at intermediate output at certain breakpoints
#endif
                for (int i = 0; i < data.Length; ++i)
                    writer.WriteLine(data[i].ToString("g17", new System.Globalization.CultureInfo("en-US", false).NumberFormat));
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int Find1DIndex(int i, int j)
        {
            return i + (j * (j + 1)) / 2;
        }
    }
}
