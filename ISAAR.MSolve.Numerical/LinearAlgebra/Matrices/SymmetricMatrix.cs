using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

//TODO: align data using mkl_malloc
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// Symmetric matrix. Only the upper triangle is stored in Packed format (only stores the n*(n+1)/2 non zeros) and column 
    /// major order. Uses MKL.
    /// </summary>
    public class SymmetricMatrix: IMatrixView
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
        /// Sets A[<see cref="i"/>, <see cref="j"/>] = A[<see cref="j"/>, <see cref="i"/>] = <see cref="value"/>. Not efficient.
        /// </summary>
        /// <param name="i">The row index: 0 &lt;= i &lt; <see cref="Order"/></param>
        /// <param name="j">The column index: 0 &lt;= j &lt; <see cref="Order"/></param>
        public void SetSymmetric(int i, int j, double value) //TODO: Should I add an efficient version without error checking.
        {
            Definiteness = DefiniteProperty.Unknown;
            if ((i < 0) || (i >= Order) || (j < 0) || (j >= Order))
            {
                throw new IndexOutOfRangeException($"Invalid indices: ({i}, {j})");
            }
            if (i <= j) data[Find1DIndex(i, j)] = value;
            else data[Find1DIndex(j, i)] = value;
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

        /// <summary>
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="Order"/> 
        /// and length(1) = <see cref="Order"/>. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="Order"/>, <see cref="Order"/>] array 
        /// with the entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            return Conversions.PackedUpperColMajorToArray2D(data, Order);
        }

        public Matrix CopyToGeneralMatrix()
        {
            double[] fullData = Conversions.PackedUpperColMajorToFullSymmColMajor(data, Order);
            return Matrix.CreateFromArray(fullData, Order, Order, false);
        }

        public IEntrywiseOperable DoEntrywise(IEntrywiseOperable other, Func<double, double, double> binaryOperation)
        {
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        public IEntrywiseOperable DoToAllEntries(Func<double, double> unaryOperation)
        {
            var result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i)
            {
                result[i] = unaryOperation(data[i]);
            }
            return new SymmetricMatrix(result, NumRows, DefiniteProperty.Unknown);
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        /// <summary>
        /// Matrix vector multiplication, with the vector on the right: matrix * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
        /// <returns></returns>
        public VectorMKL MultiplyRight(VectorMKL vector)
        {
            Preconditions.CheckMultiplicationDimensions(this.NumColumns, vector.Length);
            double[] result = new double[NumRows];
            CBlas.Dspmv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, Order,
                1.0, ref data[0], ref vector.InternalData[0], 1, 0.0, ref result[0], 1);
            return VectorMKL.CreateFromArray(result, false);
        }

        /// <summary>
        /// Calculates some factorization of the symmetric matrix.
        /// </summary>
        /// <param name="tryCholesky">True to apply the Cholesky factorization, before resorting to more expensive procedures.
        ///     Since Cholesky has a complecity of O(1/3*n^3), this will further increase the computational cost for matrices 
        ///     that are not positive definite. Therefore set it to false if you are sure that the matrix is not positive
        ///     definite.</param>
        /// <returns></returns>
        public IFactorization Factorize(bool tryCholesky = true)
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
        public CholeskyFactorization FactorCholesky()
        {
            // Copy matrix. This may exceed available memory and needs an extra O(n^2) accesses. 
            // To avoid these, use the ~InPlace version.
            double[] upper = new double[data.Length];
            Array.Copy(data, upper, data.Length);

            // Call MKL
            int n = Order;
            int[] permutation = new int[n];
            int info = MKLUtilities.DefaultInfo;
            Lapack.Dpptrf("U", ref n, ref upper[0], ref info);

            // Check MKL execution
            if (info == MKLUtilities.DefaultInfo)
            {
                // first check the default info value, since it lies in the other intervals.
                // info == dafeult => the MKL call did not succeed. 
                // info > 0 should not be returned at all by MKL, but it is here for completion.
                throw new MKLException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else if (info < 0)
            {
                string msg = $"The {-info}th parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.";
                throw new MKLException(msg);
            }
            else if (info > 0)
            {
                string msg = "The leading minor of order " + (info -1) + " (and therefore the matrix itself) is not"
                + " positive-definite, and the factorization could not be completed.";
                throw new IndefiniteMatrixException(msg);
            }

            Definiteness = DefiniteProperty.PositiveDefinite;
            return new CholeskyFactorization(upper, n);
        }

        public void WriteToConsole()
        {
            DenseStrategies.WriteToConsole(this);
        }

        public void WriteToFile(string path, bool append = false)
        {
            DenseStrategies.WriteToFile(this, path, append);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int Find1DIndex(int i, int j)
        {
            return i + (j * (j + 1)) / 2;
        }

        public IMatrixView MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            throw new NotImplementedException();
        }

        public IMatrixView MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            throw new NotImplementedException();
        }

        public IVectorView MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            throw new NotImplementedException();
        }

        public IMatrixView Transpose()
        {
            throw new NotImplementedException();
        }
    }
}
