using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;

//TODO: align data using mkl_malloc, at least when creating placeholder copies of the matrix or when it is initialized from 0.
namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    /// <summary>
    /// General square matrix. Dense (full) storage. Uses MKL. Stored as 1D column major array.
    /// </summary>
    public class SquareMatrixMKL: MatrixMKL
    {
        private SquareMatrixMKL(double[] data, int order): base(data, order, order)
        {
            this.Order = order;
        }

        /// <summary>
        /// The number of rows or columns of the matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Create a new <see cref="SquareMatrixMKL"/> from a provided array. The array will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional containing the elements of the matrix. 
        /// Its lengths in both dimensions must be the same.</param>
        /// <returns></returns>
        public static new SquareMatrixMKL CreateFromArray(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1);
            if (numRows != numCols)
            {
                string msg = string.Format("Provided array must have the same dimensions, but was ({0}x{1})", numRows, numCols);
                throw new NonMatchingDimensionsException(msg);
            }
            return new SquareMatrixMKL(Conversions.Array2DToFullColMajor(array2D), numRows);
        }

        /// <summary>
        /// Create a new <see cref="SquareMatrixMKL"/> with the specified dimensions and all entries equal to 0.
        /// </summary>
        /// <param name="order">The number of rows or columns of the matrix.</param>
        /// <returns></returns>
        public static SquareMatrixMKL CreateZero(int order)
        {
            double[] data = new double[order * order];
            return new SquareMatrixMKL(data, order);
        }

        public double CalcDeterminant()
        {
            try
            {
                LUFactorizationMKL factorization = FactorLU();
                return factorization.CalcDeterminant();
            }
            catch (SingularMatrixException)
            {
                return 0.0;
            }
        }

        public double CalcDeterminantInPlace()
        {
            throw new NotImplementedException();
        }

        public LUFactorizationMKL FactorLU()
        {
            // Copy matrix. This may exceed available memory and needs an extra O(n^2) accesses. 
            // To avoid these, use the ~InPlace version.
            double[] lowerUpper = new double[data.Length];
            Array.Copy(data, lowerUpper, data.Length);

            int n = Order;
            int[] permutation = new int[n];
            int info = MKLUtilities.DefaultInfo;
            Lapack.Dgetrf(ref n, ref n, ref lowerUpper[0], ref n, ref permutation[0], ref info);

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
                string msg = string.Format("The {0}th parameter has an illegal value.", -info)
                    + " Please contact the developer responsible for the linear algebra project.";
                throw new MKLException(msg);
            }
            else if (info > 0)
            {
                int idx = info - 1;
                string msg = "The factorization has been completed, but U is singular."
                    + string.Format(" The first zero pivot is U[{0}, {1}] = 0.", idx, idx);
                throw new SingularMatrixException(msg);
            }

            return new LUFactorizationMKL(n, lowerUpper, permutation);
        }

        public LUFactorizationMKL FactorLUInPlace()
        {
            throw new NotImplementedException();
        }

        public SquareMatrixMKL Invert()
        {
            throw new NotImplementedException();
        }

        public SquareMatrixMKL InvertInPlace()
        {
            throw new NotImplementedException();
        }
    }
}
