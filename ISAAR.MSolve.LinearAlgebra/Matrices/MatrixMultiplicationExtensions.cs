using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Most of the multiplications C = A^T * B* A will need updating if the individual matrix multiplication are updated, in 
//      order to ensure that they are optimal.
//TODO: Variations of the multiplications C = A^T * B* A, where A is IMatrixView.

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public static class MatrixMultiplicationExtensions
    {
        /// <summary>
        /// Performs the matrix-vector multiplication: oper(<paramref name="matrix"/>) * <paramref name="vector"/>.
        /// To multiply <paramref name="matrix"/> * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * <paramref name="matrix"/>, set <paramref name="transposeThis"/> to true.
        /// </summary>
        /// <param name="matrix">The matrix to multiply.</param>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to the 
        ///     <see cref="IIndexable2D.NumColumns"/> of oper(<paramref name="matrix"/>).</param>
        /// <param name="transposeThis">If true, oper(<paramref name="matrix"/>) = transpose(<paramref name="matrix"/>). 
        ///     Otherwise oper(<paramref name="matrix"/>) = <paramref name="matrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of 
        ///     oper(<paramref name="matrix"/>).</exception>
        public static double[] MultiplyRight(this CscMatrix matrix, double[] vector, bool transposeThis)
        { //TODO: delete this once legacy vectors, matrices are no longer used.
            var asVector = Vector.CreateFromArray(vector, false);
            return matrix.Multiply(asVector, transposeThis).RawData;
        }

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="csr"/>) * <paramref name="other"/> * <paramref name="csr"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="csr">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTimesOtherTimesThisTranspose(this CsrMatrix csr, IMatrixView other)
        { //TODO: Perhaps this should not be a separate method than the one where other is Matrix
            if (other is Matrix dense) csr.ThisTimesOtherTimesThisTranspose(dense);

            // As of 15 December 2018, right multiplying is more efficient than left multiplying, since the resulting matrix 
            // of right multiplying is column major and we can leverage BLAS. Therefore we want the right multiplication to 
            // happen for the larger matrix: csr * matrix.
            if (other.NumColumns < csr.NumRows)
            {
                Matrix temp = csr.MultiplyLeft(other, true, false); // temp is larger than other
                return csr.MultiplyRight(temp, false, false);
            }
            else
            {
                Matrix temp = csr.MultiplyRight(other, false, false); // other is larger than temp
                return csr.MultiplyLeft(temp, true, false);
            }
        }

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="csr"/>) * <paramref name="other"/> * <paramref name="csr"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="csr">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTimesOtherTimesThisTranspose(this CsrMatrix csr, Matrix other)
        { 
            // As of 15 December 2018, right multiplying is more efficient than left multiplying, since the resulting matrix 
            // of right multiplying is column major and we can leverage BLAS. Therefore we want the right multiplication to 
            // happen for the larger matrix: csr * matrix.
            if (other.NumColumns < csr.NumRows)
            {
                Matrix temp = csr.MultiplyLeft(other, true, false); // temp is larger than other
                return csr.MultiplyRight(temp, false, false);
            }
            else
            {
                Matrix temp = csr.MultiplyRight(other, false, false); // other is larger than temp
                return csr.MultiplyLeft(temp, true, false);
            }
        }

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="csr"/>) * <paramref name="other"/> * <paramref name="csr"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="csr">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTimesOtherTimesThisTranspose(this CsrMatrix csr, SymmetricMatrix other)
            => throw new NotImplementedException("Placeholder for when SymmetricMatrix is fully implemented");

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="csc"/>) * <paramref name="other"/> * <paramref name="csc"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="csc">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTransposeTimesOtherTimesThis(this CscMatrix csc, IMatrixView other)
        { //TODO: Perhaps this should not be a separate method than the one where other is Matrix
            if (other is Matrix dense) return csc.ThisTransposeTimesOtherTimesThis(dense);

            // As of 15 December 2018, right multiplying is more efficient than left multiplying, since the resulting matrix 
            // of right multiplying is column major and we can leverage BLAS. Therefore we want the right multiplication to 
            // happen for the larger matrix: transpose(csc) * matrix.
            if (other.NumColumns < csc.NumColumns)
            {
                Matrix temp = csc.MultiplyLeft(other, false, false); // temp is larger than other
                return csc.MultiplyRight(temp, true, false);
            }
            else
            {
                Matrix temp = csc.MultiplyRight(other, true, false); // other is larger than temp
                return csc.MultiplyLeft(temp, false, false);
            }
        }

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="csc"/>) * <paramref name="other"/> * <paramref name="csc"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="csc">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTransposeTimesOtherTimesThis(this CscMatrix csc, Matrix other)
        { 
            // As of 15 December 2018, right multiplying is more efficient than left multiplying, since the resulting matrix 
            // of right multiplying is column major and we can leverage BLAS. Therefore we want the right multiplication to 
            // happen for the larger matrix: transpose(csc) * matrix.
            if (other.NumColumns < csc.NumColumns)
            {
                Matrix temp = csc.MultiplyLeft(other, false, false); // temp is larger than other
                return csc.MultiplyRight(temp, true, false);
            }
            else
            {
                Matrix temp = csc.MultiplyRight(other, true, false); // other is larger than temp
                return csc.MultiplyLeft(temp, false, false);
            }
        }

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="csc"/>) * <paramref name="other"/> * <paramref name="csc"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="csc">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTransposeTimesOtherTimesThis(this CscMatrix csc, SymmetricMatrix other)
            => throw new NotImplementedException("Placeholder for when SymmetricMatrix is fully implemented");

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="dense"/>) * <paramref name="other"/> * <paramref name="dense"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="dense">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTransposeTimesOtherTimesThis(this Matrix dense, IMatrixView other)
        { //TODO: Perhaps this should not be a separate method than the one where other is Matrix
            if (other is Matrix otherDense) return dense.ThisTransposeTimesOtherTimesThis(otherDense);

            //TODO: perhaps I should dedice about the order of multiplications depending on the number of columns of each matrix
            Matrix temp = other.MultiplyRight(dense, false, false);
            return dense.MultiplyRight(temp, true, false);
        }

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="dense"/>) * <paramref name="other"/> * <paramref name="dense"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="dense">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTransposeTimesOtherTimesThis(this Matrix dense, Matrix other)
        { 
            // The order of operations is not important if all matrices are in full storage format.
            Matrix temp = other.MultiplyRight(dense, false, false);
            return dense.MultiplyRight(temp, true, false);
        }

        /// <summary>
        /// Performs the operation: result = transpose(<paramref name="dense"/>) * <paramref name="other"/> * <paramref name="dense"/>
        /// in an efficient way, by appropriately selecting which methods should be called for these matrices and in what order. 
        /// </summary>
        /// <param name="dense">The matrix that will be multiplied "outside".</param>
        /// <param name="other">The matrix that will be multiplied "inside". It must be square.</param>
        public static Matrix ThisTransposeTimesOtherTimesThis(this Matrix dense, SymmetricMatrix other)
            => throw new NotImplementedException("Placeholder for when SymmetricMatrix is fully implemented");
    }
}
