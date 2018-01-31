using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public class LDUFactorization
    {
        private readonly int order;

        /// <summary>
        /// Stores the lower triangle L, diagonal D and upper triangle U in the same 2D array. 
        /// The unitary diagonal of L and U are implied.
        /// </summary>
        private readonly double[,] data;

        /// <summary>
        /// Vector format of the permutation matrix P*A=L*U. Use: permutation[originalRowIndex] = permutedRowIndex
        /// </summary>
        private readonly int[] permutation;

        public LDUFactorization(double[,] factorizedMatrix, int[] permutation)
        {
            this.data = factorizedMatrix;
            this.order = factorizedMatrix.GetLength(0);
        }

        /// <summary>
        /// Returns the determinant of the original matrix. det(A) = det(L*D*U) = det(L)*det(D)*det(U). 
        /// Since all these are triangular matrices their determinants is the product of their diagonal entries:
        /// det(L) = 1*1*...*1 = 1 = det(U). Thus det(A) = det(D) = D1*D2*...*Dn
        /// </summary>
        /// <returns></returns>
        public double CalcDeterminant()
        {
            double product = 1.0;
            for (int i = 0; i < order; ++i)
            {
                product *= data[i, i];
            }
            return product;
        }

        /// <summary>
        /// Solves A*x = b with respect to x, 
        /// using forward substitution, multiplication with the inverse diagonal and backward substitution. 
        /// </summary>
        /// <param name="b">The right hand side of the original system 
        /// (before any factorizations or permutations)</param>
        /// <returns></returns>
        public IVector Solve(IVectorView b)
        {
            double[] pb = PermuteRhs(b); //A*x=b => P*A*x=P*b
            double[] y = ForwardSubstitution(pb); //P*A*x=P*b => L*D*U*x=P*b => L*y=P*b => y=inv(L)*P*b 
            DivideWithDiagonal(y);
            return BackSubstitution(y); // U*x=y => x=inv*U)*y
        }

        private double[] PermuteRhs(IVectorView b)
        {
            double[] pb = new double[b.Length];
            for (int i = 0; i < pb.Length; ++i)
            {
                pb[i] = b[permutation[i]];
            }
            return pb;
        }

        private double[] ForwardSubstitution(double[] pb)
        {
            double[] y = new double[order];
            for (int row = 0; row < order; ++row)
            {
                double nominator = pb[row];
                for (int col = 0; col < row; ++col)
                {
                    nominator -= data[row, col] * pb[col]; //row>col thus data[row,col] accesses lower triangle entries
                }
                y[row] = nominator; // denominator = diagonal entry of L = 1
            }
            return y;
        }

        private void DivideWithDiagonal(double[] y)
        {
            for (int i = 0; i < order; ++i)
            {
                y[i] /= data[i, i];
            }
        }

        private IVector BackSubstitution(double[] y)
        {
            double[] x = new double[order];
            for (int row = order - 1; row >= 0; --row)
            {
                double nominator = y[row];
                for (int col = row + 1; col < order; ++col)
                {
                    nominator -= data[row, col] * y[col]; //col>row thus data[row,col] accesses upper triangle entries
                }
                x[row] = nominator; // denominator = diagonal entry of U = 1
            }
            return DenseVector.CreateFromArray(x, false);
        }
    }
}
