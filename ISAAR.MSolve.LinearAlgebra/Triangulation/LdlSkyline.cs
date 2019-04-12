using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//also implement ISymmetricSparseMatrix.
//TODO: the heavy operations should be ported to a C dll. and called from there..
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// LDL factorization of a symmetric matrix: A = L * D * transpose(L) = transpose(U) * D * U. The matrix is stored in 
    /// skyline format. Only the active columns of the upper triangle part of the matrix is stored and factorized. 
    /// LDL factorization is unique and stable for symmetric positive definite matrices. It may also succeed even if the matrix 
    /// is indefinite, but generally pivoting is required, which is not implemented by <see cref="LdlSkyline"/>. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class LdlSkyline : SkylineFactorizationBase
    {
        /// <summary>
        /// The default value under which a diagonal entry (pivot) is considered to be 0 during LDL factorization.
        /// </summary>
        public const double PivotTolerance = 1e-15;

        private LdlSkyline(int order, double[] values, int[] diagOffsets) : base(order, values, diagOffsets)
        { }

        /// <summary>
        /// Calculates the LDL factorization of a symmetric matrix, such that A = transpose(U) * D * U. 
        /// Does not need any extra memory.
        /// </summary>
        /// <param name="order">The number of rows/ columns of the square matrix.</param>
        /// <param name="skyValues">The non-zero entries of the original <see cref="SkylineMatrix"/>. Ths array will be
        ///     overwritten during the factorization.</param>
        /// <param name="skyDiagOffsets">The indexes of the diagonal entries into <paramref name="skyValues"/>. The new 
        ///     <see cref="LdlSkyline"/> instance will hold a reference to <paramref name="skyDiagOffsets"/>. However they 
        ///     do not need copying, since they will not be altered during or after the factorization.</param>
        /// <param name="pivotTolerance">If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the 
        ///     original matrix is not invertible and an <see cref="SingularMatrixException"/> will be thrown.</param>
        ///<exception cref="SingularMatrixException">Thrown if the original skyline matrix turns out to be singular.</exception>
        public static LdlSkyline Factorize(int order, double[] skyValues, int[] skyDiagOffsets, 
            double pivotTolerance = LdlSkyline.PivotTolerance)
        {
            (List<int> dependentColumns, List<double[]> nullSpaceBasis) = 
                FactorizeInternal(order, skyValues, skyDiagOffsets, pivotTolerance);

            if (dependentColumns.Count > 0) throw new SingularMatrixException("The matrix is singular.");
            return new LdlSkyline(order, skyValues, skyDiagOffsets);
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        public override double CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Copies the upper triangular matrix U and diagonal matrix D that resulted from the Cholesky factorization: 
        /// A = transpose(U) * D * U to a new <see cref="TriangularUpper"/> matrix.
        /// </summary>
        public (Vector diagonal, TriangularUpper upper) GetFactorsDU() //TODO: not sure if this is ever needed.
        {
            double[] diag = new double[NumColumns];
            var upper = TriangularUpper.CreateZero(NumColumns);
            for (int j = 0; j < NumColumns; ++j)
            {
                int diagOffset = diagOffsets[j];
                int colTop = j - diagOffsets[j + 1] + diagOffset + 1;
                diag[j] = values[diagOffset];
                upper[j, j] = 1.0;
                for (int i = j - 1; i >= colTop; --i)
                {
                    int offset = diagOffset + j - i;
                    upper[i, j] = values[offset];
                }

            }
            return (Vector.CreateFromArray(diag, false), upper);
        }

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
        /// </summary>
        public override void SolveLinearSystem(Vector rhs, Vector solution)
        {
            Preconditions.CheckSystemSolutionDimensions(this, rhs);
            Preconditions.CheckMultiplicationDimensions(NumColumns, solution.Length);
            Solve(NumColumns, values, diagOffsets, rhs.RawData, solution.RawData);
        }

        /// <summary>
        /// Solves the linear systems A * X = B, where A is the original matrix (before the factorization), 
        /// B = <paramref name="rhsVectors"/> and X is the matrix containing the solution vectors, which will overwrite the 
        /// provided <paramref name="solutionVectors"/>.
        /// </summary>
        /// <param name="rhsVectors">
        /// A matrix that contains the right hand side vectors as its columns. Constraints: 
        /// a) Its <see cref="IIndexable2D.NumRows"/> must be equal to the <see cref="IIndexable2D.NumRows"/> of the original 
        /// matrix A. b) Its <see cref="IIndexable2D.NumColumns"/> must be equal to the <see cref="IIndexable2D.NumColumns"/> of
        /// <paramref name="solutionVectors"/>.
        /// </param>
        /// <param name="solutionVectors">
        /// Output matrix that will be overwritten with the solutions of the linear system as its columns. Constraints:
        /// a) Its <see cref="IIndexable2D.NumRows"/> must be equal to the <see cref="IIndexable2D.NumRows"/> of the original 
        /// matrix A. b) Its <see cref="IIndexable2D.NumColumns"/> must be equal to the <see cref="IIndexable2D.NumColumns"/> of
        /// <paramref name="rhsVectors"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="rhsVectors"/> or <paramref name="solutionVectors"/> violate the described constraints.
        /// </exception>
        public void SolveLinearSystems(Matrix rhsVectors, Matrix solutionVectors)
        {
            Preconditions.CheckSystemSolutionDimensions(this.NumRows, rhsVectors.NumRows);
            Preconditions.CheckMultiplicationDimensions(this.NumColumns, solutionVectors.NumRows);
            Preconditions.CheckSameColDimension(rhsVectors, solutionVectors);

            for (int j = 0; j < rhsVectors.NumColumns; ++j)
            {
                int offset = j * NumRows;
                SolveWithOffsets(NumColumns, values, diagOffsets, rhsVectors.RawData, offset, solutionVectors.RawData, offset);
            }
        }

        internal static (List<int> dependentColumns, List<double[]> nullSpaceBasis) FactorizeInternal(int order, double[] values,
            int[] diagOffsets, double pivotTolerance)
        {
            // Copied from Stavroulakis code.
            var zemCols = new List<int>();
            var zems = new List<double[]>();

            int kFix = 0;

            #region factorization of positive definite matrix
            for (int n = 0; n < order; n++)
            {
                int KN = diagOffsets[n];
                int KL = KN + 1;
                int KU = diagOffsets[n + 1] - 1;
                int KH = KU - KL; // height of current column
                if (KH < 0) continue; // TODO: Shouldn't this throw an exception? The height of Skyline columns should be at least 0

                int K;
                if (KH > 0)
                {
                    K = n - KH; // index of top entry in column n
                    //int IC = 0;
                    int KLT = KU; // Offset of non-zero entries in column n. Starts from the top (KU-1) and goes down till the diagonal.

                    // Dot product of parts of columns n, j: sum(U[k,n]*U[k,j]). The formula is usually: 
                    // sum(L[i,n]*L[k,j]/D[k,k]), but the diagonal entries will be processed in the next nested loop.
                    for (int j = 0; j < KH; j++)
                    {
                        //IC++;
                        KLT--;
                        int KI = diagOffsets[K];
                        int ND = diagOffsets[K + 1] - KI - 1;
                        if (ND > 0)
                        {
                            //int KK = Math.Min(IC, ND);
                            int KK = Math.Min(j + 1, ND);
                            double C = 0;
                            for (int l = 1; l <= KK; l++)
                                C += values[KI + l] * values[KLT + l];
                            values[KLT] -= C;
                        }
                        K++;
                    }
                }
                K = n;
                double B = 0;
                for (int KK = KL; KK <= KU; KK++)
                {
                    K--;
                    int KI = diagOffsets[K]; // offset of diagonal entry of previous column: D[k,k] in the original algorithm

                    //TODO: values[diagOffsets[K]] (K<n) is the diagonal entry of a previous column and was updated in a 
                    //      previous iteration. When that happened, if it was near-zero, it was set to 1 and column K was added 
                    //      to the nullspace. Therefore this check is useless. Moreover it doesn't make sense, since this 
                    //      semidefinite-LDL algorithm is supposed to work for singular matrices too.
                    if (Math.Abs(values[KI]) < pivotTolerance)
                    {
                        throw new SingularMatrixException($"Near-zero element in diagonal at index {KI}."); //TODO: Not sure if this happens only to singular matrices.
                    }
                    double C = values[KK] / values[KI];
                    B += C * values[KK];
                    values[KK] = C;
                }
                values[KN] -= B;
                #endregion

                #region partial extraction of zero energy modes (aka rigid body motions / nullspace) during factorization
                //double pivot = Math.Abs(values[KN]); //TODO: for debugging
                if (Math.Abs(values[KN]) < pivotTolerance)
                {
                    values[KN] = 1;
                    zemCols.Add(n);
                    int j1 = n;
                    zems.Add(new double[order]);
                    zems[kFix][j1] = 1; // the diagonal becomes 1
                    for (int i1 = KN + 1; i1 <= KU; i1++) // The other entries of column n become 0
                    {
                        j1--;
                        zems[kFix][j1] = -values[i1]; // partial calculation of the zero energy mode
                        values[i1] = 0;
                    }
                    for (int irest = n + 1; irest < order; irest++) // The other entries of row n become 0
                    {
                        int m1 = diagOffsets[irest] + irest - n;
                        if (m1 <= diagOffsets[irest + 1]) values[m1] = 0; //TODO: Why <= instead of < ?
                    }
                    kFix++;
                }
                #endregion
            }
            //isFactorized = true; // not needed here
            if (order < 2) return (zemCols, zems);

            #region back substitution to finish the calculation of the zero energy modes
            for (int ifl = 0; ifl < kFix; ifl++)
            {
                int n = order - 1;
                for (int l = 1; l < order; l++)
                {
                    int KL = diagOffsets[n] + 1;
                    int KU = diagOffsets[n + 1] - 1;
                    if (KU - KL >= 0)
                    {
                        int k = n;
                        for (int KK = KL; KK <= KU; KK++)
                        {
                            k--;
                            zems[ifl][k] -= values[KK] * zems[ifl][n];
                        }
                    }
                    n--;
                }
            }
            #endregion

            return (zemCols, zems);
        }

        internal static void Solve(int order, double[] valuesA, int[] diagOffsetsA, double[] vectorB, double[] vectorX)
        {
            // Copied from Stavroulakis code.

            // Copy the b vector
            Array.Copy(vectorB, vectorX, order);

            // RHS vector reduction
            int n;
            for (n = 0; n < order; n++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL) //TODO: Is it possible that this doesn't hold? Then something is wrong with the skyline format
                {
                    int k = n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += valuesA[KK] * vectorX[k];
                    }
                    vectorX[n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) vectorX[n] /= valuesA[diagOffsetsA[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL) //TODO: Is it possible that this doesn't hold? Then something is wrong with the skyline format
                {
                    int k = n;
                    double xn = vectorX[n];
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        vectorX[k] -= valuesA[KK] * xn;
                    }
                }
                n--;
            }
        }

        private static void SolveWithOffsets(int order, double[] valuesA, int[] diagOffsetsA,
            double[] vectorB, int offsetB, double[] vectorX, int offsetX)
        {
            // Copied from Stavroulakis code.

            // Copy the y vector
            Array.Copy(vectorB, offsetB, vectorX, offsetX, order);

            // RHS vector reduction
            int n;
            for (n = 0; n < order; n++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL) //TODO: can't I avoid this check, by accessing the entries in a smarter way?
                {
                    int k = offsetX + n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += valuesA[KK] * vectorX[k];
                    }
                    vectorX[offsetX + n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) vectorX[offsetX + n] /= valuesA[diagOffsetsA[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = diagOffsetsA[n] + 1;
                int KU = diagOffsetsA[n + 1] - 1;
                if (KU >= KL) //TODO: can't I avoid this check, by accessing the entries in a smarter way?
                {
                    int k = offsetX + n;
                    double xn = vectorX[k];
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        vectorX[k] -= valuesA[KK] * xn;
                    }
                }
                n--;
            }
        }
    }
}
