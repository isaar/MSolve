using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations
{
    public class CholeskySkyline : IIndexable2D, ISparseMatrix
    {
        public const double PivotTolerance = 1e-15;

        private readonly double[] values;
        private readonly int[] diagOffsets;

        private CholeskySkyline(int order, double[] values, int[] diagOffsets)
        {
            this.NumColumns = order;
            this.values = values;
            this.diagOffsets = diagOffsets;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return NumColumns; } }

        // TODO: should I add index bound checking?
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                SkylineMatrix.ProcessIndices(ref rowIdx, ref colIdx);
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffsets[colIdx] - 1; // excluding diagonal
                if (entryHeight > maxColumnHeight) return 0.0; // outside stored non zero pattern
                else return values[diagOffsets[colIdx] + entryHeight];
            }
        }

        /// <summary>
        /// Cholesky factorization algorithm in skyline format. Does not need any extra memory. The non zero values array will be 
        /// overwritten, while the indexing array will be used as is. Throws a <see cref="IndefiniteMatrixException"/> if the 
        /// unfactorized skyline matrix is not positive definite.
        /// </summary>
        /// <param name="order">The number of rows = number of columns of the matrix</param>
        /// <param name="skyValues">The non zero entries of the original <see cref="SkylineMatrix"/>. They will  be
        ///     overwritten.</param>
        /// <param name="skyDiagOffsets">The indexes of the diagonal entries into skyValues. They will be stored as is inside
        ///     the new <see cref="CholeskySkyline"/> object. However they do not need copying, since they will not be altered
        ///     during or after the factorization.</param>
        /// <param name="tolerance">If a diagonal entry is closer to zero than this tolerance, an 
        ///     <see cref="IndefiniteMatrixException"/> exception will be thrown.</param>
        /// <returns></returns>
        public static CholeskySkyline Factorize(int order, double[] skyValues, int[] skyDiagOffsets, 
            double tolerance = CholeskySkyline.PivotTolerance)
        {
            // Copied from Stavroulakis code.
            var zemCols = new List<int>();
            var zems = new List<double[]>();

            int kFix = 0;
            for (int n = 0; n < order; n++)
            {
                int KN = skyDiagOffsets[n];
                int KL = KN + 1;
                int KU = skyDiagOffsets[n + 1] - 1;
                int KH = KU - KL;
                if (KH < 0) continue;

                int K;
                if (KH > 0)
                {
                    K = n - KH;
                    //int IC = 0;
                    int KLT = KU;
                    for (int j = 0; j < KH; j++)
                    {
                        //IC++;
                        KLT--;
                        int KI = skyDiagOffsets[K];
                        int ND = skyDiagOffsets[K + 1] - KI - 1;
                        if (ND > 0)
                        {
                            //int KK = Math.Min(IC, ND);
                            int KK = Math.Min(j + 1, ND);
                            double C = 0;
                            for (int l = 1; l <= KK; l++)
                                C += skyValues[KI + l] * skyValues[KLT + l];
                            skyValues[KLT] -= C;
                        }
                        K++;
                    }
                }
                K = n;
                double B = 0;
                for (int KK = KL; KK <= KU; KK++)
                {
                    K--;
                    int KI = skyDiagOffsets[K];
                    //if (d[KI] == 0) throw new InvalidOperationException(String.Format("Zero element in diagonal at index {0}.", KI));
                    if (Math.Abs(skyValues[KI]) < tolerance)
                    {
                        throw new IndefiniteMatrixException($"Near-zero element in diagonal at index {KI}.");
                    }
                    double C = skyValues[KK] / skyValues[KI];
                    B += C * skyValues[KK];
                    skyValues[KK] = C;
                }
                skyValues[KN] -= B;

                if (Math.Abs(skyValues[KN]) < tolerance)
                {
                    skyValues[KN] = 1;
                    zemCols.Add(n);
                    int j1 = n;
                    zems.Add(new double[order]);
                    zems[kFix][j1] = 1;
                    for (int i1 = KN + 1; i1 <= KU; i1++)
                    {
                        j1--;
                        zems[kFix][j1] = -skyValues[i1];
                        skyValues[i1] = 0;
                    }
                    for (int irest = n + 1; irest < order; irest++)
                    {
                        int m1 = skyDiagOffsets[irest] + irest - n;
                        if (m1 <= skyDiagOffsets[irest + 1]) skyValues[m1] = 0;
                    }
                    kFix++;
                }
            }
            //isFactorized = true; // not needed here
            if (order < 2) return new CholeskySkyline(order, skyValues, skyDiagOffsets);

            for (int ifl = 0; ifl < kFix; ifl++)
            {
                int n = order - 1;
                for (int l = 1; l < order; l++)
                {
                    int KL = skyDiagOffsets[n] + 1;
                    int KU = skyDiagOffsets[n + 1] - 1;
                    if (KU - KL >= 0)
                    {
                        int k = n;
                        for (int KK = KL; KK <= KU; KK++)
                        {
                            k--;
                            zems[ifl][k] -= skyValues[KK] * zems[ifl][n];
                        }
                    }
                    n--;
                }
            }
            if (zemCols.Count > 0)
            {
                throw new IndefiniteMatrixException("Cholesky factorization can only applied to positive definite matrices.");
            }
            return new CholeskySkyline(order, skyValues, skyDiagOffsets);
        }

        public int CountNonZeros()
        {
            return values.Length;
        }

        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            return SkylineMatrix.CreateFromArrays(NumColumns, values, diagOffsets, false,  false).EnumerateNonZeros();
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return SkylineMatrix.CreateFromArrays(NumColumns, values, diagOffsets, false, false).Equals(other, tolerance);
        }

        /// <summary>
        /// A = transpose(U) * U
        /// </summary>
        /// <returns></returns>
        public TriangularUpper GetFactorU()
        {
            // The factorization A = transpose(u) * D * u, u = unit upper triangular is stored. Thus U = sqrt(D) * u.
            // Since D is diagonal, we need to scale each column j of u by sqrt(D[j,j]).
            var upper = TriangularUpper.CreateZero(NumColumns); 
            for (int j = 0; j < NumColumns; ++j)
            {
                int diagOffset = diagOffsets[j];
                double sqrtD = Math.Sqrt(values[diagOffset]); // The diagonal of D is stored instead of the unit diagonal entries of u
                int colTop = j - diagOffsets[j + 1] + diagOffset + 1;
                upper[j, j] = sqrtD; // U[j,j] = u[j,j] * sqrt(D[j,j]) = 1 * sqrt(D[j,j])
                for (int i = j - 1; i >= colTop; --i)
                {
                    int offset = diagOffset + j - i;
                    upper[i, j] = sqrtD * values[offset]; // U[:,j] = sqrt(D[j,j]) * u[:,j]
                }
            }
            int[] diagOffsetsCopy = new int[diagOffsets.Length];
            Array.Copy(diagOffsets, diagOffsetsCopy, diagOffsets.Length);
            return upper;
        }

        /// <summary>
        /// A = transpose(U) * D * U, U being unit upper triangular.
        /// </summary>
        /// <returns></returns>
        public (VectorMKL diagonal, TriangularUpper upper) GetFactorsDU()
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
            return (VectorMKL.CreateFromArray(diag), upper); 
        }

        public SparseFormat GetSparseFormat()
        {
            return SkylineMatrix.CreateFromArrays(NumColumns, values, diagOffsets, false, false).GetSparseFormat();
        }

        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            Preconditions.CheckSystemSolutionDimensions(this, rhs);

            //var e = DateTime.Now;
            //double[] result = new double[K.Rows];
            double[] result = rhs.CopyToArray();

            // RHS vector reduction
            int n;
            for (n = 0; n < NumColumns; n++)
            {
                int KL = diagOffsets[n] + 1;
                int KU = diagOffsets[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += values[KK] * result[k];
                    }
                    result[n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < NumColumns; n++) result[n] /= values[diagOffsets[n]];

            n = NumColumns - 1;
            for (int l = 1; l < NumColumns; l++)
            {
                int KL = diagOffsets[n] + 1;
                int KU = diagOffsets[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        result[k] -= values[KK] * result[n];
                    }
                }
                n--;
            }
            //var x = new List<TimeSpan>();
            //x.Add(DateTime.Now - e);

            return VectorMKL.CreateFromArray(result, false);
        }
    }
}
