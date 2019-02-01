using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// LDL-like factorization of a symmetric positive semi-definite matrix: A = L * D * transpose(L) = transpose(U) * D * U.
    /// This factorization will apply the LDL algorithm on the independent columns of the matrix, set the dependent ones 
    /// equal to columns of the identity matrix and return the nullspace of the matrix.
    /// The matrix is stored in skyline format. Only the active columns of the upper triangle part of the matrix is stored and 
    /// factorized. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SemidefiniteLdlSkyline
    {
        /// <summary>
        /// The default value under which a diagonal entry (pivot) is considered to be 0 during Cholesky factorization.
        /// For the default value see "The FETI Level 1 Method: Theory and Implementation, Kamath Chandrika, 2000"
        /// </summary>
        public const double PivotTolerance = 1e-7; //

        private readonly double[] values;
        private readonly int[] diagOffsets;
        private readonly List<double[]> nullSpaceBasis; //TODO: should these be sparse vectors?

        private SemidefiniteLdlSkyline(int order, double[] values, int[] diagOffsets, 
            List<int> dependentColumns, List<double[]> nullSpaceBasis)
        {
            this.Order = order;
            this.values = values;
            this.diagOffsets = diagOffsets;
            this.DependentColumns = dependentColumns;
            this.nullSpaceBasis = nullSpaceBasis;
        }

        public IReadOnlyList<int> DependentColumns { get; }

        public IReadOnlyList<double[]> NullSpaceBasis => nullSpaceBasis; //TODO: return IVectorView

        /// <summary>
        /// The number of rows/columns of the original square matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Applies the LDL factorization to the independent columns of a symmetric positive semi-definite matrix,
        /// sets the dependent ones equal to columns of the identity matrix and return the nullspace of the matrix. Requires 
        /// extra memory for the basis vectors of the nullspace.
        /// </summary>
        /// <param name="order">The number of rows/ columns of the square matrix.</param>
        /// <param name="skyValues">
        /// The non-zero entries of the original <see cref="SkylineMatrix"/>. This array will be overwritten during the 
        /// factorization.
        /// </param>
        /// <param name="skyDiagOffsets">
        /// The indexes of the diagonal entries into <paramref name="skyValues"/>. The new 
        /// <see cref="SemidefiniteLdlSkyline"/> instance will hold a reference to <paramref name="skyDiagOffsets"/>. 
        /// However they do not need copying, since they will not be altered during or after the factorization.
        /// </param>
        /// <param name="pivotTolerance">
        /// If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the corresponding column is dependent 
        /// on the rest. The Cholesky factorization only applies to independent column, while dependent ones are used to compute
        /// the nullspace. Therefore it is important to select a tolerance that will identify small pivots that result from 
        /// singularity, but not from ill-conditioning.
        /// </param>
        public static SemidefiniteLdlSkyline Factorize(int order, double[] skyValues, int[] skyDiagOffsets,
            double pivotTolerance = SemidefiniteLdlSkyline.PivotTolerance)
        {
            (List<int> dependentColumns, List<double[]> nullSpaceBasis) = 
                LdlSkyline.FactorizeInternal(order, skyValues, skyDiagOffsets, pivotTolerance);
            Debug.Assert(dependentColumns.Count == nullSpaceBasis.Count);
            return new SemidefiniteLdlSkyline(order, skyValues, skyDiagOffsets, dependentColumns, nullSpaceBasis);
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        public double CalcDeterminant()
        {
            if (DependentColumns.Count > 0) return 0.0;
            else throw new NotImplementedException(); //TODO: call the code of the regular CholeskyLdl
        }

        /// <summary>
        /// Performs the operation: <paramref name="result"/> = generalized_inverse(A) * <paramref name="vector"/>
        /// </summary>
        /// <param name="vector">The vector that will be multiplied. Its <see cref="IIndexable1D.Length"/> must be equal to 
        /// <see cref="IIndexable2D.NumRows"/> of the original matrix A.
        /// </param>
        /// <param name="result">
        /// Output vector that will be overwritten with the solution of the linear system. Its <see cref="IIndexable1D.Length"/>  
        /// must be equal to <see cref="IIndexable2D.NumColumns"/> of the original matrix A.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="vector"/> or <paramref name="result"/> violate the described constraints.
        /// </exception>
        /// <exception cref="IndefiniteMatrixException">
        /// Thrown if the original skyline matrix turns out to not be symmetric positive semi-definite.
        /// </exception>
        public void MultiplyGeneralizedInverseMatrixTimesVector(Vector vector, Vector result)
        {
            Preconditions.CheckSystemSolutionDimensions(Order, vector.Length);
            Preconditions.CheckMultiplicationDimensions(Order, result.Length);

            LdlSkyline.Solve(Order, values, diagOffsets, vector.RawData, result.RawData);
        }
    }
}
