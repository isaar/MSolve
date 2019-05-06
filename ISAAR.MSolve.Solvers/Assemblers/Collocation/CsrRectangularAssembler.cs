using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Solvers.Commons;

//TODO: Instead of storing the raw CSR arrays, use a reusable DOK or CsrIndexer class. That class should provide methods to 
//      assemble the values part of the global matrix more efficiently than the general purpose DOK. The general purpose DOK 
//      should only be used to assemble the first global matrix and whenever the dof ordering changes. Now it is used everytime 
//      and the indexing arrays are discarded.
namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is square and stored in CSR format, but
    /// both triangles are explicitly stored. This format is suitable for matrix/vector multiplications, therefore it can be 
    /// combined with many iterative solvers. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CsrRectangularAssembler : IGlobalMatrixRectangularAssembler<CsrMatrix>
    {
        private const string name = "CsrRectangularAssembler"; // for error messages
        private readonly bool sortColsOfEachRow;
        private ConstrainedMatricesAssembler constrainedAssembler = new ConstrainedMatricesAssembler();

        bool isIndexerCached = false;
        private int[] cachedColIndices, cachedRowOffsets;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sortColsOfEachRow">
        /// Sorting the columns of each row in the CSR storage format may increase performance of the matrix vector 
        /// multiplications. It is recommended to set it to true, especially for iterative linear system solvers.
        /// </param>
        public CsrRectangularAssembler(bool sortColsOfEachRow = true)
        {
            this.sortColsOfEachRow = sortColsOfEachRow;
        }

        public CsrMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofRowOrdering,
            ISubdomainFreeDofOrdering dofColOrdering, IEnumerable<IElement> elements,
            IElementMatrixProvider matrixProvider)
        {
            var subdomainMatrix = DokRowMajor.CreateEmpty(dofRowOrdering.NumFreeDofs, dofColOrdering.NumFreeDofs);

            foreach (IElement element in elements)
            {
                (int[] elementColumnIndices, int[] subdomainColumnIndices) =
                    dofColOrdering.MapFreeDofsElementToSubdomain(element);
                (int[] elementRowIndices, int[] subdomainRowIndices) =
                    dofRowOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrix(elementMatrix, elementRowIndices, subdomainRowIndices,
                    elementColumnIndices, subdomainColumnIndices);
            }

            (double[] values, int[] colIndices, int[] rowOffsets) = subdomainMatrix.BuildCsrArrays(sortColsOfEachRow);

            if (!isIndexerCached)
            {
                cachedColIndices = colIndices;
                cachedRowOffsets = rowOffsets;
                isIndexerCached = true;
            }
            else
            {
                Debug.Assert(Utilities.AreEqual(cachedColIndices, colIndices));
                Debug.Assert(Utilities.AreEqual(cachedRowOffsets, rowOffsets));
            }

            return CsrMatrix.CreateFromArrays(dofRowOrdering.NumFreeDofs, dofColOrdering.NumFreeDofs, values,
                cachedColIndices, cachedRowOffsets, false);
        }

        public (CsrMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree, IMatrixView
            matrixConstrConstr) BuildGlobalSubmatrices(ISubdomainFreeDofOrdering freeDofRowOrdering,
                ISubdomainFreeDofOrdering freeDofColOrdering,
                ISubdomainConstrainedDofOrdering constrainedDofRowOrdering,
                ISubdomainConstrainedDofOrdering constrainedDofColOrdering, IEnumerable<IElement> elements,
                IElementMatrixProvider matrixProvider)
        {
            int numFreeRowDofs = freeDofRowOrdering.NumFreeDofs;
            int numFreeColDofs = freeDofColOrdering.NumFreeDofs;

            var subdomainMatrix = DokRowMajor.CreateEmpty(numFreeRowDofs, numFreeColDofs);

            constrainedAssembler.InitializeNewMatrices(freeDofRowOrdering.NumFreeDofs,
                constrainedDofRowOrdering.NumConstrainedDofs);

            foreach (var element in elements)
            {
                (int[] elementRowDofsFree, int[] subdomainRowDofsFree) =
                    freeDofRowOrdering.MapFreeDofsElementToSubdomain(element);
                (int[] elementColDofsFree, int[] subdomainColDofsFree) =
                    freeDofColOrdering.MapFreeDofsElementToSubdomain(element);

                (int[] elementRowDofsConstrained, int[] subdomainRowDofsConstrained) =
                    constrainedDofRowOrdering.MapConstrainedDofsElementToSubdomain(element);
                (int[] elementColDofsConstrained, int[] subdomainColDofsConstrained) =
                    constrainedDofColOrdering.MapConstrainedDofsElementToSubdomain(element);

                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrix(elementMatrix, elementRowDofsFree, subdomainRowDofsFree,
                    elementColDofsFree, subdomainColDofsFree);
                constrainedAssembler.AddElementMatrix(elementMatrix, elementRowDofsFree, subdomainRowDofsFree,
                    elementRowDofsConstrained, subdomainRowDofsConstrained); //TODO: check validity
            }

            (double[] values, int[] colIndices, int[] rowOffsets) = subdomainMatrix.BuildCsrArrays(sortColsOfEachRow);
            if (!isIndexerCached)
            {
                cachedColIndices = colIndices;
                cachedRowOffsets = rowOffsets;
                isIndexerCached = true;
            }
            else
            {
                Debug.Assert(Utilities.AreEqual(cachedColIndices, colIndices));
                Debug.Assert(Utilities.AreEqual(cachedRowOffsets, rowOffsets));
            }

            subdomainMatrix = null;
            var matrixFreeFree =
                CsrMatrix.CreateFromArrays(numFreeRowDofs, numFreeColDofs, values, cachedColIndices, cachedRowOffsets,
                    false);
            (CsrMatrix matrixConstrFree, CsrMatrix matrixConstrConstr) =
                constrainedAssembler.BuildMatrices(); // TODO: see if this work
            return (matrixFreeFree, matrixConstrFree, matrixConstrFree.TransposeToCSC(false), matrixConstrConstr);
        }

        public void HandleDofOrderingWillBeModified()
        {
            //TODO: perhaps the indexer should be disposed altogether. Then again it could be in use by other matrices.
            cachedColIndices = null;
            cachedRowOffsets = null;
            isIndexerCached = false;
        }
    }
}