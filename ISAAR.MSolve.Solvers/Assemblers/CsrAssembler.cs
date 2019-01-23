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
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is square and stored in CSR format, but
    /// both triangles are explicitly stored. This format is suitable for matrix/vector multiplications, therefore it can be 
    /// combined with many iterative solvers. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CsrAssembler : IGlobalMatrixAssembler<CsrMatrix>
    {
        private const string name = "CsrAssembler"; // for error messages
        private readonly bool sortColsOfEachRow;

        bool isIndexerCached = false;
        private int[] cachedColIndices, cachedRowOffsets;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sortColsOfEachRow">
        /// Sorting the columns of each row in the CSR storage format may increase performance of the matrix vector 
        /// multiplications. It is recommended to set it to true, especially for iterative linear system solvers.
        /// </param>
        public CsrAssembler(bool sortColsOfEachRow = true)
        {
            this.sortColsOfEachRow = sortColsOfEachRow;
        }

        public CsrMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement_v2> elements, 
            IElementMatrixProvider_v2 matrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var subdomainMatrix = DokRowMajor.CreateEmpty(numFreeDofs, numFreeDofs);

            foreach (IElement_v2 element in elements)
            {
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrixSymmetric(elementMatrix, elementDofIndices, subdomainDofIndices);
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
            return CsrMatrix.CreateFromArrays(numFreeDofs, numFreeDofs, values, cachedColIndices, cachedRowOffsets, false);
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
