using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is full column major format. It can be 
    /// symmetric, but both triangles will be stored explicitly.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DenseMatrixAssembler: IGlobalMatrixAssembler<Matrix>
    {
        public Matrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement_v2> elements, 
            IElementMatrixProvider_v2 elementMatrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var subdomainMatrix = Matrix.CreateZero(numFreeDofs, numFreeDofs);

            foreach (IElement_v2 element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                //IReadOnlyDictionary<int, int> elementToGlobalDofs = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementMatrix = elementMatrixProvider.Matrix(element);
                AddElementToGlobalMatrix(subdomainMatrix, elementMatrix, elementDofIndices, subdomainDofIndices);
            }

            return subdomainMatrix;
        }

        public void HandleDofOrderingWillBeModified()
        {
           // Do nothing, since there are no idexing arrays to cache.
        }

        private static void AddElementToGlobalMatrix(Matrix globalMatrix, IMatrixView elementMatrix,
            int[] elementIndices, int[] globalIndices)
        {
            Debug.Assert(elementMatrix.NumRows == elementMatrix.NumColumns);
            Debug.Assert(globalIndices.Length == elementIndices.Length);

            int numRelevantRows = elementIndices.Length;
            for (int i = 0; i < numRelevantRows; ++i)
            {
                int elementRow = elementIndices[i];
                int globalRow = globalIndices[i];
                for (int j = 0; j < numRelevantRows; ++j)
                {
                    int elementCol = elementIndices[j];
                    int globalCol = globalIndices[j];

                    globalMatrix[globalRow, globalCol] += elementMatrix[elementRow, elementCol];
                }
            }
        }

        private static void AddElementToGlobalMatrix(Matrix globalMatrix, IMatrixView elementMatrix,
            IReadOnlyDictionary<int, int> elementRowsToGlobalRows, IReadOnlyDictionary<int, int> elementColsToGlobalCols)
        {
            foreach (var rowPair in elementRowsToGlobalRows)
            {
                int elementRow = rowPair.Key;
                int globalRow = rowPair.Value;
                foreach (var colPair in elementColsToGlobalCols)
                {
                    int elementCol = colPair.Key;
                    int globalCol = colPair.Value;

                    globalMatrix[globalRow, globalCol] += elementMatrix[elementRow, elementCol];
                }
            }
        }
    }
}
