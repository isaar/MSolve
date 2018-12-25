using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.Assemblers
{
    public class DenseMatrixAssembler
    {
        public IMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement_v2> elements, 
            IElementMatrixProvider_v2 elementMatrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var Kff = Matrix.CreateZero(numFreeDofs, numFreeDofs);

            foreach (IElement_v2 element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                //IReadOnlyDictionary<int, int> elementToGlobalDofs = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementK = elementMatrixProvider.Matrix(element);
                //AddElementToGlobalMatrix(Kff, elementK, mapStandard, mapStandard);
                AddElementToGlobalMatrix(Kff, elementK, elementDofIndices, subdomainDofIndices);
            }

            return Kff;
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
