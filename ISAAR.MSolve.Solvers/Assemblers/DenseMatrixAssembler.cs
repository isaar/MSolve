using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers
{
    public class DenseMatrixAssembler
    {
        public IMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements, 
            IElementMatrixProvider elementMatrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var Kff = Matrix.CreateZero(numFreeDofs, numFreeDofs);

            foreach (IElement element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                //IReadOnlyDictionary<int, int> elementToGlobalDofs = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix2D elementK = elementMatrixProvider.Matrix(element);
                //AddElementToGlobalMatrix(Kff, elementK, mapStandard, mapStandard);
                AddElementToGlobalMatrix(Kff, elementK, elementDofIndices, subdomainDofIndices);
            }

            return Kff;
        }

        private static void AddElementToGlobalMatrix(Matrix globalMatrix, IMatrix2D elementMatrix,
            int[] elementIndices, int[] globalIndices)
        {
            Debug.Assert(elementMatrix.Rows == elementMatrix.Columns);
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

        private static void AddElementToGlobalMatrix(Matrix globalMatrix, IMatrix2D elementMatrix,
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
