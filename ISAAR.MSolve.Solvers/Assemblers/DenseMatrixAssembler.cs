using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers
{
    public class DenseMatrixAssembler
    {
        public IMatrix BuildGlobalMatrix(IDofOrderer dofOrderer, IEnumerable<IElement> elements, 
            IElementMatrixProvider elementMatrixProvider)
        {
            int numFreeDofs = dofOrderer.NumFreeDofs;
            var Kff = Matrix.CreateZero(numFreeDofs, numFreeDofs);

            foreach (IElement element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                IReadOnlyDictionary<int, int> mapStandard = dofOrderer.MapFreeDofsElementToGlobal(element);
                IMatrix2D elementK = elementMatrixProvider.Matrix(element);
                AddElementToGlobalMatrix(Kff, elementK, mapStandard, mapStandard);
            }

            return Kff;
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
