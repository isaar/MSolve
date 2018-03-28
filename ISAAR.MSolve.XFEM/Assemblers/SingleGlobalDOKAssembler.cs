using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

/// <summary>
/// The matrix that will be "inverted" is a unified DOK matrix, where enriched dofs are numbered after all 
/// standard dofs. 
/// TODO: The enriched dof columns will have huge heights. A more sophisticated solver and matrix assembler are
/// needed. Also the global constrained submatrix must be sparse.
/// </summary>
namespace ISAAR.MSolve.XFEM.Assemblers
{
    static class SingleGlobalDOKAssembler
    {
        public static void BuildGlobalMatrix(Model2D model,
            out SymmetricDOKColMajor globalMatrixFreeFree, out Matrix2D globalMatrixFreeCons)
        {
            DOFEnumerator dofEnumerator = model.DofEnumerator;
            int constrainedDofsCount = dofEnumerator.ConstrainedDofsCount;
            int allFreeDofsCount = dofEnumerator.FreeDofsCount + dofEnumerator.ArtificialDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            globalMatrixFreeFree = new SymmetricDOKColMajor(allFreeDofsCount);

            // TODO: this should be in a sparse format. Only used for SpMV and perhaps transpose SpMV.
            // Row = standard free dofs + enriched dofs. Columns = standard constrained dofs. 
            globalMatrixFreeCons = new Matrix2D(dofEnumerator.FreeDofsCount + dofEnumerator.ArtificialDofsCount,
                dofEnumerator.ConstrainedDofsCount);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Element matrices
                SymmetricMatrix2D elementMatrixStdStd = element.BuildStandardStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out Matrix2D elementMatrixEnrStd, 
                    out SymmetricMatrix2D elementMatrixEnrEnr);

                // Element to global dofs mappings
                // TODO: perhaps that could be done during the assembly to avoid iterating over the dofs twice
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element, 
                    out IReadOnlyDictionary<int, int> elementToGlobalFreeDofs, 
                    out IReadOnlyDictionary<int, int> elementToGlobalConstrainedDofs);
                IReadOnlyDictionary<int, int> elementToGlobalEnrDofs =
                    dofEnumerator.MatchElementToGlobalArtificialDofsOf(element);

                // Add the element contributions to the global matrices
                AddElementToGlobalMatrix(globalMatrixFreeFree, elementMatrixStdStd, elementToGlobalFreeDofs, elementToGlobalFreeDofs);
                AddElementToGlobalMatrix(globalMatrixFreeFree, elementMatrixEnrStd, elementToGlobalEnrDofs, elementToGlobalFreeDofs);
                AddElementToGlobalMatrix(globalMatrixFreeFree, elementMatrixEnrEnr, elementToGlobalEnrDofs, elementToGlobalEnrDofs);

                AddElementToGlobalMatrix(globalMatrixFreeCons, elementMatrixStdStd, elementToGlobalFreeDofs, elementToGlobalConstrainedDofs);
                AddElementToGlobalMatrix(globalMatrixFreeCons, elementMatrixEnrStd, elementToGlobalEnrDofs, elementToGlobalConstrainedDofs);
            }
        }

        private static void AddElementToGlobalMatrix(SymmetricDOKColMajor globalMatrix, IMatrix2D elementMatrix,
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

                    // SymmetricDOKColMajor matrix stores the upper triangle
                    // Check if the entry is in the upper triangle. This check concerns global dofs only.
                    // No need to check if the entry is within the row height.
                    // TODO: This check is redundant for the enriched-standard submatrices, since they always 
                    // have row=enrichedDof > col=standardDof. Perhaps I could have an optimized method that 
                    // doesn't check that height>0. Frankly it would be too much code just for avoiding this 
                    // check. Perhaps after optimizing other things first.
                    if (globalCol >= globalRow)
                    {
                        globalMatrix.AddToEntry(globalRow, globalCol, elementMatrix[elementRow, elementCol]);
                    }
                }
            }
        }

        /// <summary>
        /// This method is for the non symmetric part of the matrix (free - constrained).
        /// TODO: The global matrix should be sparse. Probably in CSR/DOK format.
        /// </summary>
        /// <param name="globalMatrix"></param>
        /// <param name="elementMatrix"></param>
        /// <param name="elementRowsToGlobalRows"></param>
        /// <param name="elementColsToGlobalCols"></param>
        private static void AddElementToGlobalMatrix(Matrix2D globalMatrix, IMatrix2D elementMatrix,
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

                    globalMatrix[globalRow, globalCol] = elementMatrix[elementRow, elementCol];
                }
            }
        }
    }
}
