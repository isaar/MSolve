using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    /// <summary>
    /// The matrix that will be "inverted" is a unified Skyline matrix, where enriched dofs are numbered after all 
    /// standard dofs. 
    /// TODO: The enriched dof columns will have huge heights. A more sophisticated solver and matrix assembler are
    /// needed. Also the global constrained submatrix must be sparse.
    /// </summary>
    static class SingleGlobalSkylineAssembler
    {
        public static void BuildGlobalMatrix(Model2D model, 
            out SkylineMatrix2D globalMatrixLhsLhs, out Matrix2D globalMatrixLhsCon)
        {
            DOFEnumerator dofEnumerator = model.DofEnumerator;
            int constrainedDofsCount = dofEnumerator.ConstrainedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            globalMatrixLhsLhs = new SkylineMatrix2D(CalculateSkylineIndexer(model));

            // TODO: this should be in a sparse format. Only used for SpMV and perhaps transpose SpMV.
            // Row = standard free dofs + enriched dofs. Columns = standard constrained dofs. 
            globalMatrixLhsCon = new Matrix2D(dofEnumerator.FreeDofsCount + dofEnumerator.ArtificialDofsCount, 
                dofEnumerator.ConstrainedDofsCount);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Element matrices
                SymmetricMatrix2D elementMatrixStdStd = element.BuildStandardStiffnessMatrix();
                Matrix2D elementMatrixEnrStd;
                SymmetricMatrix2D elementMatrixEnrEnr;
                element.BuildEnrichedStiffnessMatrices(out elementMatrixEnrStd, out elementMatrixEnrEnr);

                // Element to global dofs mappings
                // TODO: perhaps that could be done during the assembly to avoid iterating over the dofs twice
                IReadOnlyDictionary<int, int> elementToGlobalFreeDofs, elementToGlobalConstrainedDofs;
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element, 
                    out elementToGlobalFreeDofs, out elementToGlobalConstrainedDofs);
                IReadOnlyDictionary<int, int> elementToGlobalEnrDofs = 
                    dofEnumerator.MatchElementToGlobalArtificialDofsOf(element);

                // Add the element contributions to the global matrices
                AddElementToGlobalMatrix(globalMatrixLhsLhs, elementMatrixStdStd, elementToGlobalFreeDofs, elementToGlobalFreeDofs);
                AddElementToGlobalMatrix(globalMatrixLhsLhs, elementMatrixEnrStd, elementToGlobalEnrDofs, elementToGlobalFreeDofs);
                AddElementToGlobalMatrix(globalMatrixLhsLhs, elementMatrixEnrEnr, elementToGlobalEnrDofs, elementToGlobalEnrDofs);

                AddElementToGlobalMatrix(globalMatrixLhsCon, elementMatrixStdStd, elementToGlobalFreeDofs, elementToGlobalConstrainedDofs);
                AddElementToGlobalMatrix(globalMatrixLhsCon, elementMatrixEnrStd, elementToGlobalEnrDofs, elementToGlobalConstrainedDofs);
            }
        }

        private static int[] CalculateSkylineIndexer(Model2D model)
        {
            int[] rowHeights = CalculateRowHeights(model);

            int[] rowIndex = new int[rowHeights.Length + 1]; // The indexer uses 1 extra ghost entry to facilitate some operations
            rowIndex[0] = 0; // The diagonal entries of the first 2 rows always have the same positions
            rowIndex[1] = 1;
            for (int i = 1; i < rowHeights.Length; ++i)
                rowIndex[i + 1] = rowIndex[i] + rowHeights[i] + 1; //+1 because the heights don't include the diagonal
            return rowIndex;
        }

        /// Only the entries above the diagonal are considered
        private static int[] CalculateRowHeights(Model2D model) // I vote to call them RowWidths!!!
        {
            DOFEnumerator dofEnumerator = model.DofEnumerator;
            int[] rowHeights = new int[dofEnumerator.FreeDofsCount + dofEnumerator.ArtificialDofsCount];
            foreach (XContinuumElement2D element in model.Elements)
            {
                // To determine the row height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0 entries in the element stiffness matrix.
                int minDOF = Int32.MaxValue;
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofEnumerator.GetFreeDofsOf(node)) minDOF = Math.Min(dof, minDOF);
                    foreach (int dof in dofEnumerator.GetArtificialDofsOf(node)) Math.Min(dof, minDOF);
                }

                // The height of each row is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofEnumerator.GetFreeDofsOf(node))
                    {
                        rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                    }
                    foreach (int dof in dofEnumerator.GetArtificialDofsOf(node))
                    {
                        rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                    }
                }
            }
            return rowHeights;
        }

        private static void AddElementToGlobalMatrix(SkylineMatrix2D globalMatrix, IMatrix2D elementMatrix,
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

                    int height = globalRow - globalCol; // Skyline matrix stores the lower triangle
                    // Check if the entry is in the lower triangle. This check concerns global dofs only.
                    // No need to check if the entry is within the row height.
                    // TODO: This check is redundant for the enriched-standard submatrices, since they always 
                    // have row=enrichedDof > col=standardDof. Perhaps I could have an optimized method that 
                    // doesn't check that height>0. Frankly it would be too much code just for avoiding this 
                    // check. Perhaps after optimizing other things first.
                    if (height >= 0)
                    {
                        // Shouldn't the skyline indexing logic be abstracted?
                        int skylinePos = globalMatrix.RowIndex[globalRow] + height;
                        globalMatrix.Data[skylinePos] += elementMatrix[elementRow, elementCol];
                    }
                }
            }
        }

        // TODO: The global matrix should be sparse. Probably in CSR/DOK format
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
