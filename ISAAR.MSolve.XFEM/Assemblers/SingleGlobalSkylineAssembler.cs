using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    /// <summary>
    /// A unified Skyline matrix, where enriched dofs are numbered after all standard dofs. The enriched dof columns 
    /// will have huge heights. A more sophisticated solver and matrix assembler are needed.
    /// </summary>
    static class SingleGlobalSkylineAssembler
    {
        public static SkylineMatrix2D BuildGlobalMatrix(Model2D model)
        {
            SkylineMatrix2D globalMatrix = new SkylineMatrix2D(CalculateSkylineIndexer(model));
            DOFEnumerator dofEnumerator = model.DofEnumerator;
            foreach (Element2D element in model.Elements)
            {
                // Element matrices
                SymmetricMatrix2D elementMatrixStdStd = element.ElementType.BuildStandardStiffnessMatrix();
                Matrix2D elementMatrixEnrStd;
                SymmetricMatrix2D elementMatrixEnrEnr;
                element.ElementType.BuildEnrichedStiffnessMatrices(out elementMatrixEnrStd, out elementMatrixEnrEnr);

                // Element to global dofs mappings
                IReadOnlyList<int> elementToGlobalStdDofs = dofEnumerator.MatchElementToGlobalStandardDofsOf(element);
                IReadOnlyList<int> elementToGlobalEnrDofs = dofEnumerator.MatchElementToGlobalArtificialDofsOf(element);

                // Add the element contributions to the global matrix
                AddElementToGlobalMatrix(globalMatrix, elementMatrixStdStd, elementToGlobalStdDofs, elementToGlobalStdDofs);
                AddElementToGlobalMatrix(globalMatrix, elementMatrixEnrStd, elementToGlobalEnrDofs, elementToGlobalStdDofs);
                AddElementToGlobalMatrix(globalMatrix, elementMatrixEnrEnr, elementToGlobalEnrDofs, elementToGlobalEnrDofs);
            }
            return globalMatrix;
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
            int[] rowHeights = new int[dofEnumerator.StandardDofsCount + dofEnumerator.ArtificialDofsCount];
            foreach (Element2D element in model.Elements)
            {
                // To determine the row height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0 entries in the element stiffness matrix.
                int minDOF = Int32.MaxValue;
                foreach (XNode2D node in element.ElementType.Nodes) // Should I draw the nodes from element.ElementType?
                {
                    // Both standard and artificial dofs are considered. I could have a DOFEnumerator.GetAllDofsOfNode() method.
                    foreach (int dof in dofEnumerator.GetStandardDofsOf(node))
                    {
                        if (dof != -1) minDOF = Math.Min(dof, minDOF);
                    }
                    foreach (int dof in dofEnumerator.GetArtificialDofsOf(node))
                    {
                        //if (dof != -1) minDOF = Math.Min(dof, minDOF); //Is it even possible to have constrained artificial dofs?
                        Math.Min(dof, minDOF);
                    }
                }

                // The height of each row is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (XNode2D node in element.ElementType.Nodes)
                {
                    foreach (int dof in dofEnumerator.GetStandardDofsOf(node))
                    {
                        if (dof != -1) rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                    }
                    foreach (int dof in dofEnumerator.GetArtificialDofsOf(node))
                    {
                        //if (dof != -1) rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF); //Is it even possible to have constrained artificial dofs?
                        rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                    }
                }
            }
            return rowHeights;
        }

        private static void AddElementToGlobalMatrix(
            SkylineMatrix2D globalMatrix, IMatrix2D elementMatrix,
            IReadOnlyList<int> elementRowsToGlobalRows, IReadOnlyList<int> elementColsToGlobalCols)
        {
            for (int elementRow = 0; elementRow < elementRowsToGlobalRows.Count; ++elementRow)
            {
                int globalRow = elementRowsToGlobalRows[elementRow];
                if (globalRow != -1)
                {
                    for (int elementCol = 0; elementCol < elementColsToGlobalCols.Count; ++elementCol)
                    {
                        int globalCol = elementColsToGlobalCols[elementCol];
                        if (globalCol != -1)
                        {
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
                
            }
        }
    }
}
