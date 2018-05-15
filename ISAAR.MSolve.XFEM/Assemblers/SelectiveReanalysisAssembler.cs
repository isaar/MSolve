using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: perhaps the matrix should be stored in a dedicated DOK class. DOKColMajor with Dictionary<int, Dictionary<int, double>>
//      as data seems like a good start.
// WARNING: when other types of boundary conditions are needed (e.g. distributed loads, displacements) then the rhs update
//      must be modified appropriately.
namespace ISAAR.MSolve.XFEM.Assemblers
{
    class SelectiveReanalysisAssembler
    {
        private readonly IDofOrderer dofOrderer;

        public SelectiveReanalysisAssembler(IDofOrderer dofOrderer)
        {
            this.dofOrderer = dofOrderer;
        }

        /// <summary>
        /// Build the whole stiffness matrices.
        /// </summary>
        /// <param name="allElements"></param>
        /// <returns></returns>
        public (DOKSymmetricColMajor Kff, DOKRowMajor Kfc) BuildGlobalMatrix(IEnumerable<XContinuumElement2D> allElements)
        {
            int numDofsConstrained = dofOrderer.NumConstrainedDofs;
            int numDofsFree = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;
            var Kff = DOKSymmetricColMajor.CreateEmpty(numDofsFree);
            var Kfc = DOKRowMajor.CreateEmpty(numDofsFree, numDofsConstrained);

            foreach (XContinuumElement2D element in allElements)
            {
                // Build standard element matrix and add its contributions to the global matrices
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kff.AddSubmatrixSymmetric(kss, mapStandard);
                Kfc.AddSubmatrix(kss, mapStandard, mapConstrained);

                // Build enriched element matrices and add their contributions to the global matrices
                IReadOnlyDictionary<int, int> mapEnriched = dofOrderer.MatchElementToGlobalEnrichedDofsOf(element);
                if (mapEnriched.Count > 0)
                {
                    (Matrix kee, Matrix kse) = element.BuildEnrichedStiffnessMatricesUpper();
                    Kff.AddSubmatrix(kse, mapStandard, mapEnriched);
                    Kff.AddSubmatrixSymmetric(kee, mapEnriched);
                    if (mapConstrained.Count > 0) Kfc.AddSubmatrix(kse.Transpose(), mapEnriched, mapConstrained);
                }
            }
            #region DEBUG code
            //(Matrix expectedKuu, Matrix expectedKuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, dofOrderer);
            //Console.WriteLine("Check Kuu:");
            //CheckMatrix(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //CheckMatrix(expectedKuc, Kuc);
            #endregion

            // Treat inactive/removed enriched dofs. I used to initialize the DOK to identity, but that was incorrect, since  
            // I am adding the non zero diagonals to 1.0, instead of replacing the 1.0.
            Kff.SetStructuralZeroDiagonalEntriesToUnity(); 
            return (Kff, Kfc);
        }

        /// <summary>
        /// Only build the columns that change.
        /// </summary>
        /// <param name="modifiedElements"></param>
        /// <param name="wantedColumns">Should not include the columns that will be identity.</param>
        /// <param name="previousRhs">It will be modified necessary.</param>
        /// <param name="prescribedNodalDisplacements"></param>
        /// <returns></returns>
        public PartialMatrixColumns BuildGlobalMatrixColumns(IEnumerable<XContinuumElement2D> modifiedElements, 
            ISet<int> wantedColumns, Vector previousRhs, Vector prescribedNodalDisplacements)
        {
            int numFreeDofs = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;
            var Kff = new PartialMatrixColumns(numFreeDofs);
            var Kfc = new PartialMatrixRows(numFreeDofs, dofOrderer.NumConstrainedDofs);
            bool isRhsChanged = false;

            foreach (XContinuumElement2D element in modifiedElements)
            {
                dofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                IReadOnlyDictionary<int, int> mapEnriched = dofOrderer.MatchElementToGlobalEnrichedDofsOf(element);

                // Is it possible that a standard element is included in the modifiedElements? If using fixed enrichment
                // area for the tip enrichments, then once it moves, some tip enriched/blending elements will be far enough from
                // the crack to avoid the Heaviside enrichment. However, should these elements be included in modifiedElements in
                // the first place? I do not care about the columns of their dofs as they will be identity. On the other hand it
                // doesn't cost much (only the enriched dofs map) and simplifies bookkeeping, as I would otherwise need yet 
                // another set of tracked elements.
                if (mapEnriched.Count > 0) 
                {
                    // Update modified stiffness matrix columns (free dofs)
                    (Matrix kee, Matrix kse) = element.BuildEnrichedStiffnessMatricesUpper();
                    AddElementStiffnessToGlobal(wantedColumns, Kff, kee, kse, mapEnriched, mapStandard);

                    // Update modified stiffness matrix rows for constrained dofs, if any
                    if (mapConstrained.Count > 0)
                    {
                        isRhsChanged = true;
                        AddElementStiffnessToGlobalConstrained(Kfc, kse.Transpose(), mapEnriched, mapConstrained);
                    }
                }
            }

            if (isRhsChanged)
            {
                Console.WriteLine("The rhs vector is modified. This should happen only if an element on the boundary has at"
                    + " least one node with changed enrichment or body level set");
                SparseVector rhsChanges = Kfc.MultiplyRight(prescribedNodalDisplacements);
                foreach ((int idx, double val) in rhsChanges.EnumerateNonZeros()) previousRhs[idx] = val;
            }

            return Kff;
        }

        //TODO: I process each enriched column in one loop, which saves some time. Provide this functionallity in the DOKs too
        //      and benchmark the gain.
        private void AddElementStiffnessToGlobal(ISet<int> wantedColumns,
            PartialMatrixColumns globalMatrix, Matrix elementMatrixEnrEnr, Matrix elementMatrixStdEnr,
            IReadOnlyDictionary<int, int> enrDofsElementToGlobal, IReadOnlyDictionary<int, int> stdDofsElementToGlobal)
        {
            foreach (var colPair in enrDofsElementToGlobal)
            {
                int elementCol = colPair.Key;
                int globalCol = colPair.Value;
                if (!wantedColumns.Contains(globalCol)) continue;

                foreach (var rowPair in enrDofsElementToGlobal)
                {
                    double elementValue = elementMatrixEnrEnr[rowPair.Key, elementCol];
                    globalMatrix.AddToEntry(rowPair.Value, globalCol, elementValue);
                }

                foreach (var rowPair in stdDofsElementToGlobal)
                {
                    double elementValue = elementMatrixStdEnr[rowPair.Key, elementCol];
                    globalMatrix.AddToEntry(rowPair.Value, globalCol, elementValue);
                }
            }
        }

        private void AddElementStiffnessToGlobalConstrained(PartialMatrixRows globalMatrixFreeConstr, Matrix elementMatrixEnrStd,
            IReadOnlyDictionary<int, int> enrDofsElementToGlobal, IReadOnlyDictionary<int, int> stdDofsElementToGlobal)
        {
            foreach (var colPair in stdDofsElementToGlobal)
            {
                int elementCol = colPair.Key;
                int globalCol = colPair.Value;

                // Some dofs miss from wantedCols, beacuse they are identity and are handled differently. Not applicable here. 
                //if (!wantedColumns.Contains(globalRow)) continue; 

                foreach (var rowPair in enrDofsElementToGlobal)
                {
                    double elementValue = elementMatrixEnrStd[rowPair.Key, elementCol];
                    globalMatrixFreeConstr.AddToEntry(rowPair.Value, globalCol, elementValue);
                }
            }
        }
    }
}
