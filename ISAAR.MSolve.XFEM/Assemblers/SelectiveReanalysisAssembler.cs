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
// WARNING: when other types of boundary conditions are needed (e.g. distributed loads, displacements) then the creationf rhs
//      must be updated.
namespace ISAAR.MSolve.XFEM.Assemblers
{
    class SelectiveReanalysisAssembler
    {
        private readonly int order;
        private Vector rhs;

        public SelectiveReanalysisAssembler(int order)
        {
            this.order = order;
        }

        //TODO: should I return it without copying it? Actually that needs SolveLinearSystem() that use IVectorView and do not 
        //      mutate it.
        public Vector Rhs { get { return rhs.Copy(); } } 

        /// <summary>
        /// 
        /// </summary>
        /// <param name="dofOrderer"></param>
        /// <param name="modifiedElements"></param>
        /// <param name="wantedColumns">Should not include the columns that will be identity.</param>
        public PartialMatrixColumns BuildGlobalMatrixColumns(IDofOrderer dofOrderer, ISet<int> wantedColumns, 
            IEnumerable<XContinuumElement2D> modifiedElements)
        {
            //TODO: also update Rhs if needed

            var globalK = new PartialMatrixColumns(order);
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
                    element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);
                    Matrix kse = kes.Transpose(); //TODO: use a version of BuildEnrichedStiffnessMatrices() that returns kee, kse
                    AddElementStiffnessToGlobal(wantedColumns, globalK, kee, kse, mapEnriched, mapStandard);

                    //TODO: also update Rhs if needed
                }
            }

            return globalK;
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
    }
}
