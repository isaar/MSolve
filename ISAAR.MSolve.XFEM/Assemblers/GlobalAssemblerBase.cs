using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

// TODO: No need for inheritance here. Just different methods for different matrix types. Although that might violate Open-Closed
namespace ISAAR.MSolve.XFEM.Assemblers
{
    /// <summary>
    /// It assumes that element matrices are symmetric!!!
    /// </summary>
    abstract class GlobalAssemblerBase<TMatrix> where TMatrix: ISymmetricMatrixBuilder //TODO: do not abstract the returned types
    {
        public (TMatrix Kuu, CSRMatrix Kuc) BuildGlobalMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int numDofsConstrained = dofEnumerator.ConstrainedDofsCount;
            int numDofsUnconstrained = dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            TMatrix Kuu = InitializeGlobalUncontrainedMatrix(model, dofEnumerator);

            // TODO: perhaps I should return a CSC matrix and do the transposed multiplication. This way I will not have to 
            // transpose the element matrix. Another approach is to add an AddTransposed() method to the DOK.
            var Kuc = DOKRowMajor.CreateEmpty(numDofsUnconstrained, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Element matrices
                Matrix kss = element.BuildStandardStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                // TODO: options: 1) Only work with upper triangle in all symmetric matrices. Same applies to Elements.
                // 2) The Elements have two versions of BuildStiffness(). 
                // 3) The Elements return both (redundant; If someone needs it he can make it himself like here) 
                Matrix kse = kes.Transpose();

                // Element to global dofs mappings
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapFree, out IReadOnlyDictionary<int, int> mapConstrained);
                IReadOnlyDictionary<int, int> mapEnriched =
                    dofEnumerator.MatchElementToGlobalEnrichedDofsOf(element);

                // Add the element contributions to the global matrices
                Kuu.AddSubmatrixSymmetric(kss, mapFree);
                Kuu.AddSubmatrix(kse, mapFree, mapEnriched);
                Kuu.AddSubmatrixSymmetric(kee, mapEnriched);

                Kuc.AddSubmatrix(kss, mapFree, mapConstrained);
                Kuc.AddSubmatrix(kes, mapEnriched, mapConstrained);
            }

            #region DEBUG code
            //(Matrix expectedKuu, Matrix expectedKuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, dofEnumerator);
            //Console.WriteLine("Check Kuu:");
            //CheckMatrix(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //CheckMatrix(expectedKuc, Kuc);
            #endregion

            //TODO: perhaps I should filter the matrices in the concrete class before returning (e.g. sok.Build())
            return (Kuu, Kuc.BuildCSRMatrix(true)); 
        }

        protected abstract TMatrix InitializeGlobalUncontrainedMatrix(Model2D model, IDOFEnumerator dofEnumerator);

        // TODO: The global matrix should be sparse. Probably in CSR/DOK format
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

        private static void CheckMatrix(IIndexable2D expected, IIndexable2D computed)
        {
            bool isCorrect = true;
            ValueComparer comparer = new ValueComparer(1e-6);
            if ((computed.NumRows != expected.NumRows) || (computed.NumColumns != expected.NumColumns))
            {
                Console.WriteLine("Invalid dimensions");
            }
            for (int i = 0; i < computed.NumRows; ++i)
            {
                for (int j = 0; j < computed.NumColumns; ++j)
                {
                    if (!comparer.AreEqual(computed[i, j], expected[i, j]))
                    {
                        Console.WriteLine($"Computed[{i}, {j}] = {computed[i, j]}   -   Expected[{i}, {j}] = {expected[i, j]}");
                        isCorrect = false;
                    }
                }
            }
            if (isCorrect) Console.WriteLine("Correct");
        }
    }
}
