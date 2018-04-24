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

/// <summary>
/// The matrix that will be "inverted" is a unified DOK matrix, where enriched dofs are numbered after all 
/// standard dofs. 
/// TODO: The enriched dof columns will have huge heights. A more sophisticated solver and matrix assembler are
/// needed. Also the global constrained submatrix must be sparse.
/// </summary>
namespace ISAAR.MSolve.XFEM.Assemblers
{
    class GlobalDOKAssembler
    {
        public (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) BuildGlobalMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int numDofsConstrained = dofEnumerator.ConstrainedDofsCount;
            int numDofsUnconstrained = dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            var Kuu = DOKSymmetricColMajor.CreateEmpty(dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount);

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
    }
}
