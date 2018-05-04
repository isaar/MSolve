using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Tests.Tools;

//TODO: perhaps use a symmetric row major DOK and then build full CSR matrix.
namespace ISAAR.MSolve.XFEM.Assemblers
{
    class GlobalCSRAssembler
    {
        public (DOKRowMajor Kuu, CSRMatrix Kuc) BuildGlobalMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsConstrained = dofOrderer.ConstrainedDofsCount;
            int numDofsUnconstrained = dofOrderer.StandardDofsCount + dofOrderer.EnrichedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            DOKRowMajor Kuu = DOKRowMajor.CreateEmpty(numDofsUnconstrained, numDofsUnconstrained);

            // TODO: perhaps I should return a CSC matrix and do the transposed multiplication. This way I will not have to 
            // transpose the element matrix. Another approach is to add an AddTransposed() method to the DOK.
            var Kuc = DOKRowMajor.CreateEmpty(numDofsUnconstrained, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrices and add it contributions to the global matrices
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapFree, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kuu.AddSubmatrix(kss, mapFree, mapFree);
                Kuc.AddSubmatrix(kss, mapFree, mapConstrained);

                // Build enriched element matrices and add it contributions to the global matrices
                IReadOnlyDictionary<int, int> mapEnriched = dofOrderer.MatchElementToGlobalEnrichedDofsOf(element);
                if (mapEnriched.Count > 0)
                {
                    element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                    // TODO: options: 1) Only work with upper triangle in all symmetric matrices. Same applies to Elements.
                    // 2) The Elements have two versions of BuildStiffness(). 
                    // 3) The Elements return both (redundant; If someone needs it he can make it himself like here) 
                    Matrix kse = kes.Transpose();
                    Kuu.AddSubmatrix(kse, mapFree, mapEnriched);
                    Kuu.AddSubmatrix(kes, mapEnriched, mapFree);
                    Kuu.AddSubmatrix(kee, mapEnriched, mapEnriched);
                    Kuc.AddSubmatrix(kes, mapEnriched, mapConstrained);
                }                
            }
            #region DEBUG code
            //var checker = new SubmatrixChecker();
            //(Matrix expectedKuu, Matrix expectedKuc) = GlobalDenseAssembler.BuildGlobalMatrix(model, dofOrderer);
            //Console.WriteLine("Check Kuu:");
            //checker.Check(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //checker.Check(expectedKuc, Kuc);
            #endregion
            return (Kuu, Kuc.BuildCSRMatrix(true));
        }
    }
}
