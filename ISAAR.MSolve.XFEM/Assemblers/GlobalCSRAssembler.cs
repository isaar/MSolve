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
        public (DOKRowMajor Kff, CSRMatrix Kfc) BuildGlobalMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsConstrained = dofOrderer.NumConstrainedDofs;
            int numDofsFree = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;
            var Kff = DOKRowMajor.CreateEmpty(numDofsFree, numDofsFree);
            var Kfc = DOKRowMajor.CreateEmpty(numDofsFree, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrix and add its contributions to the global matrices.
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kff.AddSubmatrix(kss, mapStandard, mapStandard);
                Kfc.AddSubmatrix(kss, mapStandard, mapConstrained);

                // Build enriched element matrices and add their contributions to the global matrices
                IReadOnlyDictionary<int, int> mapEnriched = dofOrderer.MatchElementToGlobalEnrichedDofsOf(element);
                if (mapEnriched.Count > 0)
                {
                    (Matrix kee, Matrix kes) = element.BuildEnrichedStiffnessMatricesLower();
                    Matrix kse = kes.Transpose();
                    Kff.AddSubmatrix(kes, mapEnriched, mapStandard);
                    Kff.AddSubmatrix(kse, mapStandard, mapEnriched);
                    Kff.AddSubmatrix(kee, mapEnriched, mapEnriched);
                    Kfc.AddSubmatrix(kes, mapEnriched, mapConstrained);
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
            return (Kff, Kfc.BuildCSRMatrix(true));
        }
    }
}
