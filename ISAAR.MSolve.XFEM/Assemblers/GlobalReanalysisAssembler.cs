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

//TODO: Consider adding the functionality of initializing the global matrix as identity in the GlobalDOKAssembler. At the time 
//      of writing, the only difference is in the line that created the DOK, which can easily be controlled with a flag.
namespace ISAAR.MSolve.XFEM.Assemblers
{
    /// <summary>
    /// The global matrix will also contain entries related to inactive enriched degrees of freedom. These entries are taken 
    /// from the corresponding rows and columns of the identity matrix, and thus do not affect the solution. They do increase 
    /// the order of the matrix and thus the lengths of the rhs and solution vectors. 
    /// </summary>
    class GlobalReanalysisAssembler
    {
        public (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) BuildGlobalMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsConstrained = dofOrderer.NumConstrainedDofs;
            int numDofsUnconstrained = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            var Kuu = DOKSymmetricColMajor.CreateIdentity(dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs);

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
                Kuu.AddSubmatrixSymmetric(kss, mapFree);
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
                    Kuu.AddSubmatrixSymmetric(kee, mapEnriched);
                    Kuc.AddSubmatrix(kes, mapEnriched, mapConstrained);
                }
            }
            #region DEBUG code
            //(Matrix expectedKuu, Matrix expectedKuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, dofOrderer);
            //Console.WriteLine("Check Kuu:");
            //CheckMatrix(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //CheckMatrix(expectedKuc, Kuc);
            #endregion
            return (Kuu, Kuc.BuildCSRMatrix(true));
        }
    }
}
