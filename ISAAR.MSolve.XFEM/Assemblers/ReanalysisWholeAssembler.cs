using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Elements;
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
    class ReanalysisWholeAssembler
    {
        public (DokSymmetric Kff, DokRowMajor Kfc) BuildGlobalMatrix(IEnumerable<XContinuumElement2D> allElements, 
            IDofOrderer dofOrderer)
        {
            int numDofsConstrained = dofOrderer.NumConstrainedDofs;
            int numDofsFree = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;
            var Kff = DokSymmetric.CreateEmpty(numDofsFree);
            var Kfc = DokRowMajor.CreateEmpty(numDofsFree, numDofsConstrained);

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
    }
}
