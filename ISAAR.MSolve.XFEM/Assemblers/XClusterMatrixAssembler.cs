using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    class XClusterMatrixAssembler
    {
        public (DOKRowMajor Kss, DOKRowMajor Ksc) BuildStandardMatrices(Model2D model, XClusterDofOrderer globalDofOrderer)
        {
            int numDofsConstrained = globalDofOrderer.NumConstrainedDofs;
            int numDofsStandard = globalDofOrderer.NumStandardDofs;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            DOKRowMajor Kss = DOKRowMajor.CreateEmpty(numDofsStandard, numDofsStandard);

            // TODO: perhaps I should return a CSC matrix and do the transposed multiplication. This way I will not have to 
            // transpose the element matrix. Another approach is to add an AddTransposed() method to the DOK.
            var Ksc = DOKRowMajor.CreateEmpty(numDofsStandard, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrices and add it contributions to the global matrices
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                globalDofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kss.AddSubmatrix(kss, mapStandard, mapStandard);
                Ksc.AddSubmatrix(kss, mapStandard, mapConstrained);
            }

            return (Kss, Ksc);
        }

        public (DOKSymmetricColMajor Kee, DOKRowMajor Kes, DOKRowMajor Kec) BuildSubdomainMatrices(XSubdomain2D subdomain, 
            XClusterDofOrderer globalDofOrderer)
        {
            int numDofsEnriched = subdomain.DofOrderer.NumEnrichedDofs;
            int numDofsStandard = globalDofOrderer.NumStandardDofs;
            int numDofsConstrained = globalDofOrderer.NumConstrainedDofs;

            var Kee = DOKSymmetricColMajor.CreateEmpty(numDofsEnriched);
            var Kes = DOKRowMajor.CreateEmpty(numDofsEnriched, numDofsStandard);
            var Kec = DOKRowMajor.CreateEmpty(numDofsEnriched, numDofsConstrained);

            foreach (XContinuumElement2D element in subdomain.Elements)
            {
                // Build enriched element matrices and add their contributions to the global matrices
                Dictionary<int, int> enrichedMap = subdomain.DofOrderer.MatchElementToSubdomainEnrichedDofs(element);
                globalDofOrderer.MatchElementToGlobalStandardDofsOf(element, out IReadOnlyDictionary<int, int> standardMap,
                    out IReadOnlyDictionary<int, int> constrainedMap);

                // Not all elements are enriched necessarily. The domain decomposition might be done only at the start.
                if (enrichedMap.Count > 0)
                {
                    element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                    Kee.AddSubmatrixSymmetric(kee, enrichedMap);
                    Kes.AddSubmatrix(kes, enrichedMap, standardMap);
                    Kec.AddSubmatrix(kes, enrichedMap, constrainedMap);
                }
            }

            return (Kee, Kes, Kec);
        }
    }
}
