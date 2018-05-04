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
    class XSubdomainMatrixAssembler
    {
        public DOKSymmetricColMajor BuildSubdomainMatrices(XSubdomain2D subdomain)
        {
            XSubdomainDofOrderer dofOrder = subdomain.DofOrder;

            int numDofsStd = subdomain.Nodes.Count * 2;
            int numDofsEnr = dofOrder.NumEnrichedDofs;

            var Kee = DOKSymmetricColMajor.CreateEmpty(numDofsEnr);
            //var Kes = DOKRowMajor.CreateEmpty(numDofsEnr, numDofsStd);

            foreach (XContinuumElement2D element in subdomain.Elements)
            {
                // Build enriched element matrices and add their contributions to the global matrices
                IReadOnlyDictionary<int, int> mapEnriched = dofOrder.MatchElementToSubdomainEnrichedDofs(element);

                // Not all elements are enriched necessarily. The domain decomposition might be done only at the start.
                if (mapEnriched.Count > 0) 
                {
                    element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                    Kee.AddSubmatrixSymmetric(kee, mapEnriched);
                    //Kes.AddSubmatrix(kes, mapEnriched, mapFree);
                }
            }

            return Kee;
            //return (Kee, Kes);
        }
    }
}
