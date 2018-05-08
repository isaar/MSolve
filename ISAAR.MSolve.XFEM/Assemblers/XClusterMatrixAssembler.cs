using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
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

        // Create the signed boolean matrix/matrices to ensure continuity of displacements. 
        public SignedBooleanMatrix BuildGlobalSignedBooleanMatrix(XCluster2D cluster)
        {
            // Find the number of rows = number of continuity equations. TODO: It could be computed while each row is processed
            int numContinuityEquations = 0;
            var nodeMembership = cluster.FindEnrichedBoundaryNodeMembership();
            foreach (var nodeSubdomains in nodeMembership)
            {
                int nodeMultiplicity = nodeSubdomains.Value.Count;
                int numEnrichedDofs = nodeSubdomains.Key.EnrichedDofsCount;
                numContinuityEquations += (nodeMultiplicity - 1) * numEnrichedDofs;
            }

            //TODO: numCols could have been computed during matrix subdomain matrix assembly. numRows
            int numColumns = cluster.DofOrderer.NumStandardDofs + cluster.DofOrderer.NumEnrichedDofs; // Not sure about the std dofs
            var signedBoolean = new SignedBooleanMatrix(numContinuityEquations, numColumns);

            int globalEquation = 0; // index of continuity equation onto the global signed boolean matrix
            foreach (var nodeSubdomains in nodeMembership)
            {
                XNode2D node = nodeSubdomains.Key;

                // All enriched dofs of this node will have the same [1 -1] pattern in B.
                XSubdomain2D[] subdomains = nodeSubdomains.Value.ToArray();
                var numEquationsPerDof = subdomains.Length - 1;
                var positiveSubdomains = new XSubdomain2D[numEquationsPerDof];
                var negativeSubdomains = new XSubdomain2D[numEquationsPerDof];
                for (int equation = 0; equation < numEquationsPerDof; ++equation)
                {
                    positiveSubdomains[equation] = subdomains[equation];
                    negativeSubdomains[equation] = subdomains[equation + 1];
                }

                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (var dof in enrichment.Dofs)
                    {
                        for (int equation = 0; equation < positiveSubdomains.Length; ++equation)
                        {
                            int row = globalEquation + equation;
                            int posCol = positiveSubdomains[equation].DofOrderer.GetGlobalEnrichedDofOf(node, dof);
                            int negCol = negativeSubdomains[equation].DofOrderer.GetGlobalEnrichedDofOf(node, dof);
                            signedBoolean.AddEntry(row, posCol, true);
                            signedBoolean.AddEntry(row, negCol, false);
                        }
                        globalEquation += numEquationsPerDof;
                    }
                }
            }

            return signedBoolean;
        }
    }
}
