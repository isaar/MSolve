using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    class XClusterMatrixAssembler
    {
        public (DokSymmetric Kss, DokRowMajor Ksc) BuildStandardMatrices(Model2D model, XClusterDofOrderer globalDofOrderer)
        {
            int numDofsConstrained = globalDofOrderer.NumConstrainedDofs;
            int numDofsStandard = globalDofOrderer.NumStandardDofs;
            var Kss = DokSymmetric.CreateEmpty(numDofsStandard);
            var Ksc = DokRowMajor.CreateEmpty(numDofsStandard, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrix and add its contributions to the global matrices
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                globalDofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kss.AddSubmatrixSymmetric(kss, mapStandard);
                Ksc.AddSubmatrix(kss, mapStandard, mapConstrained);
            }

            return (Kss, Ksc);
        }

        public (DokSymmetric Kee, DokRowMajor Kes, DokRowMajor Kec) BuildSubdomainMatrices(XSubdomain2D subdomain, 
            XClusterDofOrderer globalDofOrderer)
        {
            int numDofsConstrained = globalDofOrderer.NumConstrainedDofs;
            int numDofsEnriched = subdomain.DofOrderer.NumEnrichedDofs;
            int numDofsStandard = globalDofOrderer.NumStandardDofs;
            var Kee = DokSymmetric.CreateEmpty(numDofsEnriched);
            var Kes = DokRowMajor.CreateEmpty(numDofsEnriched, numDofsStandard);
            var Kec = DokRowMajor.CreateEmpty(numDofsEnriched, numDofsConstrained);

            foreach (XContinuumElement2D element in subdomain.Elements)
            {
                // Build enriched element matrices and add their contributions to the global matrices
                Dictionary<int, int> enrichedMap = subdomain.DofOrderer.MatchElementToSubdomainEnrichedDofs(element);
                
                // Not all elements are enriched necessarily. The domain decomposition might be done only at the start.
                if (enrichedMap.Count > 0)
                {
                    //TODO: This was already done in BuildStandardMatrices(). Find a way to reuse it.
                    globalDofOrderer.MatchElementToGlobalStandardDofsOf(element, 
                        out IReadOnlyDictionary<int, int> standardMap, out IReadOnlyDictionary<int, int> constrainedMap);

                    (Matrix kee, Matrix kes) = element.BuildEnrichedStiffnessMatricesLower();
                    Kee.AddSubmatrixSymmetric(kee, enrichedMap);
                    Kes.AddSubmatrix(kes, enrichedMap, standardMap);
                    Kec.AddSubmatrix(kes, enrichedMap, constrainedMap);
                }
            }

            return (Kee, Kes, Kec);
        }

        /// <summary>
        /// Create the signed boolean matrix with columns corresponding to all dofs in the system, to ensure continuity of 
        /// displacements. 
        /// </summary>
        /// <param name="cluster"></param>
        /// <returns></returns>
        public SignedBooleanMatrix BuildGlobalSignedBooleanMatrix(XCluster2D cluster)
        {
            Dictionary<XNode2D, SortedSet<XSubdomain2D>> nodeMembership = cluster.FindEnrichedBoundaryNodeMembership();
            int numContinuityEquations = CountContinuityEquations(cluster, nodeMembership);

            //TODO: numCols could have been computed during matrix subdomain matrix assembly.
            int numColumns = cluster.DofOrderer.NumStandardDofs + cluster.DofOrderer.NumEnrichedDofs; // Not sure about the std dofs
            var booleanMatrix = new SignedBooleanMatrix(numContinuityEquations, numColumns);

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
                            booleanMatrix.AddEntry(row, posCol, true);
                            booleanMatrix.AddEntry(row, negCol, false);
                        }
                        globalEquation += numEquationsPerDof;
                    }
                }
            }

            return booleanMatrix;
        }

        /// <summary>
        /// Create the signed boolean matrix of each subdomain to ensure continuity of displacements. If there are no boundary
        /// enriched dofs, the returned dictionary will be empty.
        /// TODO: perhaps it is more efficient to process each subdomain separately.
        /// </summary>
        /// <param name="cluster"></param>
        /// <returns></returns>
        public Dictionary<XSubdomain2D, SignedBooleanMatrix> BuildSubdomainSignedBooleanMatrices(XCluster2D cluster)
        {
            Dictionary<XNode2D, SortedSet<XSubdomain2D>> nodeMembership = cluster.FindEnrichedBoundaryNodeMembership();
            if (nodeMembership.Count == 0) return new Dictionary<XSubdomain2D, SignedBooleanMatrix>();
            var enrichedSubdomains = new HashSet<XSubdomain2D>();
            foreach (var nodeSubdomains in nodeMembership) enrichedSubdomains.UnionWith(nodeSubdomains.Value);

            ContinuityEquations eqs = FindContinuityEquations(nodeMembership);
            Console.WriteLine("Number of continuity equations = " + eqs.NumEquations);
            if (eqs.NumEquations == 0) return new Dictionary<XSubdomain2D, SignedBooleanMatrix>(); // Not sure if possible

            //TODO: numCols could have been computed during matrix subdomain matrix assembly.
            var booleanMatrices = new Dictionary<XSubdomain2D, SignedBooleanMatrix>();
            foreach (var subdomain in enrichedSubdomains)
            {
                booleanMatrices.Add(subdomain, 
                    new SignedBooleanMatrix(eqs.NumEquations, subdomain.DofOrderer.NumEnrichedDofs));
            }

            for (int i = 0; i < eqs.NumEquations; ++i)
            {
                XNode2D node = eqs.Nodes[i];
                EnrichedDof dof = eqs.Dofs[i];

                XSubdomain2D subdomainPos = eqs.PositiveSubdomains[i];
                // Presumably this will not try to access singular boundary Heaviside dofs.
                int dofIdxPos = subdomainPos.DofOrderer.GetSubdomainEnrichedDofOf(node, dof);
                booleanMatrices[subdomainPos].AddEntry(i, dofIdxPos, true);

                XSubdomain2D subdomainNeg = eqs.NegativeSubdomains[i];
                int dofIdxNeg = subdomainNeg.DofOrderer.GetSubdomainEnrichedDofOf(node, dof);
                booleanMatrices[subdomainNeg].AddEntry(i, dofIdxNeg, false);

            }

            return booleanMatrices;
        }

        /// <summary>
        /// Create the signed boolean matrix of each subdomain to ensure continuity of displacements.
        /// TODO: perhaps it is more efficient to process each subdomain separately.
        /// </summary>
        /// <param name="cluster"></param>
        /// <returns></returns>
        public Dictionary<XSubdomain2D, SignedBooleanMatrix> BuildSubdomainSignedBooleanMatricesOLD(XCluster2D cluster)
        {
            Dictionary<XNode2D, SortedSet<XSubdomain2D>> nodeMembership = cluster.FindEnrichedBoundaryNodeMembership();
            var enrichedSubdomains = new HashSet<XSubdomain2D>();
            foreach (var nodeSubdomains in nodeMembership) enrichedSubdomains.UnionWith(nodeSubdomains.Value);
            int numContinuityEquations = CountContinuityEquations(cluster, nodeMembership);

            //TODO: numCols could have been computed during matrix subdomain matrix assembly.
            var booleanMatrices = new Dictionary<XSubdomain2D, SignedBooleanMatrix>();
            foreach (var subdomain in enrichedSubdomains)
            {
                booleanMatrices.Add(subdomain,
                    new SignedBooleanMatrix(numContinuityEquations, subdomain.DofOrderer.NumEnrichedDofs));
            }

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

                            XSubdomain2D posSubdomain = positiveSubdomains[equation];
                            int posCol = posSubdomain.DofOrderer.GetSubdomainEnrichedDofOf(node, dof);
                            booleanMatrices[posSubdomain].AddEntry(row, posCol, true);

                            XSubdomain2D negSubdomain = negativeSubdomains[equation];
                            int negCol = negSubdomain.DofOrderer.GetSubdomainEnrichedDofOf(node, dof);
                            booleanMatrices[negSubdomain].AddEntry(row, negCol, false);
                        }
                        globalEquation += numEquationsPerDof;
                    }
                }
            }

            return booleanMatrices;
        }

        /// <summary>
        /// Find the number of continuity equations = number of rows in the signed boolean matrices. 
        /// TODO: It could be computed while each row is processed
        /// </summary>
        /// <param name="cluster"></param>
        /// <returns></returns>
        public int CountContinuityEquations(XCluster2D cluster, Dictionary<XNode2D, SortedSet<XSubdomain2D>> nodeMembership)
        {
            int numContinuityEquations = 0;

            bool singularHeavisideDofsExists = true;
            //bool singularHeavisideDofsExists = false; //TODO: comment out when done debugging
            //foreach (var subdomain in cluster.Subdomains)
            //{
            //    if (subdomain.DofOrderer.SingularHeavisideNodes.Count > 0)
            //    {
            //        singularHeavisideDofsExists = true;
            //        break;
            //    }
            //}
            
            if (singularHeavisideDofsExists)
            {
                // The following works even if a node is not enriched with Heaviside in one subdomain, but is enriched in others
                foreach (var nodeSubdomains in nodeMembership)
                {
                    XNode2D node = nodeSubdomains.Key;
                    foreach (EnrichedDof dof in node.EnrichedDofs) //TODO: perhaps this should be done for enrichment items instead to reduce cost
                    {
                        int dofMultiplicity = 0;
                        foreach (XSubdomain2D subdomain in nodeSubdomains.Value)
                        {
                            if (subdomain.DofOrderer.BoundaryDofs[node].Contains(dof)) ++dofMultiplicity;
                        }
                        numContinuityEquations += dofMultiplicity - 1;
                    }
                }
            }
            else
            {
                // The following works if dofs of the nodes are enriched in all subdomains
                foreach (var nodeSubdomains in nodeMembership)
                {
                    int nodeMultiplicity = nodeSubdomains.Value.Count;
                    int numEnrichedDofs = nodeSubdomains.Key.EnrichedDofsCount;
                    numContinuityEquations += (nodeMultiplicity - 1) * numEnrichedDofs;
                }
            }

            //Console.WriteLine("Num continuity equations = " + numContinuityEquations);

            return numContinuityEquations;
        }

        public ContinuityEquations FindContinuityEquations(Dictionary<XNode2D, SortedSet<XSubdomain2D>> nodeMembership)
        {
            var eqsBuilder = new ContinuityEquations.Builder();
            foreach (var nodeSubdomains in nodeMembership)
            {
                XNode2D node = nodeSubdomains.Key;
                foreach (EnrichedDof dof in node.EnrichedDofs)
                {
                    foreach (XSubdomain2D subdomain in nodeSubdomains.Value)
                    {
                        if (subdomain.DofOrderer.BoundaryDofs[node].Contains(dof)) eqsBuilder.Register(node, dof, subdomain);
                    }
                }
            }
            return eqsBuilder.Build(ContinuityEquations.EquationOrder.DofMajor); // TODO: Subdomain major is probably better for the QR factorization
        }
    }
}
