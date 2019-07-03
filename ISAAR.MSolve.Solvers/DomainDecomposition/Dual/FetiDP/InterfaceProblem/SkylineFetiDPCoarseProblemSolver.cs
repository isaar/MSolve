using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Use Skyline assembler instead of reimplementing it.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public class SkylineFetiDPCoarseProblemSolver : IFetiDPCoarseProblemSolver
    {
        private readonly IReorderingAlgorithm reordering;
        private readonly IReadOnlyList<ISubdomain> subdomains;
        private LdlSkyline inverseGlobalKccStar;

        public SkylineFetiDPCoarseProblemSolver(IReadOnlyList<ISubdomain> subdomains, IReorderingAlgorithm reordering)
        {
            this.subdomains = subdomains;
            this.reordering = reordering;
        }

        public void ClearCoarseProblemMatrix()
        {
            inverseGlobalKccStar = null;
        }

        public void CreateAndInvertCoarseProblemMatrix(Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains, 
            FetiDPDofSeparator dofSeparator, Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers)
        {
            SkylineMatrix globalKccStar = CreateGlobalKccStar(cornerNodesOfSubdomains, dofSeparator, matrixManagers);
            this.inverseGlobalKccStar = globalKccStar.FactorLdl(true);
        }

        public Vector CreateCoarseProblemRhs(FetiDPDofSeparator dofSeparator,
            Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers,
            Dictionary<int, Vector> fr, Dictionary<int, Vector> fbc)
        {
            // Static condensation for the force vectors
            var globalFcStar = Vector.CreateZero(dofSeparator.NumGlobalCornerDofs);
            for (int s = 0; s < subdomains.Count; ++s)
            {
                IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];

                // fcStar[s] = fbc[s] - Krc[s]^T * inv(Krr[s]) * fr[s]
                // globalFcStar = sum_over_s(Lc[s]^T * fcStar[s])
                UnsignedBooleanMatrix Lc = dofSeparator.CornerBooleanMatrices[s];
                Vector temp = matrices.MultiplyInverseKrrTimes(fr[s]);
                temp = matrices.MultiplyKcrTimes(temp);
                Vector fcStar = fbc[s] - temp;
                globalFcStar.AddIntoThis(Lc.Multiply(fcStar, true));
            }
            return globalFcStar;
        }

        public Vector MultiplyInverseCoarseProblemMatrixTimes(Vector vector) => inverseGlobalKccStar.SolveLinearSystem(vector);

        public void ReorderCornerDofs(FetiDPDofSeparator dofSeparator)
        {
            if (reordering == null) return; // Use the natural ordering and do not modify any stored dof data
            var pattern = SparsityPatternSymmetric.CreateEmpty(dofSeparator.NumGlobalCornerDofs);
            for (int s = 0; s < subdomains.Count; ++s) 
            {
                // Treat each subdomain as a superelement with only its corner nodes.
                var localCornerDofOrdering = dofSeparator.SubdomainCornerDofOrderings[s];
                int numLocalCornerDofs = localCornerDofOrdering.EntryCount;
                var subdomainToGlobalDofs = new int[numLocalCornerDofs];
                foreach ((INode node, IDofType dofType, int localIdx) in localCornerDofOrdering)
                {
                    int globalIdx = dofSeparator.GlobalCornerDofOrdering[node, dofType];
                    subdomainToGlobalDofs[localIdx] = globalIdx;
                }
                pattern.ConnectIndices(subdomainToGlobalDofs, false);
            }
            (int[] permutation, bool oldToNew) = reordering.FindPermutation(pattern);
            dofSeparator.GlobalCornerDofOrdering.Reorder(permutation, oldToNew);
            dofSeparator.GlobalCornerToFreeDofMap =
                ReorderingUtilities.ReorderKeysOfDofIndicesMap(dofSeparator.GlobalCornerToFreeDofMap, permutation, oldToNew);
        }

        //TODO: Use Skyline assembler
        private SkylineMatrix CreateGlobalKccStar(Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains, 
            FetiDPDofSeparator dofSeparator, Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers)
        {
            int[] skylineColHeights = FindSkylineColumnHeights(cornerNodesOfSubdomains, dofSeparator);
            var skylineBuilder = SkylineBuilder.Create(dofSeparator.NumGlobalCornerDofs, skylineColHeights);

            for (int s = 0; s < subdomains.Count; ++s)
            {
                IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];

                // KccStar[s] = Kcc[s] - Krc[s]^T * inv(Krr[s]) * Krc[s]
                if (subdomains[s].StiffnessModified)
                {
                    Debug.WriteLine($"{this.GetType().Name}: Calculating Schur complement of remainder dofs"
                        + " for the stiffness of subdomain {s}");
                    matrices.CalcSchurComplementOfRemainderDofs(); //TODO: At this point Kcc and Krc can be cleared. Maybe Krr too.
                }

                int[] subdomainToGlobalIndices = dofSeparator.CornerBooleanMatrices[s].GetRowsToColumnsMap();
                IMatrixView subdomainMatrix = matrices.SchurComplementOfRemainderDofs;
                skylineBuilder.AddSubmatrixSymmetric(subdomainMatrix, subdomainToGlobalIndices);
            }

            return skylineBuilder.BuildSkylineMatrix();
        }

        private int[] FindSkylineColumnHeights(Dictionary<int, HashSet<INode>> cornerNodesOfSubdomains, 
            FetiDPDofSeparator dofSeparator)
        {
            //only entries above the diagonal count towards the column height
            int[] colHeights = new int[dofSeparator.NumGlobalCornerDofs]; 
            for (int s = 0; s < subdomains.Count; ++s)
            {
                HashSet<INode> cornerNodes = cornerNodesOfSubdomains[s];

                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0.0 entries in the element stiffness matrix.
                int minDof = Int32.MaxValue;
                foreach (INode node in cornerNodes)
                {
                    foreach (int dof in dofSeparator.GlobalCornerDofOrdering.GetValuesOfRow(node)) minDof = Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (INode node in cornerNodes)
                {
                    foreach (int dof in dofSeparator.GlobalCornerDofOrdering.GetValuesOfRow(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return colHeights;
        }
    }
}
