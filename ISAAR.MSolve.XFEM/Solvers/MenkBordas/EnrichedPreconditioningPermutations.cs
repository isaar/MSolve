using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class EnrichedPreconditioningPermutations : IEnrichedPreconditioning
    {
        public EnrichedPreconditioningPermutations()
        {
            this.Ordering = new EnrichedOrderingCamd();
        }

        public IEnrichedOrdering Ordering { get; }

        public IFactorizationLQ CreateContinuityEquationsPreconditioner(MenkBordasSystem.Dimensions dimensions,
            IReadOnlyDictionary<XSubdomain2D, SignedBooleanMatrix> allB,
            IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> allPe)
        {
            if (allB.Count < 2) throw new ArgumentException("There must be at least 2 subdomains.");
            IEnumerable<XSubdomain2D> subdomains = allB.Keys; // WARNING: Use this collection in all subsequent loops


            //Go through the boundary dofs of each subdomain, store it original index in Q, find its new index, so that it is 
            // moved to the top, after all other boundary dofs of all subdomains, and store it in a permutation matrix 


            // For each subdomain find the number of boundary dofs, the offset of the first into the subdomain's B, into the 
            // original Q and into the permuted Q
            int totalEnrichedDofs = dimensions.NumDofsEnr; // Subdomains with no boundary dofs are also included in here
            int totalBoundaryDofs = 0;
            var numBoundaryDofs = new Dictionary<XSubdomain2D, int>();
            var numSubdomainDofs = new Dictionary<XSubdomain2D, int>();
            var offsetB = new Dictionary<XSubdomain2D, int>();
            var offsetQoriginal = new Dictionary<XSubdomain2D, int>();
            var offsetQpermuted = new Dictionary<XSubdomain2D, int>();
            foreach (XSubdomain2D sub in subdomains)
            {
                int counter = 0;
                foreach (var nodeDofsPair in sub.DofOrderer.BoundaryDofs) counter += nodeDofsPair.Value.Count;
                numBoundaryDofs[sub] = counter;
                numSubdomainDofs[sub] = allB[sub].NumColumns;

                offsetB[sub] = numSubdomainDofs[sub] - numBoundaryDofs[sub]; // the boundary dofs are last
                offsetQoriginal[sub] = (dimensions.SubdomainStarts[sub] - dimensions.NumDofsStd) + offsetB[sub];
                offsetQpermuted[sub] = totalBoundaryDofs; // Bring all boundary dofs to the beginning

                totalBoundaryDofs += numBoundaryDofs[sub]; // this must be done after offsetQpermuted
            }

            // Build the permutation matrix
            var permutation = PermutationMatrix.CreateIdentity(totalEnrichedDofs);
            int nextNonzeroRowQ = 0;
            foreach (XSubdomain2D sub in subdomains)
            {
                for (int i = 0; i < numBoundaryDofs[sub]; ++i)
                {
                    // Move each boundary dof under the rest boundary dofs
                    int originalIdx = offsetQoriginal[sub] + i;
                    int permutedIdx = nextNonzeroRowQ++;
                    permutation.ExchangeRows(originalIdx, permutedIdx);
                }
            }

            // Matrix that will undergo QR:
            // (B*Pe)^T = [B1*inv(U1) B2*inv(U2) ...]^T = [inv(L1)*B1^T inv(L2)*B2^T ...]
            // Dims: B = numEquations -by- numDofsEnr, Pe = numDofsEnr -by- numDofsEnr, (B*Pe)^T = numDofsEnr -by- numEquations
            // All non zero rows correspong to boundary dofs and are moved to the top. The rest of (B*Pe)^T is 0 and will remain 
            // so after factorization. Thus only totalBoundaryDofs -by -numEquations is needed
            var BPeTransp = Matrix.CreateZero(totalBoundaryDofs, dimensions.NumEquations);

            foreach (XSubdomain2D sub in subdomains)
            {
                SignedBooleanMatrix B = allB[sub];
                CholeskySuiteSparse Pe = allPe[sub];

                // Dimensions and offsets
                //int numDofsSubdomain = B.NumColumns;
                //int numDofsBoundary = 0;
                //foreach (var nodeDofsPair in subdomain.DofOrderer.BoundaryDofs) numDofsBoundary += nodeDofsPair.Value.Count;
                //int boundaryDofsStart = numDofsSubdomain - numDofsBoundary;
                //int globalStart = dimensions.SubdomainStarts[subdomain] - dimensions.NumDofsStd;

                // Contribution to the matrix that will undergo QR
                foreach (var row in B.CopyNonZeroRowsToVectors())
                {
                    int rowIdx = row.Key;
                    Vector rowVector = row.Value;
                    Vector contribution = allPe[sub].ForwardSubstitution(rowVector); //TODO: only apply forward substitution for the last dofs
                    //BPeTransp.SetColumn(rowIdx, offsetQoriginal[sub],
                    BPeTransp.SetSubcolumn(rowIdx, contribution.GetSubvector(offsetB[sub], numSubdomainDofs[sub]), 
                        offsetQpermuted[sub]); //TODO: directly copy them. No need to create a temporary Vector first
                    //Debug.Assert(contribution.Equals(contribution2, new LinearAlgebra.Testing.Utilities.ValueComparer(double.Epsilon)));
                }

                #region debug
                //CheckBPSparsity(subdomain, subB.Value, allPe[subdomain]);
                //TODO: Solving each column separately gives slightly different results, than the following. This is a SuiteSparse thing. Investigate. 
                //Matrix contributionOLD = allPe[subdomain].ForwardSubstitution(subB.Value.CopyToFullMatrix(true)); 
                //Matrix part = contributionOLD.Slice(boundaryDofsStart, numDofsSubdomain, 0, dimensions.NumEquations);
                //BPeTransp.SetSubmatrix(globalStart + boundaryDofsStart, 0, part);
                #endregion

                //TODO: most columns correspond to continuity equations between dofs of other subdomains. Do not solve for them, since they are 0. 
            }

            // LQ factorization 
            //TODO: various optimizations might be possible here
            var qr = BPeTransp.FactorQR();
            Matrix blockQ1 = qr.GetEconomyFactorQ();
            TriangularUpper R = qr.GetEconomyFactorR();

            #region Debug
            //FullMatrixWriter.NumericFormat = new ExponentialFormat { NumDecimalDigits = 4 };
            //Console.WriteLine("Before QR: Pe^T * B^T = ");
            //(new FullMatrixWriter(BPeTransp)).WriteToConsole();
            //Console.WriteLine();

            //Console.WriteLine("Q = ");
            //(new FullMatrixWriter(qr.GetFactorQ().Slice(0, numDofsEnr, 0, numEquations))).WriteToConsole();
            //Console.WriteLine();

            //Console.WriteLine("R = ");
            //(new FullMatrixWriter(qr.GetFactorR().Slice(0, numEquations, 0, numEquations))).WriteToConsole();
            //Console.WriteLine();
            #endregion

            return new PQR(blockQ1, R, permutation);
        }

        public CholeskySuiteSparse CreateEnrichedPreconditioner(DokSymmetric Kee)
        {
            // Enriched preconditioner = cholesky factor U
            var (valuesEnr, rowIndicesEnr, colOffsetsEnr) = Kee.BuildSymmetricCscArrays(true);
            return CholeskySuiteSparse.Factorize(Kee.NumRows, valuesEnr.Length, valuesEnr, rowIndicesEnr,
                colOffsetsEnr, true, SuiteSparseOrdering.Natural);
            //TODO: perhaps I should discard Kee here, instead of in MenkBordasSystem.
        }

        private static void CheckBPSparsity(XSubdomain2D subdomain, SignedBooleanMatrix Be, CholeskySuiteSparse Pe)
        {
            int numDofsSubdomain = Be.NumColumns;
            int numDofsBoundary = 0;
            foreach (var nodeDofsPair in subdomain.DofOrderer.BoundaryDofs) numDofsBoundary += nodeDofsPair.Value.Count;
            int boundaryDofsStart = numDofsSubdomain - numDofsBoundary;
            HashSet<int> nonZerosRows = new HashSet<int>(Be.FindNonZeroRows());

            for (int i = 0; i < Be.NumRows; ++i)
            {
                Vector row = Be.GetRow(i);
                Vector solution = Pe.ForwardSubstitution(row);
                if (nonZerosRows.Contains(i)) 
                {
                    // Check that only the last columns, which correspond to boundary dofs are non zero.
                    Vector internalColsB = row.GetSubvector(0, boundaryDofsStart);
                    Vector boundaryColsB = row.GetSubvector(boundaryDofsStart, numDofsSubdomain);
                    Vector internalRowsSolution = solution.GetSubvector(0, boundaryDofsStart);
                    Vector boundaryRowsSolution = solution.GetSubvector(boundaryDofsStart, numDofsSubdomain);

                    Debug.Assert(internalColsB.IsZero(0.0), 
                        $"In row {i} of Be the columns corresponding to internal dofs are not 0 as they must be");
                    Debug.Assert(!boundaryColsB.IsZero(0.0),
                        $"In row {i} of Be the columns corresponding to boundary dofs are all 0. At least one must be non zero");
                    Debug.Assert(internalRowsSolution.IsZero(0.0),
                        $"In col {i} of inv(Pe^T) * Be^T the rows corresponding to internal dofs are not 0 as they must be");
                    Debug.Assert(!boundaryRowsSolution.IsZero(0.0),
                        $"In row {i} of inv(Pe^T) * Be^T the rows corresponding to boundary dofs are all 0. At least one must be non zero");

                }
                else  // TODO: find these continuity equations instead of using Be.FindNonZeroRows()
                {
                    // Check that rows corresponding to continuity equations between dofs of other subdomains are zero
                    Debug.Assert(row.IsZero(0.0),
                        $"Rows of Be corresponding to continuity dofs of other subdomains must be 0, but row {i} isn't");
                    Debug.Assert(solution.IsZero(0.0),
                        $"Cols of inv(Pe^T) * Be^T corresponding to continuity dofs of other subdomains must be 0, but col {i} isn't");
        }
            }

        }

        private class PQR : IFactorizationLQ
        {
            private readonly PermutationMatrix permutation;
            private readonly Matrix blockQ1;
            private readonly TriangularUpper R1;

            internal PQR(Matrix blockQ1, TriangularUpper R, PermutationMatrix permutation)
            {
                this.blockQ1 = blockQ1;
                this.R1 = R;
                this.permutation = permutation;
            }

            public Vector InverseLTimesVector(Vector x, bool transposePreconditioner)
            {
                if (transposePreconditioner) // (L^(-T))^T * x = R^T \ x
                {
                    return R1.SolveLinearSystem(x, true);
                }
                else // L^(-T) * x = R \ x
                {
                    return R1.SolveLinearSystem(x, false);
                }
            }

            public Vector QTimesVector(Vector x, bool transposeQ)
            {
                // TODO: in a lot of steps, I could avoid zero
                // Since I work with QR instead of LQ, the transposed flags are opposite
                if (transposeQ) 
                {
                    // Q1 * x = [blockQ1 * x ; 0] = upper block vector
                    //TODO: The following could be done by multiply matrix * subVector1 into subVector2
                    Vector blockQ1x = blockQ1.MultiplyRight(x, false);
                    var Q1x = Vector.CreateZero(permutation.NumRows);
                    Q1x.SetSubvector(blockQ1x, 0);

                    // P^T * (Q1*x)
                    return permutation.MultiplyRight(Q1x, true);
                }
                else
                {
                    // P * x
                    Vector Px = permutation.MultiplyRight(x, false);

                    // (P^T * Q1)^T = Q1^T * P*x = Q1^T * block(Px) = dense vector
                    Vector blockPx = Px.GetSubvector(0, blockQ1.NumRows);
                    return blockQ1.MultiplyRight(blockPx, true);
                }
            }
        }
    }
}
