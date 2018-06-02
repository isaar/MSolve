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

            // Matrix that will undergo QR:
            // (B*Pe)^T = [B1*inv(U1) B2*inv(U2) ...]^T = [inv(L1)*B1^T inv(L2)*B2^T ...]
            // Dims: B = numEquations -by- numDofsEnr, Pe = numDofsEnr -by- numDofsEnr, (B*Pe)^T = numDofsEnr -by- numEquations
            var BPeTransp = Matrix.CreateZero(dimensions.NumDofsEnr, dimensions.NumEquations);

            foreach (var subB in allB)
            {
                var subdomain = subB.Key;
                SignedBooleanMatrix B = subB.Value;
                CholeskySuiteSparse Pe = allPe[subdomain];

                // Dimensions and offsets
                int numDofsSubdomain = B.NumColumns;
                int numDofsBoundary = 0;
                foreach (var nodeDofsPair in subdomain.DofOrderer.BoundaryDofs) numDofsBoundary += nodeDofsPair.Value.Count;
                int boundaryDofsStart = numDofsSubdomain - numDofsBoundary;
                int globalStart = dimensions.SubdomainStarts[subdomain] - dimensions.NumDofsStd;

                // Contribution to the matrix that will undergo QR
                foreach (var row in B.CopyNonZeroRowsToVectors())
                {
                    int rowIdx = row.Key;
                    Vector rowVector = row.Value;
                    Vector contribution = allPe[subdomain].ForwardSubstitution(rowVector); //TODO: only apply forward substitution for the last dofs
                    BPeTransp.SetColumn(rowIdx, globalStart + boundaryDofsStart,
                        contribution.Slice(boundaryDofsStart, numDofsSubdomain)); //TODO: directly copy them. No need to create a temporary Vector first
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
            Matrix Q = qr.GetEconomyFactorQ();
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

            return new QR(Q, R);
        }

        public CholeskySuiteSparse CreateEnrichedPreconditioner(DOKSymmetricColMajor Kee)
        {
            // Enriched preconditioner = cholesky factor U
            var (valuesEnr, rowIndicesEnr, colOffsetsEnr) = Kee.BuildSymmetricCSCArrays(true);
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
                Vector row = Be.CopyRowToVector(i);
                Vector solution = Pe.ForwardSubstitution(row);
                if (nonZerosRows.Contains(i)) 
                {
                    // Check that only the last columns, which correspond to boundary dofs are non zero.
                    Vector internalColsB = row.Slice(0, boundaryDofsStart);
                    Vector boundaryColsB = row.Slice(boundaryDofsStart, numDofsSubdomain);
                    Vector internalRowsSolution = solution.Slice(0, boundaryDofsStart);
                    Vector boundaryRowsSolution = solution.Slice(boundaryDofsStart, numDofsSubdomain);

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

        private class QR : IFactorizationLQ
        {
            private readonly Matrix Q1;
            private readonly TriangularUpper R1;

            internal QR(Matrix Q, TriangularUpper R)
            {
                this.Q1 = Q;
                this.R1 = R;
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
                return Q1.MultiplyRight(x, !transposeQ); //The stored Q1 of QR is the transpose of the Q1 of LQ
            }
        }
    }
}
