using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

//TODO: Allow it to update some submatrices, while keeping others the same. 
//TODO: Break up Q and L into submatrices or store the submatrices that assemble the intermediate matrix before QR factorization.
//TODO: I need a new DOK class:
//      1) The same DOK should produce the full CSR (for multiplications) and the upper CSC (for factorization).
//          Fortunately the arrays for the upper CSC are identical to the ones for lower CSR.
//      2) Sorting the DOK should only be done once. It actually costs substantially. Also a lot of other operations can only 
//          once, when creating the full and lower CSR.
//      3) Ok having the full CSR/CSC matrix is good for multiplications. However can't I factorize using just that? Do I need
//          a separate upper copy? I think SuiteSparse can ignore the lower/upper entries
//TODO: Q and L are dense now. Try to take advantage of any sparsity. Perhaps create some by appropriate DD or permutations.
//      Also try to take advantage of the fact that B are boolean matrices, with at most 2 entries per column
// TODO: Delete this class after resolving the other TODOs or transfering the comments to the new OOP design.
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    //TODO: apply AMD before factorizing. Perhaps this should be done by the solver itself, so that the dofOrders are also updated.
    static class MenkBordasPreconditioner
    {
        
        public static CholeskySuiteSparse CreateStandardPreconditioner(DokSymmetric Kss)
        {
            // Standard preconditioner = cholesky factor U
            var (valuesStd, rowIndicesStd, colOffsetsStd) = Kss.BuildSymmetricCscArrays(true);
            return CholeskySuiteSparse.Factorize(Kss.NumRows, valuesStd.Length, valuesStd, rowIndicesStd, colOffsetsStd,
                true, SuiteSparseOrdering.Natural);
        }

        public static CholeskySuiteSparse[] CreateEnrichedPreconditioners(MenkBordasSystem.Dimensions dimensions,
             IReadOnlyList<DokSymmetric> Kee)
        {
            var Pe = new CholeskySuiteSparse[dimensions.NumSubdomains];
            for (int i = 0; i < dimensions.NumSubdomains; ++i)
            {
                // Enriched preconditioner = cholesky factor U
                var (valuesEnr, rowIndicesEnr, colOffsetsEnr) = Kee[i].BuildSymmetricCscArrays(true);
                Pe[i] = CholeskySuiteSparse.Factorize(Kee[i].NumRows, valuesEnr.Length, valuesEnr, rowIndicesEnr, colOffsetsEnr,
                    true, SuiteSparseOrdering.Natural);
            }
            return Pe;
        }

        public static CholeskySuiteSparse CreateEnrichedPreconditioner(DokSymmetric Kee)
        {
            // Enriched preconditioner = cholesky factor U
            var (valuesEnr, rowIndicesEnr, colOffsetsEnr) = Kee.BuildSymmetricCscArrays(true);
            return CholeskySuiteSparse.Factorize(Kee.NumRows, valuesEnr.Length, valuesEnr, rowIndicesEnr, 
                colOffsetsEnr, true, SuiteSparseOrdering.Natural);
        }

        public static (Matrix L, Matrix Q) CreateContinuityEquationsPreconditioners(MenkBordasSystem.Dimensions dimensions,
             IReadOnlyDictionary<XSubdomain2D, SignedBooleanMatrix> B, IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> Pe)
        {
            if (B.Count < 2 ) throw new ArgumentException("There must be at least 2 subdomains.");

            // Matrix that will undergo QR:
            // (B*Pe)^T = [B1*inv(U1) B2*inv(U2) ...]^T = [inv(L1)*B1^T inv(L2)*B2^T ...]
            // Dims: B = numEquations -by- numDofsEnr, Pe = numDofsEnr -by- numDofsEnr, (B*Pe)^T = numDofsEnr -by- numEquations
            var BPeTransp = Matrix.CreateZero(dimensions.NumDofsEnr, dimensions.NumEquations);

            foreach (var sub in B.Keys)
            {
                // Contribution to the matrix that will undergo QR
                Matrix contribution = Pe[sub].ForwardSubstitutions(B[sub].CopyToFullMatrix(true));
                BPeTransp.SetSubmatrix(dimensions.SubdomainStarts[sub] - dimensions.NumDofsStd, 0, contribution);
            }

            // LQ factorization 
            //TODO: various optimizations might be possible here
            var qr = BPeTransp.FactorQR();
            Matrix L = qr.GetFactorR().GetSubmatrix(0, dimensions.NumEquations, 0, dimensions.NumEquations).Transpose(); //TODO: should probably use a packed UpperTriangular R
            Matrix Q = qr.GetFactorQ().GetSubmatrix(0, dimensions.NumDofsEnr, 0, dimensions.NumEquations).Transpose(); //TODO: MKL has routines that only build some columns of Q!!!

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

            return (L, Q);
        }
    }
}
