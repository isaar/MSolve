using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class EnrichedPreconditioningNaive : IEnrichedPreconditioning
    {
        public EnrichedPreconditioningNaive()
        {
            this.Ordering = new EnrichedOrderingNatural();
        }

        public IEnrichedOrdering Ordering { get; } 

        public IFactorizationLQ CreateContinuityEquationsPreconditioner(MenkBordasSystem.Dimensions dimensions, 
            IReadOnlyDictionary<XSubdomain2D, SignedBooleanMatrix> B, 
            IReadOnlyDictionary<XSubdomain2D, CholeskySuiteSparse> Pe)
        {
            if (B.Count < 2) throw new ArgumentException("There must be at least 2 subdomains.");

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

            return new LQ(L, Q);
        }

        public CholeskySuiteSparse CreateEnrichedPreconditioner(DokSymmetric Kee)
        {
            // Enriched preconditioner = cholesky factor U
            var (valuesEnr, rowIndicesEnr, colOffsetsEnr) = Kee.BuildSymmetricCscArrays(true);
            return CholeskySuiteSparse.Factorize(Kee.NumRows, valuesEnr.Length, valuesEnr, rowIndicesEnr,
                colOffsetsEnr, true);
            //TODO: perhaps I should discard Kee here, instead of in MenkBordasSystem.
        }

        private class LQ : IFactorizationLQ
        {
            private readonly Matrix L;
            private readonly Matrix Q;

            internal LQ(Matrix L, Matrix Q)
            {
                this.L = L;
                this.Q = Q;
            }

            public Vector InverseLTimesVector(Vector x, bool transposePreconditioner)
            {
                // The preconditioner has L^(-T) = R^(-1), thus the transpositions are:
                if (transposePreconditioner) 
                {
                    return L.Invert() * x;
                }
                else
                {
                    return L.Transpose().Invert() * x;
                }
            }

            public Vector QTimesVector(Vector x, bool transposeQ)
            {
                return Q.Multiply(x, transposeQ);
            }
        }
    }
}
