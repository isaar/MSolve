using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    internal class ProjectorMatrix
    {
        private readonly Dictionary<int, SignedBooleanMatrix> booleanMatrices;
        private readonly Dictionary<int, List<Vector>> rigidBodyModes;
        private readonly Matrix matrixQ;
        private Matrix matrixG;
        private CholeskyFull factorGQG; // Because Karmath suggests using POTRF, POTRS.

        internal ProjectorMatrix(Dictionary<int, SignedBooleanMatrix> booleanMatrices, 
            Dictionary<int, List<Vector>> rigidBodyModes, Matrix matrixQ)
        {
            this.booleanMatrices = booleanMatrices;
            this.rigidBodyModes = rigidBodyModes;
            this.matrixQ = matrixQ;
        }

        /// <summary>
        /// inv(G^T * Q * G)
        /// </summary>
        internal void InvertCoarseProblemMatrix()
        {
            CalculateMatrixG();
            Matrix GQG = matrixG.ThisTransposeTimesOtherTimesThis(matrixQ);
            factorGQG = GQG.FactorCholesky(true);
        }

        /// <summary>
        /// a = inv(G^T*Q*G) * G^T * Q * (Fe * λ - d)
        /// </summary>
        /// <param name="flexibilityTimeslagrangeMultipliers"></param>
        /// <param name="boundaryDisplacements"></param>
        /// <returns></returns>
        internal Vector CalculateRigidBodyMotionsCoefficients(Vector flexibilityTimeslagrangeMultipliers, 
            Vector boundaryDisplacements)
        {
            Vector x = flexibilityTimeslagrangeMultipliers - boundaryDisplacements;
            return factorGQG.SolveLinearSystem(matrixG.Multiply(matrixQ * x, true));
        }

        /// <summary>
        /// λ0 = Q * G * inv(G^T * Q * G) * e
        /// </summary>
        internal void InitializeLagrangeMultipliers(Vector rigidBodyMotionsWork, Vector lagrange)
        {
            matrixQ.MultiplyIntoResult(matrixG * (factorGQG.SolveLinearSystem(rigidBodyMotionsWork)), lagrange);
        }

        /// <summary>
        ///  P = I - Q * G * inv(G^T*Q*G) * G^T,  projected = P * original
        /// </summary>
        internal void ProjectVector(Vector original, Vector projected)
        {
            // P * x = x - Q * (G * (inv(G^T*Q*G) * (G^T * x))
            projected.CopyFrom(original);
            projected.SubtractIntoThis(matrixQ * (matrixG * (factorGQG.SolveLinearSystem(matrixG.Multiply(original, true)))));
        }

        /// <summary>
        /// G = [B(1)*R(1) B(2)*R(2) ... B(ns)*R(ns)
        /// </summary>
        private void CalculateMatrixG()
        {
            int numEquations = booleanMatrices.First().Value.NumRows;
            int numRbms = 0;
            foreach (int subdomain in rigidBodyModes.Keys) numRbms += rigidBodyModes[subdomain].Count;

            matrixG = Matrix.CreateZero(numEquations, numRbms);
            int colCounter = 0;
            foreach (int subdomain in booleanMatrices.Keys)
            {
                SignedBooleanMatrix matrixB = booleanMatrices[subdomain];
                List<Vector> matrixR = rigidBodyModes[subdomain];
                foreach (Vector columnR in matrixR)
                {
                    Vector columnG = matrixB.MultiplyRight(columnR, false);
                    matrixG.SetSubcolumn(colCounter++, columnG);
                }
            }
        }
    }
}
