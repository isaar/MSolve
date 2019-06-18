using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

//TODO: If Q is identity, optimizations will be possible: E.g. multiplications with Q do not have to copy the vector/matrix
//      (for lagranges they do), G^T*G might be optimizable, etc. Should there be a dedicated class for such cases.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection
{
    public class Feti1Projection : IInterfaceProjection
    {
        private readonly Dictionary<int, SignedBooleanMatrixColMajor> booleanMatrices;
        private readonly Dictionary<int, IFeti1SubdomainMatrixManager> matrixManagers;
        private readonly IMatrixQ matrixQ;
        private Matrix matrixG;
        private CholeskyFull factorGQG; // Because Karmath suggests using POTRF, POTRS.

        internal Feti1Projection(Dictionary<int, SignedBooleanMatrixColMajor> booleanMatrices,
            Dictionary<int, IFeti1SubdomainMatrixManager> matrixManagers, IMatrixQ matrixQ)
        {
            this.booleanMatrices = booleanMatrices;
            this.matrixManagers = matrixManagers;
            this.matrixQ = matrixQ;
        }

        /// <summary>
        /// λ0 = Q * G * inv(G^T * Q * G) * e
        /// </summary>
        public Vector CalcParticularLagrangeMultipliers(Vector rigidBodyModesWork)
            => matrixQ.Multiply(matrixG * (factorGQG.SolveLinearSystem(rigidBodyModesWork)));

        /// <summary>
        /// a = inv(G^T*Q*G) * G^T * Q * (Fe * λ - d)
        /// </summary>
        /// <param name="flexibilityTimeslagrangeMultipliers"></param>
        /// <param name="boundaryDisplacements"></param>
        public Vector CalcRigidBodyModesCoefficients(Vector flexibilityTimeslagrangeMultipliers, 
            Vector boundaryDisplacements)
        {
            Vector x = flexibilityTimeslagrangeMultipliers - boundaryDisplacements;
            return factorGQG.SolveLinearSystem(matrixG.Multiply(matrixQ.Multiply(x), true));
        }

        /// <summary>
        /// inv(G^T * Q * G)
        /// </summary>
        public void InvertCoarseProblemMatrix()
        {
            CalculateMatrixG();
            //Matrix GQG = matrixG.ThisTransposeTimesOtherTimesThis(matrixQ);
            Matrix GQG = matrixG.MultiplyRight(matrixQ.Multiply(matrixG), true);

            var writer = new LinearAlgebra.Output.FullMatrixWriter();
            //writer.WriteToFile(GQG, @"C:\Users\Serafeim\Desktop\output.txt");


            factorGQG = GQG.FactorCholesky(true);
        }

        /// <summary>
        ///  P = I - Q * G * inv(G^T*Q*G) * G^T,  projected = P * original
        /// </summary>
        public void ProjectVector(Vector original, Vector projected, bool transposeProjector)
        {
            if (transposeProjector)
            {
                // I, Q and inv(G^T*Q*G) are symmetric thus
                // P^T * x = I^T * x - G * inv(G^T*Q*G)^T * G^T * Q^T * x <=>
                // P^T * x = x - G * (inv(G^T*Q*G) * (G^T * (Q * x)))
                projected.CopyFrom(original);
                projected.SubtractIntoThis(matrixG * factorGQG.SolveLinearSystem(
                    matrixG.Multiply(matrixQ.Multiply(original), true)));
            }
            else
            {
                // P * x = x - Q * (G * (inv(G^T*Q*G) * (G^T * x)))
                projected.CopyFrom(original);
                projected.SubtractIntoThis(matrixQ.Multiply(matrixG * (
                    factorGQG.SolveLinearSystem(matrixG.Multiply(original, true)))));
            }
        }

        /// <summary>
        /// G = [B(1)*R(1) B(2)*R(2) ... B(ns)*R(ns)]
        /// </summary>
        private void CalculateMatrixG()
        {
            int numEquations = booleanMatrices.First().Value.NumRows;
            int numRbms = 0;
            foreach (IFeti1SubdomainMatrixManager matrixManager in matrixManagers.Values)
            {
                numRbms += matrixManager.RigidBodyModes.Count;
            }

            matrixG = Matrix.CreateZero(numEquations, numRbms);
            int colCounter = 0;
            foreach (int s in booleanMatrices.Keys)
            {
                SignedBooleanMatrixColMajor matrixB = booleanMatrices[s];
                List<Vector> matrixR = matrixManagers[s].RigidBodyModes;
                foreach (Vector columnR in matrixR)
                {
                    Vector columnG = matrixB.Multiply(columnR, false); //TODO: Perhaps I can cache some of these for reanalysis.
                    matrixG.SetSubcolumn(colCounter++, columnG);
                }
            }
        }
    }
}
