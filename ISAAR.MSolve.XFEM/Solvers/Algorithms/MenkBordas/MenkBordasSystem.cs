using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.Algorithms.MenkBordas
{
    class MenkBordasSystem
    {
        private readonly int numSubdomains;
        private readonly int numContinuityEquations;

        public MenkBordasSystem(int numSubdomains, int numContinuityEquations)
        {
            if (numSubdomains < 1) throw new ArgumentException(
                $"There must be at least 1 subdomain, but {numSubdomains} were declared");
            if (numContinuityEquations < 0) throw new ArgumentException(
               $"There must be at least 0 continuity equations, but {numSubdomains} were declared");

            this.numSubdomains = numSubdomains;
            this.numContinuityEquations = numContinuityEquations;

        }

        // Matrices
        public IMatrixView Kss;
        public readonly List<IMatrixView> Kee;
        public readonly List<IMatrixView> Kes;
        public readonly List<IMatrixView> Kse;
        public readonly List<SignedBooleanMatrix> B;

        // Right hand side vectors
        public Vector bs;
        public readonly List<Vector> be;

        // Initial guesses (optional) that will be overwritten by solutions
        public Vector xs;
        public readonly List<Vector> xe;

        public void CheckDimensions()
        {
            // Standard dofs
            Preconditions.CheckSystemSolutionDimensions(Kss, bs);

            // Number of subdomains
            if ((Kee.Count != numSubdomains) || (Kes.Count != numSubdomains) || (Kse.Count != numSubdomains) 
                || (B.Count != numSubdomains) || (be.Count != numSubdomains))
            {
                throw new ArgumentException($"Mismatch: {numSubdomains} subdomains were declared, but there are"
                    + $" {Kee.Count} enriched-enriched stiffness matrices, {Kes.Count} enriched-standard stiffness matrices,"
                    + $" {Kse.Count} standard-enriched stiffness matrices, {B.Count} signed boolean matrices and"
                    + $" {be.Count} right hand side vectors.");
            }

            // Enriched dofs
            for (int i = 0; i < numSubdomains; ++i)
            {
                //TODO: All dimensions could be checked.
                Preconditions.CheckSystemSolutionDimensions(Kee[i], be[i]);
                Preconditions.CheckSystemSolutionDimensions(Kes[i], be[i]);
                Preconditions.CheckSystemSolutionDimensions(Kse[i], bs);
                Preconditions.CheckMultiplicationDimensions(B[i].NumColumns, be[i].Length); //rhs.Length = lhs.Length
            }
        }

        public MenkBordasMatrix BuildMatrix()
        {
            return new MenkBordasMatrix(numSubdomains, numContinuityEquations, Kss, Kee.ToArray(), Kes.ToArray(), Kse.ToArray(),
                B.ToArray());
        }

        public MenkBordasVector BuildRhsVector()
        {
            return new MenkBordasVector(numSubdomains, numContinuityEquations, bs, be.ToArray(), 
                Vector.CreateZero(numContinuityEquations)); //TODO: avoid creating the 0 vector, after validating it will stay 0
        }
    }
}
