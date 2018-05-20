using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasSystem
    {
        public readonly int numSubdomains;

        public MenkBordasSystem(int numSubdomains)
        {
            if (numSubdomains < 1) throw new ArgumentException(
                $"There must be at least 1 subdomain, but {numSubdomains} were declared");
            this.numSubdomains = numSubdomains;
            this.Kee = new List<CSRMatrix>();
            this.Kes = new List<CSRMatrix>();
            this.Kse = new List<CSRMatrix>();
            this.B = new List<SignedBooleanMatrix>();
            this.be = new List<Vector>();
        }

        // Matrices
        public CSRMatrix Kss;
        public readonly List<CSRMatrix> Kee;
        public readonly List<CSRMatrix> Kes;
        public readonly List<CSRMatrix> Kse;
        public readonly List<SignedBooleanMatrix> B;

        // Right hand side vectors
        public Vector bs;
        public readonly List<Vector> be;

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

        public (MenkBordasMatrix matrix, Vector rhs) BuildSystem()
        { 
            // Count various dimensions
            int numDofsStd = Kss.NumRows;
            int[] subdomainStarts = new int[numSubdomains];
            int[] subdomainEnds = new int[numSubdomains];
            int nextStart = numDofsStd;
            for (int i = 0; i < numSubdomains; ++i)
            {
                subdomainStarts[i] = nextStart;
                nextStart = nextStart + Kee[i].NumRows;
                subdomainEnds[i] = nextStart;
            }
            int equationsStart = nextStart;
            int numEquations = B[0].NumRows;
            int numDofsAll = nextStart + numEquations;

            // Create the matrix
            var matrix = new MenkBordasMatrix(numSubdomains, numEquations, numDofsStd, numDofsAll, subdomainStarts, subdomainEnds,
                equationsStart, Kss, Kee.ToArray(), Kes.ToArray(), Kse.ToArray(), B.ToArray());

            // Assemble the rhs vector
            var rhs = Vector.CreateZero(numDofsAll);
            rhs.CopyFromVector(0, bs, 0, numDofsStd);
            for (int i = 0; i < numSubdomains; ++i)
            {
                rhs.CopyFromVector(subdomainStarts[i], be[i], 0, be[i].Length);
            }

            return (matrix, rhs);
        }

        public (MenkBordasPrecondMatrix matrix, MenkBordasPreconditioner precond, Vector rhs) BuildPreconditionedSystem()
        {
            // Count various dimensions
            int numDofsStd = Kss.NumRows;
            int[] subdomainStarts = new int[numSubdomains];
            int[] subdomainEnds = new int[numSubdomains];
            int nextStart = numDofsStd;
            for (int i = 0; i < numSubdomains; ++i)
            {
                subdomainStarts[i] = nextStart;
                nextStart = nextStart + Kee[i].NumRows;
                subdomainEnds[i] = nextStart;
            }
            int equationsStart = nextStart;
            int numEquations = B[0].NumRows;
            int numDofsAll = nextStart + numEquations;


            // Assemble the rhs vector
            var rhs = Vector.CreateZero(numDofsAll);
            rhs.CopyFromVector(0, bs, 0, numDofsStd);
            for (int i = 0; i < numSubdomains; ++i)
            {
                rhs.CopyFromVector(subdomainStarts[i], be[i], 0, be[i].Length);
            }

            // Create the preconditioner
            // TODO: for now just copy the DOKs, but a custom DOK is needed to avoid duplicate operations and memory.
            //      actually Kee are not needed explicitly, just their factorization. Same holds for Kss if a cholesky 
            //      preconditioner is used.
            var copyKss = DOKSymmetricColMajor.CreateFromSparseMatrix(Kss);
            var copyKee = new DOKSymmetricColMajor[numSubdomains];
            for (int i = 0; i < numSubdomains; ++i) copyKee[i] = DOKSymmetricColMajor.CreateFromSparseMatrix(Kee[i]);
            var precond = MenkBordasPreconditioner.Create(numSubdomains, numEquations, numDofsStd, numDofsAll,
                subdomainStarts, subdomainEnds, equationsStart, 
                copyKss, copyKee, B);

            // Create the matrix
            var matrix = new MenkBordasPrecondMatrix(numSubdomains, numEquations, numDofsStd, numDofsAll,
                subdomainStarts, subdomainEnds, equationsStart, 
                Kss, Kee.ToArray(), Kes.ToArray(), Kse.ToArray(), B.ToArray(), precond);

            return (matrix, precond, rhs);
        }
    }
}
