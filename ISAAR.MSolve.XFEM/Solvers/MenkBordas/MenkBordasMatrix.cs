using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: remove MenkBordasVector, MenkBordasCG, MenkBordasMINRES
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasMatrix: ILinearTransformation<Vector>
    {
        private readonly int numSubdomains;
        private readonly int numEquations;
        private readonly int numDofsAll;
        private readonly int numDofsStd;
        private readonly int[] subdomainStarts;
        private readonly int[] subdomainEnds;
        private readonly int equationsStart;
        
        public readonly CSRMatrix Kss;
        public readonly CSRMatrix[] Kee;
        public readonly CSRMatrix[] Kes;
        public readonly CSRMatrix[] Kse;
        public readonly SignedBooleanMatrix[] B;

        public MenkBordasMatrix(int numSubdomains, int numEquations, int numDofsStd, int numDofsAll, 
            int[] subdomainStarts, int[] subdomainEnds, int equationsStart,
            CSRMatrix Kss, CSRMatrix[] Kee, CSRMatrix[] Kes, CSRMatrix[] Kse, SignedBooleanMatrix[] B)
        {
            this.numSubdomains = numSubdomains;
            this.numEquations = numEquations;
            this.numDofsStd = numDofsStd;
            this.numDofsAll = numDofsAll;
            this.subdomainStarts = subdomainStarts;
            this.subdomainEnds = subdomainEnds;
            this.equationsStart = equationsStart;

            this.Kss = Kss;
            this.Kee = Kee;
            this.Kes = Kes;
            this.Kse = Kse;
            this.B = B;
        }

        public Matrix CopyToDense()
        {
            // Dimensions
            int numStdDofs = Kss.NumRows;
            int order = numStdDofs + numEquations;
            for (int sub = 0; sub < numSubdomains; ++sub) order += Kee[sub].NumRows;
            int startBDofs = order - numEquations;
            Matrix K = Matrix.CreateZero(order, order);

            // Kss
            for (int i = 0; i < numStdDofs; ++i)
            {
                for (int j = 0; j < numStdDofs; ++j) K[i, j] = Kss[i, j];
            }

            int dofOffset = numStdDofs;
            for (int sub = 0; sub < numSubdomains; ++sub)
            {
                int numEnrDofs = Kee[sub].NumRows;

                // Kee
                for (int i = 0; i < numEnrDofs; ++i)
                {
                    for (int j = 0; j < numEnrDofs; ++j) K[dofOffset + i, dofOffset + j] = Kee[sub][i, j];
                }

                // Kes
                for (int i = 0; i < numEnrDofs; ++i)
                {
                    for (int j = 0; j < numStdDofs; ++j) K[dofOffset + i, j] = Kes[sub][i, j];
                }

                // Kse
                for (int i = 0; i < numStdDofs; ++i)
                {
                    for (int j = 0; j < numEnrDofs; ++j) K[i, dofOffset + j] = Kse[sub][i, j];
                }

                // B, Β^Τ
                for (int i = 0; i < numEquations; ++i)
                {
                    for (int j = 0; j < numEnrDofs; ++j)
                    {
                        double sign = B[sub][i, j];
                        K[startBDofs + i, dofOffset + j] = sign;
                        K[dofOffset + j, startBDofs + i] = sign;
                    }
                }

                dofOffset += numEnrDofs;
            }

            return K;
        }

        public Vector Multiply(Vector x)
        {
            var y = Vector.CreateZero(numDofsAll);
            var xs = x.Slice(0, numDofsStd);
            var xc = x.Slice(equationsStart, numDofsAll);
            Vector ys = Kss.MultiplyRight(xs); // Rows correspond to global standard dofs
            var yc = Vector.CreateZero(numEquations); // Rows correspond to the continuity equations. TODO: try to avoid this
            for (int i = 0; i < numSubdomains; ++i)
            {
                var xe = x.Slice(subdomainStarts[i], subdomainEnds[i]);
                ys.AddIntoThis(Kse[i].MultiplyRight(xe)); 
                Vector ye = Kes[i].MultiplyRight(xs);
                ye.AddIntoThis(Kee[i].MultiplyRight(xe));
                ye.AddIntoThis(B[i].MultiplyRight(xc, true));
                yc.AddIntoThis(B[i].MultiplyRight(xe, false));
                y.CopyFromVector(subdomainStarts[i], ye, 0, ye.Length);
            }
            y.CopyFromVector(0, ys, 0, numDofsStd);
            y.CopyFromVector(equationsStart, yc, 0, numEquations);
            return y;
        }

        public MenkBordasVector MultiplyRight(MenkBordasVector vector)
        {
            Vector ys = Kss.MultiplyRight(vector.Vs); // Rows correspond to global standard dofs
            var ye = new Vector[numSubdomains]; // Rows correspond to subdomain enriched dofs
            var yc = Vector.CreateZero(numEquations); // Rows correspond to the continuity equations. TODO: try to avoid this
            for (int i = 0; i < numSubdomains; ++i)
            {
                ys.AddIntoThis(Kse[i].MultiplyRight(vector.Ve[i]));
                ye[i] = Kes[i].MultiplyRight(vector.Vs);
                ye[i].AddIntoThis(Kee[i].MultiplyRight(vector.Ve[i]));
                ye[i].AddIntoThis(B[i].MultiplyRight(vector.Vc, true)); //TODO: verify that it is not needed
                yc.AddIntoThis(B[i].MultiplyRight(vector.Ve[i], false)); // Pretty sure this will not be 0
            }
            return new MenkBordasVector(numSubdomains, numEquations, ys, ye, yc);
        }

        public void WriteToConsole()
        {
            FullMatrixWriter.NumericFormat = new GeneralNumericFormat();

            Console.WriteLine("Kss: ");
            (new FullMatrixWriter(Kss)).WriteToConsole();
            Console.WriteLine();
            for (int i = 0; i < numSubdomains; ++i)
            {
                Console.WriteLine($"Subdomain {i} - Kee: ");
                (new FullMatrixWriter(Kee[i])).WriteToConsole();
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - Kes: ");
                (new FullMatrixWriter(Kes[i])).WriteToConsole();
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - Kse: ");
                (new FullMatrixWriter(Kse[i])).WriteToConsole();
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - B: ");
                (new FullMatrixWriter(B[i])).WriteToConsole();
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - transpose B: ");
                (new FullMatrixWriter(B[i].Transpose())).WriteToConsole();
                Console.WriteLine();
            }
        }

        public void WriteToFiles(string directoryPath)
        {
            var writer = new MatlabWriter();
            writer.WriteSparseMatrix(Kss, directoryPath + "Kss.txt");
            for (int i = 0; i < numSubdomains; ++i)
            {
                writer.WriteSparseMatrix(Kee[i], directoryPath + $"Kee{i + 1}.txt");
                writer.WriteSparseMatrix(Kes[i], directoryPath + $"Kes{i + 1}.txt");
                writer.WriteSparseMatrix(B[i], directoryPath + $"B{i + 1}.txt");
            }
        }
    }
}
