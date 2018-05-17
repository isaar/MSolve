using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasMatrix
    {
        private readonly int numSubdomains;
        private readonly int numEquations;

        public readonly IMatrixView Kss;
        public readonly IMatrixView[] Kee;
        public readonly IMatrixView[] Kes;
        public readonly IMatrixView[] Kse;
        public readonly SignedBooleanMatrix[] B;

        public MenkBordasMatrix(int numSubdomains, int numEquations,
            IMatrixView Kss, IMatrixView[] Kee, IMatrixView[] Kes, IMatrixView[] Kse, SignedBooleanMatrix[] B)
        {
            this.numSubdomains = numSubdomains;
            this.numEquations = numEquations;
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
    }
}
