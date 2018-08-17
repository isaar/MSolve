using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: remove MenkBordasVector, MenkBordasCG, MenkBordasMINRES
namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasMatrix: ILinearTransformation<Vector>
    {
        public readonly MenkBordasSystem.Dimensions dim;
        public readonly CsrMatrix Kss;
        public readonly CsrMatrix[] Kee;
        public readonly CsrMatrix[] Kes;
        public readonly CsrMatrix[] Kse;
        public readonly SignedBooleanMatrix[] B;

        public MenkBordasMatrix(MenkBordasSystem.Dimensions dim,
            CsrMatrix Kss, CsrMatrix[] Kee, CsrMatrix[] Kes, CsrMatrix[] Kse, SignedBooleanMatrix[] B)
        {
            this.dim = dim;
            this.Kss = Kss;
            this.Kee = Kee;
            this.Kes = Kes;
            this.Kse = Kse;
            this.B = B;
        }

        public Matrix CopyToDense()
        {
            // Dimensions
            int order = dim.NumDofsStd + dim.NumEquations;
            for (int sub = 0; sub < dim.NumSubdomains; ++sub) order += Kee[sub].NumRows;
            int startBDofs = order - dim.NumEquations;
            Matrix K = Matrix.CreateZero(order, order);

            // Kss
            for (int i = 0; i < dim.NumDofsStd; ++i)
            {
                for (int j = 0; j < dim.NumDofsStd; ++j) K[i, j] = Kss[i, j];
            }

            int dofOffset = dim.NumDofsStd;
            for (int sub = 0; sub < dim.NumSubdomains; ++sub)
            {
                // Kee
                for (int i = 0; i < dim.NumDofsEnr; ++i)
                {
                    for (int j = 0; j < dim.NumDofsEnr; ++j) K[dofOffset + i, dofOffset + j] = Kee[sub][i, j];
                }

                // Kes
                for (int i = 0; i < dim.NumDofsEnr; ++i)
                {
                    for (int j = 0; j < dim.NumDofsStd; ++j) K[dofOffset + i, j] = Kes[sub][i, j];
                }

                // Kse
                for (int i = 0; i < dim.NumDofsStd; ++i)
                {
                    for (int j = 0; j < dim.NumDofsEnr; ++j) K[i, dofOffset + j] = Kse[sub][i, j];
                }

                // B, Β^Τ
                for (int i = 0; i < dim.NumEquations; ++i)
                {
                    for (int j = 0; j < dim.NumDofsEnr; ++j)
                    {
                        double sign = B[sub][i, j];
                        K[startBDofs + i, dofOffset + j] = sign;
                        K[dofOffset + j, startBDofs + i] = sign;
                    }
                }

                dofOffset += dim.NumDofsEnr;
            }

            return K;
        }

        public Vector Multiply(Vector x)
        {
            var y = Vector.CreateZero(dim.NumDofsAll);
            var xs = x.GetSubvector(0, dim.NumDofsStd);
            var xc = x.GetSubvector(dim.EquationsStart, dim.NumDofsAll);
            Vector ys = Kss.MultiplyRight(xs); // Rows correspond to global standard dofs
            var yc = Vector.CreateZero(dim.NumEquations); // Rows correspond to the continuity equations. TODO: try to avoid this

            foreach (var sub in dim.Subdomains)
            {
                int i = sub.ID;
                var xe = x.GetSubvector(dim.SubdomainStarts[sub], dim.SubdomainEnds[sub]);
                ys.AddIntoThis(Kse[i].MultiplyRight(xe)); 
                Vector ye = Kes[i].MultiplyRight(xs);
                ye.AddIntoThis(Kee[i].MultiplyRight(xe));
                ye.AddIntoThis(B[i].MultiplyRight(xc, true));
                yc.AddIntoThis(B[i].MultiplyRight(xe, false));
                y.SetSubvector(ye, dim.SubdomainStarts[sub]);
            }
            y.SetSubvector(ys, 0);
            y.SetSubvector(yc, dim.EquationsStart);
            return y;
        }

        public MenkBordasVector MultiplyRight(MenkBordasVector vector)
        {
            Vector ys = Kss.MultiplyRight(vector.Vs); // Rows correspond to global standard dofs
            var ye = new Vector[dim.NumSubdomains]; // Rows correspond to subdomain enriched dofs
            var yc = Vector.CreateZero(dim.NumEquations); // Rows correspond to the continuity equations. TODO: try to avoid this
            for (int i = 0; i < dim.NumSubdomains; ++i)
            {
                ys.AddIntoThis(Kse[i].MultiplyRight(vector.Ve[i]));
                ye[i] = Kes[i].MultiplyRight(vector.Vs);
                ye[i].AddIntoThis(Kee[i].MultiplyRight(vector.Ve[i]));
                ye[i].AddIntoThis(B[i].MultiplyRight(vector.Vc, true)); //TODO: verify that it is not needed
                yc.AddIntoThis(B[i].MultiplyRight(vector.Ve[i], false)); // Pretty sure this will not be 0
            }
            return new MenkBordasVector(dim.NumSubdomains, dim.NumEquations, ys, ye, yc);
        }

        public void WriteToConsole()
        {
            var writer = new FullMatrixWriter();
            writer.NumericFormat = new GeneralNumericFormat();

            Console.WriteLine("Kss: ");
            writer.WriteToConsole(Kss);
            Console.WriteLine();
            for (int i = 0; i < dim.NumSubdomains; ++i)
            {
                Console.WriteLine($"Subdomain {i} - Kee: ");
                writer.WriteToConsole(Kee[i]);
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - Kes: ");
                writer.WriteToConsole(Kes[i]);
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - Kse: ");
                writer.WriteToConsole(Kse[i]);
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - B: ");
                writer.WriteToConsole(B[i]);
                Console.WriteLine();
                Console.WriteLine($"Subdomain {i} - transpose B: ");
                writer.WriteToConsole(B[i].Transpose());
                Console.WriteLine();
            }
        }

        public void WriteToFiles(string directoryPath)
        {
            var writer = new MatlabWriter();
            writer.WriteToFile(Kss, directoryPath + "Kss.txt");
            for (int i = 0; i < dim.NumSubdomains; ++i)
            {
                writer.WriteToFile(Kee[i], directoryPath + $"Kee{i + 1}.txt");
                writer.WriteToFile(Kes[i], directoryPath + $"Kes{i + 1}.txt");
                writer.WriteToFile(B[i], directoryPath + $"B{i + 1}.txt");
            }
        }
    }
}
