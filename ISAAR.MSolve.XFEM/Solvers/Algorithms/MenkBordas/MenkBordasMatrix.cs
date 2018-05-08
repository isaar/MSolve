using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.Algorithms.MenkBordas
{
    class MenkBordasMatrix
    {
        private readonly int numSubdomains;
        private readonly int numEquations;

        private readonly IMatrixView Kss;
        private readonly IMatrixView[] Kee;
        private readonly IMatrixView[] Kes;
        private readonly IMatrixView[] Kse;
        private readonly SignedBooleanMatrix[] B;

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
    }
}
