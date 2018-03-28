using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.LinearAlgebra;
//using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Analysis
{
    class LinearStaticAnalysisCSC: ILinearStaticAnalysis
    {
        private readonly Model2D model;
        public IVectorOLD Solution { get; private set; }

        public LinearStaticAnalysisCSC(Model2D model)
        {
            this.model = model;
        }

        public void Initialize()
        {
        }

        public void Solve()
        {
            (SymmetricDOKColMajor matrix, VectorMKL rhs) = ReduceToSimpleLinearSystem();
            using (SuiteSparseCholesky factorization = matrix.ToSymmetricCSC().FactorCholesky())
            {
                VectorMKL solution = factorization.SolveLinearSystem(rhs);
                Solution = new Vector(solution.CopyToArray()); //TODO: add method VectorMKL.ToOld() temporarily.
            }
        }

        public void PrintSolution()
        {
            Console.WriteLine("Displacements: ");
            for (int n = 0; n < model.Nodes.Count; ++n)
            {
                int xDof = model.DofEnumerator.GetFreeDofOf(model.Nodes[n], StandardDOFType.X);
                double dx = (xDof < 0) ? 0 : Solution[xDof];
                int yDof = model.DofEnumerator.GetFreeDofOf(model.Nodes[n], StandardDOFType.Y);
                double dy = (yDof < 0) ? 0 : Solution[yDof];

                Console.WriteLine("Node " + n + ": dx = " + dx + "\t\t , dy = " + dy);
            }
        }

        private (SymmetricDOKColMajor matrix, VectorMKL rhs) ReduceToSimpleLinearSystem()
        {
            /// The extended linear system is:
            /// [Kcc Kcu; Kuc Kuu] * [uc; uu] = [Fc; Fu]
            /// where c are the standard constrained dofs, f are the standard free dofs, e are the enriched dofs and 
            /// u = Union(f,c) are both the dofs with unknown left hand side vectors: uu = [uf; ue].
            /// To solve the system (for the unknowns ul):
            /// i) Kuu * uu = Fu - Kuc * uc = Feff
            /// ii) uu = Kuu \ Feff
            /// 
            SingleGlobalDOKAssembler.BuildGlobalMatrix(model, out SymmetricDOKColMajor Kuu, out Matrix2D Kuc);

            // TODO: Perhaps a dedicated class should be responsible for these vectors
            VectorMKL Fu = VectorMKL.CreateFromArray(model.CalculateFreeForces(), false); //TODO fix MKL dlls
            double[] uc = model.CalculateConstrainedDisplacements();

            // TODO: The linear algebra library is ridiculously cumbersome and limited.
            double[] KlcTimesUc = new double[Kuc.Rows];
            Kuc.Multiply(new Vector(uc), KlcTimesUc);
            VectorMKL Feff = Fu - VectorMKL.CreateFromArray(KlcTimesUc, false);

            return (Kuu, Feff);
        }
    }
}
