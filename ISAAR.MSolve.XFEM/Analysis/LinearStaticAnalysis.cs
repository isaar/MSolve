using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.LinearAlgebra;


namespace ISAAR.MSolve.XFEM.Analysis
{
    /// <summary>
    /// The linear analyzer and the solver should not compose the linear system. It should be passed as an argument to Solve.
    /// The linear system should be instantiated with the Matrix and the rhs.
    /// </summary>
    class LinearStaticAnalysis
    {
        private readonly Model2D model;
        private ISolver solver; // A solver devoid of linear system must be passed into the constructor
        public IVectorOLD Solution { get; private set; }

        public LinearStaticAnalysis(Model2D model)
        {
            this.model = model;
        }

        public void Initialize()
        {
        }

        public void Solve()
        {
            SkylineLinearSystem ls = ReduceToSimpleLinearSystem();
            solver = new SolverFBSubstitution(ls);
            solver.Initialize();
            solver.Solve();
            Solution = ls.Solution;
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

        private SkylineLinearSystem ReduceToSimpleLinearSystem()
        {
            /// The extended linear system is:
            /// [Kcc Kcu; Kuc Kuu] * [uc; uu] = [Fc; Fu]
            /// where c are the standard constrained dofs, f are the standard free dofs, e are the enriched dofs and 
            /// u = Union(f,c) are both the dofs with unknown left hand side vectors: uu = [uf; ue].
            /// To solve the system (for the unknowns ul):
            /// i) Kuu * uu = Fu - Kuc * uc = Feff
            /// ii) uu = Kuu \ Feff
            SkylineMatrix2D Kuu;
            Matrix2D Kuc;
            SingleGlobalSkylineAssembler.BuildGlobalMatrix(model, out Kuu, out Kuc);

            // TODO: Perhaps a dedicated class should be responsible for these vectors
            double[] Fu = model.CalculateFreeForces(); 
            double[] uc = model.CalculateConstrainedDisplacements();

            // TODO: The linear algebra library is ridiculously cumbersome and limited.
            double[] KlcTimesUc = new double[Kuc.Rows];
            Kuc.Multiply(new Vector(uc), KlcTimesUc);
            double[] Feff = MatrixUtilities.Substract(Fu, KlcTimesUc);

            SkylineLinearSystem ls = new SkylineLinearSystem(0, Feff);
            ls.Matrix = Kuu;
            return ls;
        }
    }
}
