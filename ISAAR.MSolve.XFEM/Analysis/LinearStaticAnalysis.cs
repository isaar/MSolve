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


namespace ISAAR.MSolve.XFEM.Analysis
{
    /// <summary>
    /// The linear analyzer and the solver should not compose the linear system. It should be passed as an argument to Solve.
    /// The linear system should be instantiated with the Matrix and the rhs. PERIOD!
    /// </summary>
    class LinearStaticAnalysis
    {
        private readonly Model2D model;
        private ISolver solver; // A solver devoid of linear system must be passed into the constructor
        private IVector solution;

        public LinearStaticAnalysis(Model2D model)
        {
            this.model = model;
        }

        public void Initialize()
        {
        }

        public void Solve()
        {
            SkylineLinearSystem ls = new SkylineLinearSystem(0, model.CalculateForces()); // Model should not be responsible for the rhs too
            ls.Matrix = SingleGlobalSkylineAssembler.BuildGlobalMatrix(model);
            solver = new SolverFBSubstitution(ls);

            solver.Initialize();
            solver.Solve();

            solution = ls.Solution;
        }

        public void PrintSolution()
        {
            Console.WriteLine("Displacements: ");
            for (int n = 0; n < model.Nodes.Count; ++n)
            {
                int xDof = model.DofEnumerator.GetStandardDofOf(model.Nodes[n], StandardDOFType.X);
                double dx = (xDof < 0) ? 0 : solution[xDof];
                int yDof = model.DofEnumerator.GetStandardDofOf(model.Nodes[n], StandardDOFType.Y);
                double dy = (yDof < 0) ? 0 : solution[yDof];

                Console.WriteLine("Node " + n + ": dx = " + dx + "\t\t , dy = " + dy);
            }
        }
    }
}
