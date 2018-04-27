using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class SkylineSolverOLD: ISolver
    {
        private readonly Model2D model;

        public SkylineSolverOLD(Model2D model)
        {
            this.model = model;
            this.Logger = new SolverLogger();
        }

        public IDOFEnumerator DOFEnumerator { get; private set; }

        public SolverLogger Logger { get; }

        public Vector Solution { get; private set; }

        public void Initialize()
        {
            Logger.InitializationTime = 0;
        }

        public void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            DOFEnumerator = DOFEnumeratorSeparate.Create(model); // Actually this enumeration will result in huge bandwidths.
            //DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            SkylineLinearSystem ls = ReduceToSimpleLinearSystem();
            var solver = new SolverFBSubstitution(ls); // A solver devoid of linear system must be passed into the constructor
            solver.Initialize();
            solver.Solve();
            double[] solutionArray = new double[ls.Solution.Length];
            ls.Solution.CopyTo(solutionArray, 0);
            Solution = Vector.CreateFromArray(solutionArray);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }

        /// <summary>
        /// The extended linear system is:
        /// [Kcc Kcu; Kuc Kuu] * [uc; uu] = [Fc; Fu]
        /// where c are the standard constrained dofs, f are the standard free dofs, e are the enriched dofs and 
        /// u = Union(f,c) are both the dofs with unknown left hand side vectors: uu = [uf; ue].
        /// To solve the system (for the unknowns ul):
        /// i) Kuu * uu = Fu - Kuc * uc = Feff
        /// ii) uu = Kuu \ Feff 
        /// </summary>
        /// <returns></returns>
        private SkylineLinearSystem ReduceToSimpleLinearSystem()
        {
            (Numerical.LinearAlgebra.SkylineMatrix2D Kuu, Matrix Kuc) = SingleGlobalSkylineAssemblerOLD.BuildGlobalMatrix(model, DOFEnumerator);

            // TODO: Perhaps a dedicated class should be responsible for these vectors
            Vector Fu = model.CalculateFreeForces(DOFEnumerator);
            Vector uc = model.CalculateConstrainedDisplacements(DOFEnumerator);
            Vector Feff = Fu - Kuc * uc;

            #region DEBUG code
            //Console.WriteLine("uc:");
            //(new FullVectorWriter(uc)).WriteToConsole();
            ////Console.WriteLine("Kuu:");
            //Console.WriteLine("Kcu:");
            //(new FullMatrixWriter(Kuc)).WriteToConsole();
            //Console.WriteLine("Feff:");
            //(new FullVectorWriter(Feff)).WriteToConsole();
            #endregion

            var ls = new SkylineLinearSystem(0, Feff.CopyToArray());
            ls.Matrix = Kuu;
            return ls;
        }
    }
}
